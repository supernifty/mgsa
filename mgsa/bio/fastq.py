
import sys

import bio
import fasta

class FastqGenerator(object):
  '''
    generates a complete fastq with various mutation options from a fasta
  '''
  def __init__( self, fasta_fh, cfg, vcf, variation_map=None, target_fh=sys.stdout, log=bio.log_stderr ):
    probabilistic = False

    read_probability = 1. * cfg['coverage'] / cfg['read_length']
    read_every = 1. / read_probability
    inversion_prob = cfg['inversion_prob']
    if inversion_prob > 0:
      inversion_every = 1. / inversion_prob
    else:
      inversion_every = 1e12

    error_generator = fasta.ErrorGenerator( fasta.ErrorGenerator.create_uniform_error_profile( cfg['error_prob'] ) )

    # read candidate fasta
    dna = ''
    base_candidate_position = 0
    time_to_next = read_every
    time_to_inversion = inversion_every
    max_pos = 0
    lines = 0
    #if vcf is not None:
    #  bio.log_stderr( "snp: %s; indel: %s" % ( vcf.snp_list, vcf.indel_list ) )
    for line in fasta_fh:
      if lines < 10 or lines % 1000 == 0:
        log( 'generate_reads: %i lines processed' % lines )
      lines += 1

      if line.startswith( '>' ):
        continue

      line = line.strip()
      dna += line

      process = len(dna) - cfg['read_length'] + 1
      if process > 0:
        for i in xrange(0, process):
          time_to_next -= 1
          if probabilistic and random.random() < read_probability or \
             not probabilistic and time_to_next <= 0:
            time_to_inversion -= 1
            candidate_position = base_candidate_position + i # reference == candidate
            if vcf is None: # no variations
              reference_position = candidate_position # reference == candidate
              variations = set()
              offset = 0
              length = 0
            else: # variations
              reference_position, offset, length = vcf.candidate_position_to_reference_position( candidate_position )
              variations = vcf.variations( reference_position, reference_position + cfg['read_length'], offset=offset )

            dna_with_errors = error_generator.apply_errors( dna[i:i+cfg['read_length']] ) # add errors to read
    
            inversion = probabilistic and random.random() < inversion_prob or not probabilistic and time_to_inversion <= 0
            if inversion:
              time_to_inversion += inversion_every
    
            self._write( reference_position, dna_with_errors, '~' * cfg['read_length'], variation_map=variation_map, target_fh=target_fh, variations=variations, offset=offset, debug='', inversion=inversion, length=length)
    
            max_pos = max(max_pos, reference_position + cfg['read_length'])
            time_to_next += read_every
        base_candidate_position += process
        dna = dna[process:]

  def _write( self, pos, sequence, quality, variation_map, target_fh, variations=set(), offset=0, length=0, debug='', inversion=False ):
    '''
      write a sequence of given quality
      debug is not used and for humans only
    '''
    if inversion:
      inversion_text = 'inversion_'
      sequence = bio.reverse_complement( sequence )
      quality = quality[::-1]
    else:
      inversion_text = ''
  
    target_fh.write( '@mgsa_seq_%i~%i~%i\n' % ( pos, offset, length ) ) # sam is 0 indexed
    if len(variations) > 0:
      variation_map.write( '@mgsa_seq_%i~%i~%i: %s_%s%s\n' % ( pos, offset, length, ','.join([ v for v in variations ]), inversion_text, debug ) ) # sam is 0 indexed

    target_fh.write( sequence )
    target_fh.write( '\n+\n' )
    target_fh.write( quality ) # quality 
    target_fh.write( '\n' )
  
   
class FastqPosGenerator(object):
  '''
    generate reads around a specific base
  '''

  def __init__( self, fasta_fh, log=bio.log_stderr ):
    self.fasta = fasta_fh
    self.log = log
    self.sequence = ''
    self.finished_read = False
    self._update_sequence()

  def _update_sequence( self, required_length=-1 ):
    if self.finished_read:
      return
    # for the moment pull in the whole sequence
    for line in self.fasta:
      if line.startswith( '>' ):
        continue
      self.sequence += line.strip()
    self.finished_read = True
    
  def __len__(self):
    return len(self.sequence)

  def _variation_rule( self, ref ):
    '''
      default snp is A, else T
    '''
    if ref == 'A':
      return 'T'
    else:
      return 'A'

  def write( self, target_fh=sys.stdout, pos=0, read_length=100, variation=None ):
    '''
      fastq file to write to
      @target_fh: fastq file
      @pos: pos of fasta to cover
    ''' 
    self._update_sequence( pos + read_length )
    start_pos = max(0, pos - read_length + 1)
    subsequence = self.sequence[start_pos : pos + read_length]
    if variation == 'snp':
      snp_pos = pos - start_pos
      ref = subsequence[snp_pos]
      alt = self._variation_rule( ref )
      subsequence = ''.join( ( subsequence[:snp_pos], alt, subsequence[snp_pos+1:] ) )
    written = 0
    #self.log('subsequence: %s' % subsequence )
    for offset in xrange(0, len(subsequence) - read_length + 1):
      target_fh.write( '@mgsa_seq_%i~%i~0\n' % ( start_pos, offset ) ) # sam is 0 indexed
      read = subsequence[offset:offset+read_length]
      quality = '~' * len(read)
      target_fh.write( read )
      target_fh.write( '\n+\n' )
      target_fh.write( quality ) # quality 
      target_fh.write( '\n' )
      written += 1
    return written   
     
