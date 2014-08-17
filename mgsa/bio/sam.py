
import datetime
import random
import re
import sys

import bio
import fasta

class SamToVCF(object):
  '''
    similar to SamToFasta but use a reference to only get variants
  '''
  def __init__( self, sam, reference, target_vcf, log, sam_sorted=False ):
    '''
      @sam: file handle to sam file
      @reference: file handle to fasta file
      @target: VCF object
    '''
    sam_to_fasta = SamToFasta( sam, log ) # builds fasta
    candidate = sam_to_fasta.fasta
    reference = fasta.Fasta( fasta.FastaReader( reference ) )
    pos = 0
    while pos < candidate.length:
      actual_base = reference.base_at( pos )
      if actual_base is None:
        break # reference is done
      move, candidate_base = candidate.consensus_at( pos )
      if candidate_base != actual_base and candidate_base != 'N':
        target_vcf.snp( pos, actual_base, candidate_base ) # mutation
      pos += 1
    # check for sequence finishing early
    log( "pos %i candidate length %i reference length %i" % ( pos, candidate.length, reference.length ) )
    if pos < candidate.length: # reference finished early
      target_vcf.indel( pos=pos-1, before=reference.base_at( pos-1 ), after=candidate.consensus( start=pos-1 ) )
    elif reference.length > candidate.length: # candidate finished early
      target_vcf.indel( pos=pos-1, before=reference.base_at( pos-1 ), after=reference.fragment( start=pos-1 ) )

class SamToFasta(object):
  '''
    read and evaluate assembler data
    takes a sam file and builds an imputed genome
  '''  

  def __init__( self, sam, log, allow_indels=True ):
    self.allow_indels = allow_indels
    self.log = log
    self.fasta = fasta.ProbabilisticFasta(log)
    self.stats = { 'mapped': 0, 'unmapped': 0, 'unknown_mapping': 0, 'lines': 0, 'inversions': 0 }
    for line in sam:
      self.parse_line( line.strip() )
      self.stats['lines'] += 1
      if self.stats['lines'] < 5000 and self.stats['lines'] % 1000 == 0 or self.stats['lines'] % 10000 == 0:
        self.log( 'SamToFasta: %i lines processed' % self.stats['lines'] )
    self.log( self.stats )

  def write( self, fh ):
    self.log( 'length is %i' % self.fasta.length )

    fh.write( '>reassembled\n')
    pos = 1
    while pos < self.fasta.length:
      fh.write( '%s\n' % self.fasta.consensus( pos, pos + 70 ) )
      pos += 70

  def parse_line( self, line ):
    #0 SimSeq_1101 => qname
    #1 163 => flag
    #2 gi|205355313|ref|NZ_ABGQ01000001.1| => rname
    #3 1 => pos
    #4 44 => mapq
    #5 100M => cigar
    #6 =
    #7 227
    #8 326
    #9 AAATGTCTTCATTACTTACTTTATATATAAATCCTATGTTTATTTTTATTGTTGTTTTAGATGATACTTAGAATCGTGTTTAAAAAAAAGTTTCCTGCTG
    #10 BGG=77FGBGG?GDGB7GED@CEGECDA?EG3D8:D.??06GB?-GCDGCA9G#AG?=AAEFBFFG=@DAAD#EBD;EC5#GD#5DE##A#BF#B#?=##
    #AS:i:192 #XN:i:0 #XM:i:2 #XO:i:0 #XG:i:0 #NM:i:2 #MD:Z:87G5G6 #YS:i:200 #YT:Z:CP

    if line.startswith( '@' ): # skip header
      pass
    else:
      fields = line.split()
      flag = int(fields[1])
      #if self.stats['lines'] < 100:
      #  log( 'flag is %i' % flag )
      if flag & 0x10 != 0:
        self.stats['inversions'] += 1
      if flag & 0x04 != 0:
        self.stats['unmapped'] += 1
      else:
        if flag & 0x02 != 0:
          self.stats['mapped'] += 1
        else:
          self.stats['unknown_mapping'] += 1 # still try to map
        pos = int(fields[3]) - 1 # pos in genome; 0 based
        fragment_pos = 0
        genome_pos = 0
        cigar = fields[5]
        for cigar_match in re.findall( '([0-9]+)([MIDNSHP=X])', cigar ):
          cigar_len = int(cigar_match[0])

          if cigar_match[1] == 'M': # match
            self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + genome_pos )
            genome_pos += cigar_len
            fragment_pos += cigar_len

          if self.allow_indels and cigar_match[1] == 'I': # insertion to the reference
            self.fasta.insert( fields[9][fragment_pos:fragment_pos+cigar_len], pos + genome_pos )
            fragment_pos += cigar_len

          if self.allow_indels and cigar_match[1] == 'D': # deletion from the reference
            self.fasta.delete( pos + genome_pos, cigar_len )
            genome_pos += cigar_len

          if cigar_match[1] == 'N': # skipped region from the reference
            genome_pos += cigar_len

          if cigar_match[1] == 'S': # soft clipping
            #self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + fragment_pos, confidence=0.5 )
            self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + fragment_pos - cigar_len, confidence=0.5 )
            #genome_pos += cigar_len
            fragment_pos += cigar_len

          if cigar_match[1] == 'H': # hard clipping
            fragment_pos += cigar_len

          if cigar_match[1] == 'P': # padding
            genome_pos += cigar_len

          if cigar_match[1] == '=': # sequence match
            self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + genome_pos )
            genome_pos += cigar_len
            fragment_pos += cigar_len

          if cigar_match[1] == 'X': # sequence mismatch
            genome_pos += cigar_len
            fragment_pos += cigar_len

          key = 'cigar_%s' % cigar_match[1]
          if key not in self.stats:
            self.stats[key] = 0
          self.stats[key] += cigar_len

class SamAccuracyEvaluator(object):
  '''
    evaluate sam accuracy given the correct sequence marker
  '''  
  def __init__( self, sam, variation_map=None, log=bio.log_stderr, verbose=False ):
    self.variation_map = self.parse_variation_map( variation_map )
    self.log = log
    self.verbose = verbose
    self.stats = { 'mapped': 0, 'unmapped': 0, 'unknown_mapping': 0, 'lines': 0, 'correct': 0, 'incorrect': 0, 'soft_clipping': 0, 'hard_clipping': 0, 'matched': 0, 'correct_mapq': 0.0, 'incorrect_mapq': 0.0 }
    self.incorrect = []
    self.incorrect_diff = {}
    for line in sam:
      self.parse_line( line.strip() )
      self.stats['lines'] += 1
      if self.stats['lines'] < 5000 and self.stats['lines'] % 1000 == 0 or self.stats['lines'] % 10000 == 0:
        self.log( 'SamAccuracyEvaluator: %i lines processed' % self.stats['lines'] )
    if self.stats['correct'] > 0:
      self.stats['correct_mapq'] /= self.stats['correct']
    if self.stats['incorrect'] > 0:
      self.stats['incorrect_mapq'] /= self.stats['incorrect']

  def parse_variation_map( self, vmap ):
    if vmap is None:
      return None
    result = {}
    for line in vmap:
      key, value = line.strip().split()
    return result

  def parse_line( self, line ):
    #0 SimSeq_1101 => qname
    #1 163 => flag
    #2 gi|205355313|ref|NZ_ABGQ01000001.1| => rname
    #3 1 => pos
    #4 44 => mapq
    #5 100M => cigar
    #6 =
    #7 227
    #8 326
    #9 AAATGTCTTCATTACTTACTTTATATATAAATCCTATGTTTATTTTTATTGTTGTTTTAGATGATACTTAGAATCGTGTTTAAAAAAAAGTTTCCTGCTG
    #10 BGG=77FGBGG?GDGB7GED@CEGECDA?EG3D8:D.??06GB?-GCDGCA9G#AG?=AAEFBFFG=@DAAD#EBD;EC5#GD#5DE##A#BF#B#?=##
    #AS:i:192 #XN:i:0 #XM:i:2 #XO:i:0 #XG:i:0 #NM:i:2 #MD:Z:87G5G6 #YS:i:200 #YT:Z:CP
    if line.startswith( '@' ): # skip header
      pass
    else:
      fields = line.split()
      flag = int(fields[1])
      #if self.stats['lines'] < 100:
      #  log( 'flag is %i' % flag )
      if flag & 0x04 != 0: # unmapped
        self.stats['unmapped'] += 1
      else:
        if flag & 0x02 != 0: # mapped
          self.stats['mapped'] += 1
        else:
          self.stats['unknown_mapping'] += 1 # but still try to map
        pos = int(fields[3]) - 1 # pos in genome; 0 based
        cigar = fields[5]
        prematch = True # clipping at start
        correct_pos = int(re.sub( 'mgsa_seq_([0-9]*).*', '\\1', fields[0] ))
        if 'variation_' in fields[0]:
          variations = re.sub( '.*variation_([SID0-9,-]*).*', '\\1', fields[0] ).split(',')
          #print "field %s -> %s" % ( fields[0], variations )
        elif self.variation_map is not None and fields[0] in self.variation_map:
          variations = re.sub( '([SID0-9,-]*).*', '\\1', self.variation_map[fields[0]] ).split(',')
        else:
          variations = ()
        for cigar_match in re.findall( '([0-9]+)([MIDNSHP=X])', cigar ): # match description
          cigar_len = int(cigar_match[0])
          if cigar_match[1] == 'M': # match
            if prematch:
              indel_offset = 0
              for variation in variations:
                variation_pos, variation_length = [ int(v) for v in variation[1:].split( '-' ) ]
                if variation_pos <= 0:
                  if variation[0] == 'I':
                    indel_offset -= variation_length
                  if variation[0] == 'D':
                    indel_offset += variation_length
              pos -= indel_offset

            prematch = False # no pre clipping
            self.stats['matched'] += cigar_len

          if cigar_match[1] == 'S': # soft clipping
            self.stats['soft_clipping'] += cigar_len
            if prematch:
              indel_offset = 0
              for variation in variations:
                variation_pos, variation_length = [ int(v) for v in variation[1:].split( '-' ) ]
                # migrate insertion to start of pattern
                while variation_pos > 0 and fields[9][variation_pos-1] == fields[9][variation_pos]:
                  variation_pos -= 1
                  #print "variation %s at read %i migrated from %s to %i slen %i" % ( variation, correct_pos, variation[1:], variation_pos, cigar_len )
                if variation_pos <= cigar_len: # variant in clipped region TODO variant could be partially in clipped area
                  #print "included", variation
                  if variation[0] == 'I':
                    indel_offset -= variation_length
                  if variation[0] == 'D':
                    indel_offset += variation_length
                else:
                  pass #print "ignored", variation
              pos -= cigar_len + indel_offset
              #print "newpos %i cigarlen %i offset %i" % ( pos, cigar_len, indel_offset )

          if cigar_match[1] == 'H': # hard clipping
            self.stats['hard_clipping'] += cigar_len
            if prematch:
              indel_offset = 0
              for variation in variations:
                variation_pos, variation_length = [ int(v) for v in variation[1:].split( '-' ) ]
                if variation_pos <= cigar_len:
                  if variation[0] == 'I':
                    indel_offset -= variation_length
                  if variation[0] == 'D':
                    indel_offset += variation_length
              pos -= cigar_len + indel_offset

        #print "variations", variations
        if correct_pos == pos:
          self.stats['correct'] += 1
          self.stats['correct_mapq'] += int(fields[4])
          for c in variations:
            if 'correct_%s' % c[0] not in self.stats:
              self.stats['correct_%s' % c[0] ] = 0
            self.stats['correct_%s' % c[0] ] += 1
        else:
          self.stats['incorrect'] += 1
          self.stats['incorrect_mapq'] += int(fields[4])
          for c in variations:
            if 'incorrect_%s' % c[0] not in self.stats:
              self.stats['incorrect_%s' % c[0] ] = 0
            self.stats['incorrect_%s' % c[0] ] += 1
          if self.verbose:
            self.incorrect.append( { 
              'correct_pos': correct_pos, 
              'provided_pos': pos, 
              'mapq': fields[4],
              'cigar': cigar,
              'label': fields[0],
              'read': fields[9][:10]
            } )
          else:
            self.incorrect.append( { 'correct_pos': correct_pos, 'provided_pos': pos, 'mapq': fields[4] } )
          diff = pos - correct_pos
          if diff not in self.incorrect_diff:
            self.incorrect_diff[ diff ] = 0
          self.incorrect_diff[ diff ] += 1
 
class SamWriter(object):
  def __init__( self, fh ):
    self._fh = fh

  def write_header( self, fasta ):
    self._fh.write( '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:generated\tLN:%i\n@PG\tID:baseline\tPN:baseline\tVN:1.0.0\n' % len(fasta) )

  def write_no_match( self, name, dna, extra, confidence ):
    # mge  4   *            0   0   *   *   0   0  CTGAAG ~~~ YT:Z:UU
    self._fh.write( '%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n' % ( name[1:], 4, '*', 0, 0, '*', '*', 0, 0, dna, confidence, '' ) )

  def write_match( self, name, dna, extra, confidence, match ):
    '''
      match is 0-based
    '''
    # mge  0   generated 3177 42 40M3D10M * 0 0 CGC ~~~ AS:i:86 XN:i:0  XM:i:0  XO:i:1  XG:i:3  NM:i:3  MD:Z:40^CAA10   YT:Z:UU
    self._fh.write( '%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n' % ( name[1:], 2, 'generated', match + 1, 44, '%iM' % len(dna), '*', 0, 0, dna, confidence, '' ) )

  def write_multi_match( self, name, dna, extra, confidence, matches ):
    # choose a random match
    match = matches[ random.randint( 0, len(matches) - 1 ) ]
    self._fh.write( '%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n' % ( name[1:], 0, 'generated', match + 1, 1, '%iM' % len(dna), '*', 0, 0, dna, confidence, '' ) )

