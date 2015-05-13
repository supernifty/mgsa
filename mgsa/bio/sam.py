
import collections
import datetime
import random
import re
import struct
import subprocess
import sys
import zlib

import bio
import fasta

SOFT_CLIP_CONFIDENCE = 0.0

class SamToMultiChromosomeVCF(object):
  def __init__( self, sam, multi_fasta_reference, target_vcf, log, sam_sorted=False, call_strategy='consensus' ):
    '''
      @sam: file like sam file with multiple chromosomes
      @reference: file like fasta reference with multiple fastas
      @target_vcf: where to write the variations
    '''
    self.sam_to_multi_chromosome_fasta = SamToMultiChromosomeFasta( sam, log ) # builds sam_to_fasta
    self.multi_fasta_reference = fasta.MultiFastaReaderContainer( multi_fasta_reference )
    self.sam_to_vcf = {}
    for chromosome in self.sam_to_multi_chromosome_fasta.fastas: # sam_to_fasta for each chromosome 
      target_vcf.chromosome = chromosome
      self.sam_to_vcf[chromosome] = SamToVCF( self.sam_to_multi_chromosome_fasta.fastas[chromosome], self.multi_fasta_reference.find_chromosome( chromosome ), target_vcf, log, sam_sorted, call_strategy ) # sam_to_vcf for each fasta
  

class SamToVCF(object):
  '''
    similar to SamToFasta but use a reference to only get variants
    deals with a single chromosome
  '''

  @staticmethod
  def instance( sam, reference, target_vcf, log, sam_sorted=False, call_strategy='consensus' ):
    '''
      @sam: file handle to sam file
      @reference: file handle to fasta file
      @target: VCF object
      @call_strategy: consensus|conservative|aggressive TODO not implemented
    '''
    return SamToVCF( SamToFasta( sam, log ), fasta.FastaReader( reference ), target_vcf, log, sam_sorted, call_strategy )

  def __init__( self, sam_to_fasta, fasta_reader_reference, target_vcf, log, sam_sorted=False, call_strategy='consensus' ):
    '''
      @sam_to_fasta: sam to fasta object
      @fasta_reader_reference: fasta_reader
      @target: VCF object
      @call_strategy: consensus|conservative|aggressive TODO not implemented
    '''
    self.sam_to_fasta = sam_to_fasta
    candidate = self.sam_to_fasta.fasta # ProbabilisticFasta
    reference = fasta.Fasta( fasta_reader_reference )
    candidate_pos = 0
    delete_start = None
    # TODO ploidy; make heterozygous calls and pass on to target vcf
    while candidate_pos < candidate.length:
      actual_base = reference.base_at( candidate_pos )
      if actual_base is None:
        break # reference is done
      candidate_move, candidate_variation, evidence, coverage = candidate.consensus_at( candidate_pos )
      #log( 'SamToVCF: candidate_pos %i move %i base %s evidence %s coverage %f actual_base %s' % ( candidate_pos, candidate_move, candidate_variation, evidence, coverage, actual_base ) ) # testing

      if candidate_move > 0 and delete_start: # the end of a delete
         #target_vcf.indel( pos=delete_start - 1, before='%s%s' % ( reference.base_at( delete_start-1 ), delete_variation ), after='%s' % ( reference.base_at( delete_start-1 ) ), coverage=coverage )
         #log( 'SamToVCF: deletion at %i before %s/%s after %s coverage %f' % ( delete_start - 1, reference.base_at( delete_start-1 ), delete_variation, reference.base_at( delete_start-1 ), coverage ) )
         target_vcf.indel( pos=delete_start - 1, before='%s%s' % ( reference.base_at( delete_start-1 ), delete_variation ), after='%s' % ( reference.base_at( delete_start-1 ) ), coverage=coverage )
         log( 'SamToVCF: deletion at %i before %s/%s after %s coverage %f' % ( delete_start - 1, reference.base_at( delete_start-1 ), delete_variation, reference.base_at( delete_start-1 ), coverage ) )
         delete_start = None
      
      if candidate_move > 1: # insertion
        # note that the current base is included in candidate_variation
        target_vcf.indel( pos=candidate_pos, before='%s%s' % ( reference.base_at( candidate_pos-1 ), reference.base_at( candidate_pos ) ), after='%s%s' % ( reference.base_at( candidate_pos-1 ), candidate_variation ), coverage=coverage )
        log( 'SamToVCF: insertion %i before %s/%s after %s/%s coverage %f' % ( candidate_pos, reference.base_at( candidate_pos-1 ), reference.base_at( candidate_pos ), reference.base_at( candidate_pos-1 ), reference.base_at( candidate_pos ), coverage ) )
      elif candidate_move < 1: # start of deletion
        if delete_start is None:
          delete_start = candidate_pos
          delete_variation = reference.base_at(candidate_pos)
        else:
          delete_variation += reference.base_at(candidate_pos)
      else: # candidate_move is 1 => snp
        if candidate_variation != actual_base and candidate_variation != 'N':
          confidence = 1. * evidence / coverage
          target_vcf.snp( pos=candidate_pos, ref=actual_base, alt=candidate_variation, coverage=coverage, confidence=confidence ) # mutation
          log( 'SamToVCF: snp at %i with confidence %f due to %s' % ( candidate_pos, confidence, candidate.count(candidate_pos) ) )
      candidate_pos += 1 # next base

    # check for sequence finishing early
    #log( "pos %i candidate length %i reference length %i" % ( pos, candidate.length, reference.length ) )
    if candidate_pos < candidate.length: # reference finished early
      #print "pos %i rl %i cl %i" % ( pos, reference.length, candidate.length )
      target_vcf.indel( pos=candidate_pos-1, before=reference.base_at( candidate_pos-1 ), after=candidate.consensus( start=candidate_pos-1 ) )
    # TODO this doesn't work
    #elif reference.length > candidate.length: # candidate finished early
    #  target_vcf.indel( pos=pos-1, before=reference.base_at( pos-1 ), after=reference.fragment( start=pos-1 ) )

class SamToMultiChromosomeFasta(object):
  '''
    processes sam files with multiple chromosomes to generate independent pileups
  '''
  def __init__( self, sam, log, allow_indels=True ):
    self.fastas = {} # collection of sam_to_fastas by chromosome
    self.stats = { 'total_lines': 0 }
    for line in sam:
      if line.startswith( '@' ): # skip header
        pass
      else:
        line = line.strip()
        fields = line.split()
        chromosome = fields[2] # rline
        if chromosome not in self.fastas:
          self.fastas[chromosome] = SamToFasta( sam=None, log=log, allow_indels=allow_indels )
        self.fastas[chromosome].parse_line( line )
      self.stats['total_lines'] += 1
      if self.stats['total_lines'] < 10 or self.stats['total_lines'] < 5000 and self.stats['total_lines'] % 1000 == 0 or self.stats['total_lines'] % 10000 == 0:
        log( 'SamToMultiChromosomeFasta: %i lines processed' % ( self.stats['total_lines'] ) )
    # combine stats
    for chromosome in self.fastas:
      for stat in self.fastas[chromosome].stats:
        if stat not in self.stats:
          self.stats[stat] = 0
        self.stats[stat] += self.fastas[chromosome].stats[stat]
    log( self.stats )

class SamToFasta(object):
  '''
    read and evaluate assembler data
    takes a sam file and builds an imputed genome
  '''  

  def __init__( self, sam, log, allow_indels=True ):
    '''
      @sam: fh to sam file or other iterable
    '''
    self.allow_indels = allow_indels
    self.log = log
    self.fasta = fasta.ProbabilisticFasta( log )
    self.stats = { 'mapped': 0, 'unmapped': 0, 'unknown_mapping': 0, 'lines': 0, 'inversions': 0 }
    if sam is not None:
      for line in sam:
        self.parse_line( line.strip() )
        if self.stats['lines'] < 5000 and self.stats['lines'] % 1000 == 0 or self.stats['lines'] % 10000 == 0:
          self.log( 'SamToFasta: %i lines processed' % self.stats['lines'] )
          #if self.stats['lines'] == 1000000: # testing
          #  break
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

    self.stats['lines'] += 1
    if line.startswith( '@' ): # skip header
      pass
    else:
      fields = line.split()
      if len(fields) < 4:
        self.log( 'WARN: %i: unexpected format: %s' % ( self.stats['lines'], line.strip() ) )
      else:
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
          pos = int(fields[3]) - 1 # pos on reference genome at start of fragment; 0 based
          fragment_pos = 0 # position on the read
          genome_pos = 0 # position on reference
          cigar = fields[5]
          for cigar_match in re.findall( '([0-9]+)([MIDNSHP=X])', cigar ):
            cigar_len = int(cigar_match[0])
  
            if cigar_match[1] == 'M': # match
              self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + genome_pos ) #, debug=line )
              genome_pos += cigar_len # move forward on reference
              fragment_pos += cigar_len # move forward in read
  
            if self.allow_indels and cigar_match[1] == 'I': # insertion to the reference
              self.fasta.insert( fields[9][fragment_pos:fragment_pos+cigar_len], pos + genome_pos )
              fragment_pos += cigar_len # move on in fragment, don't move on in genome since it's an insertion
  
            if self.allow_indels and cigar_match[1] == 'D': # deletion from the reference
              self.fasta.delete( pos + genome_pos, cigar_len )
              genome_pos += cigar_len # don't move on in fragment, move on in reference
  
            if cigar_match[1] == 'N': # skipped region from the reference
              genome_pos += cigar_len
  
            if cigar_match[1] == 'S': # soft clipping
              #self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + fragment_pos, confidence=0.5 )
              self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + fragment_pos - cigar_len, confidence=SOFT_CLIP_CONFIDENCE ) #, debug=line )
              #genome_pos += cigar_len
              fragment_pos += cigar_len
  
            if cigar_match[1] == 'H': # hard clipping
              fragment_pos += cigar_len
  
            if cigar_match[1] == 'P': # padding
              genome_pos += cigar_len
  
            if cigar_match[1] == '=': # sequence match
              self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + genome_pos ) #, debug=line )
              genome_pos += cigar_len
              fragment_pos += cigar_len
  
            if cigar_match[1] == 'X': # sequence mismatch
              genome_pos += cigar_len
              fragment_pos += cigar_len
  
            key = 'cigar_%s' % cigar_match[1]
            if key not in self.stats:
              self.stats[key] = 0
            self.stats[key] += cigar_len

  def __repr__(self):
    return "stats: %s" % self.stats

class SamAccuracyEvaluator(object):
  '''
    evaluate sam accuracy given the correct sequence marker
  '''  
  def __init__( self, sam, variation_map=None, log=bio.log_stderr, verbose=False ):
    '''
      @sam: file like object
    '''
    self.verbose = verbose
    self.log = log
    self.variation_map = self.parse_variation_map( variation_map )
    if verbose:
      self.log( "SamAccuracyEvaluator: variation map contains %i items" % ( len(self.variation_map) ) )
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
    count = 0
    for line in vmap:
      key, value = line.strip().split()
      # key is of the form "@field:"
      actual_key = re.sub( '.*(mgsa_seq_[0-9~]*).*', '\\1', key )
      result[actual_key] = value
      if self.verbose and count < 20:
        self.log( "SamAccuracyEvaluator: added %s: %s to variation map" % ( actual_key, value ) )
      count += 1
    if self.verbose:
      self.log( "SamAccuracyEvaluator: added %i items to variation map" % count )
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
    #print "parsing %s" % line
    if line.startswith( '@' ): # skip header
      pass
    else:
      fields = line.split()
      if len(fields) < 4:
        self.log( 'WARN: %i: unexpected format: %s' % ( self.stats['lines'], line.strip() ) )
      else:
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
          original_pos = int(fields[3]) - 1 # pos in genome; conver to 0 based
          pos = original_pos
          cigar = fields[5]
          prematch = True # clipping at start
          correct_pos, correct_offset, insertion_len = [ int(x) for x in re.sub( 'mgsa_seq_([0-9~]*).*', '\\1', fields[0] ).split('~') ]
          if correct_offset != 0:
            correct_offset = correct_offset - insertion_len
          variation_map_field_name = re.sub( '(mgsa_seq_[0-9~]*).*', '\\1', fields[0] )
          if 'variation_' in fields[0]:
            variations = re.sub( '.*variation_([SID0-9,-]*).*', '\\1', fields[0] ).split(',')
            #print "field %s -> %s" % ( fields[0], variations )
          elif self.variation_map is not None and variation_map_field_name in self.variation_map:
            variations = re.sub( '([SID0-9,-]*).*', '\\1', self.variation_map[variation_map_field_name] ).split(',')
          else:
            #print "SamAccuracyEvaluator: warning: variation %s not found" % variation_map_field_name
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
          if correct_pos + correct_offset == pos:
            self.stats['correct'] += 1
            self.stats['correct_mapq'] += int(fields[4])
            for c in variations:
              if 'correct_%s' % c[0] not in self.stats:
                self.stats['correct_%s' % c[0] ] = 0
              self.stats['correct_%s' % c[0] ] += 1
          else:
            #print "incorrect: pos", pos, " correct pos", correct_pos, " correct offset", correct_offset
            self.stats['incorrect'] += 1
            self.stats['incorrect_mapq'] += int(fields[4])
            for c in variations:
              if 'incorrect_%s' % c[0] not in self.stats:
                self.stats['incorrect_%s' % c[0] ] = 0
              self.stats['incorrect_%s' % c[0] ] += 1
            if self.verbose:
              self.incorrect.append( { 
                'original_pos': original_pos, 
                'correct_pos': correct_pos, 
                'correct_offset': correct_offset, 
                'provided_pos': pos, 
                'mapq': fields[4],
                'cigar': cigar,
                'label': fields[0],
                'read': fields[9][:10],
                'variations': variations
              } )
            else:
              self.incorrect.append( { 'correct_pos': correct_pos, 'provided_pos': pos, 'mapq': fields[4] } )
            diff = pos - ( correct_pos + correct_offset )
            if diff not in self.incorrect_diff:
              self.incorrect_diff[ diff ] = 0
            self.incorrect_diff[ diff ] += 1
 
class SamDiff(object):
  def __init__( self, sam_fhs, log=bio.log_stderr ):
    self.log = log
    self.stats = collections.defaultdict(int)

    for idx, reader in enumerate(sam_fhs):
      log( 'processing %i...' % idx )
      bit_pos = 2 ** idx
      for pos, line in enumerate(reader):
        self.parse_line( pos, bit_pos, line.strip() )
        if ( pos + 1 ) % 100000 == 0:
          log( 'processed %i...' % pos )
      log( 'processing %i: read %i lines' % ( idx, pos ) )

    log( 'analyzing...' )
    self.totals = collections.defaultdict(int)
    for key, value in self.stats.items(): # readname, distribution
      self.totals[value] += 1
 
  def parse_line( self, pos, bit_pos, line ):
    fields = line.split()
    if len(fields) < 4:
      pass #self.log( 'WARN: %i: unexpected format: %s' % ( pos, line.strip() ) )
    else:
      flag = int(fields[1])
      if flag & 0x04 != 0: # unmapped
        pass 
      else:
        if flag & 0x02 != 0: # mapped
          self.stats[fields[0]] |= bit_pos
        else: # unknown mapping (assume mapped)
          self.stats[fields[0]] |= bit_pos

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

class BamReaderExternal(object):
  '''
    sam handle interface using an external BamToSam
  '''
  def __init__( self, cmd, sam_file ):
    p = subprocess.Popen(cmd % sam_file, shell=True, bufsize=0, stdin=subprocess.PIPE, stdout=subprocess.PIPE) #, close_fds=True)
    self.stdout = p.stdout

  def __iter__(self):
    return self.stdout.__iter__()
  
# class BamReader(object):
#   '''
#     provides a sam handle like interface given a bam file
#     use as an iterator.
#     basically this class returns lines from a compressed (zlib) file
#     TODO not working, incomplete
#   '''
#   MAGIC = '\x1f\x8b\x08\x04'
#   HEADER = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
#   EOF = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
#   BC = b"BC"
# 
#   def __init__( self, bam_fh ):
#     self.fh = bam_fh
#     self.buffered_lines = []
#     self.block = 0
# 
#   def __iter__(self):
#     return self
# 
#   def next(self):
#     if len(self.buffered_lines) == 0:
#       self._get_lines()
#       if len(self.buffered_lines) == 0:
#         raise StopIteration
#     return self.buffered_lines.pop(0)
# 
#   def _get_lines(self):
#     # read block header
#     magic = self.fh.read(4)
#     if not magic:
#       return # no more data
#     if magic != BamReader.MAGIC:
#       raise Exception( 'invalid data' )
#     mod_time, extra_flags, os, extra_len = struct.unpack("<LBBH", self.fh.read(8))
#     block_size = None
#     x_len = 0
#     while x_len < extra_len:
#       subfield_id = self.fh.read(2)
#       subfield_len = struct.unpack("<H", self.fh.read(2))[0]  # uint16_t
#       subfield_data = self.fh.read(subfield_len)
#       x_len += subfield_len + 4
#       if subfield_id == BamReader.BC:
#         #assert subfield_len == 2, "Wrong BC payload length"
#         #assert block_size is None, "Two BC subfields?"
#         block_size = struct.unpack("<H", subfield_data)[0] + 1  # uint16_t
#     #assert x_len == extra_len, (x_len, extra_len)
#     #assert block_size is not None, "Missing BC, this isn't a BGZF file!"
# 
#     # Now comes the compressed data, CRC, and length of uncompressed data.
#     deflate_size = block_size - 1 - extra_len - 19
#     d = zlib.decompressobj(-15)  # Negative window size means no headers
#     data = d.decompress(self.fh.read(deflate_size)) + d.flush()
#     expected_crc = self.fh.read(4)
#     expected_size = struct.unpack("<I", self.fh.read(4))[0]
#     #assert expected_size == len(data), \
#     #       "Decompressed to %i, not %i" % (len(data), expected_size)
#     # Should cope with a mix of Python platforms...
#     crc = zlib.crc32(data)
#     if crc < 0:
#         crc = struct.pack("<i", crc)
#     else:
#         crc = struct.pack("<I", crc)
# 
#     self._process_data( data )
#     #assert expected_crc == crc, \
#     #       "CRC is %s, not %s" % (crc, expected_crc)
#     #if text_mode:
#     #    return block_size, _as_string(data)
#     #else:
#     #    return block_size, data
# 
#   def find_null_terminated( self, data, start ):
#     end = start
#     while data[end] != '\0':
#       end += 1
#     return end + 1, str(data[start:end])
# 
#   def _process_data( self, data ):
#     print "processing %i bytes" % len(data)
#     if self.block == 0:
#       magic, l_text = struct.unpack( '<4sl', data[0:8] ) # uses 4 + 4
#       # header
#       offset = 8
#       #print "magic, l_text", magic, l_text
#       header, self.n_ref = struct.unpack( '<%isl' % l_text, data[offset:offset + l_text + 4] ) # uses l_text + 4
#       #print "header, n_ref", header, self.n_ref
#       if len(header) > 0:
#         self.buffered_lines = header.split( '\n' )
#       offset = offset + l_text + 4 # n_ref
#       l_name = struct.unpack( '<l', data[offset:offset+4] )[0]
#       offset += 4 # l_name
#       #print 'l_name', l_name
#       name, self.l_ref = struct.unpack( '<%isl' % l_name, data[offset:offset + l_name + 4] )
#       offset += l_name + 4
#       #self.buffered_lines.append( name )
#       #print 'offset', offset
#     else: 
#       # sequences
#       offset = 0
#       while self.l_ref > 0:
#         print 'l_ref', self.l_ref, 'offset', offset
#         if offset >= len(data):
#           break
#         start_offset = offset
#         block_size, refID, pos, bin_mq_nl, flag_nc, l_seq, next_refID, next_pos, tlen = struct.unpack( '<lllLLllll', data[offset:offset + 36] )
#         offset += 36
#         offset, read_name = self.find_null_terminated( data, offset )
#         #print 'found', read_name, 'new offset', offset
#         self.buffered_lines.append( '%s' % read_name )
# 
#         self.l_ref -= 1
#         offset = start_offset + block_size + 4
# 
#     self.block += 1
