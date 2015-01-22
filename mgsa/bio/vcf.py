
import datetime
import math
import random
import re
import sys

import bio
import variation

# constants
DEFAULT_VCF_COVERAGE = 10
DEFAULT_VCF_ZYGOSITY = "1" # choose alt

class VCF(object):
  '''
    manage a VCF set
    @reader: fh like object
  '''

  def __init__( self, reader=None, writer=None, log=bio.log_stderr ):
    '''
      @reader: file handle
      @writer: VCFWriter object
    '''
    self.log = log
    self.manager = variation.VariationManager()
    self.snp_list = []
    self.snp_map = {} # maps pos snp
    self.writer = writer
    if reader is not None:
      self.load( reader )

  def load( self, reader ):
    '''
      @reader: file handle to vcf file
    '''
    count = 0
    for line in reader:
      if line.startswith( '#' ):
        continue
      fields = line.split()
      if len(fields) > 6:
        pos = int( fields[1] ) - 1 # vcf file is 1-based; bio uses 0-based
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        if qual == '' or qual == '.':
          confidence = 0.5
        else:
          confidence = 1. - 10 ** ( float(qual) / -10. )
        if len(ref) == 1 and len(alt) == 1:
          self.snp_map[pos] = len(self.snp_list)
          self.snp_list.append( { 'pos': pos, 'ref': ref, 'alt': alt, 'conf': confidence } )
        else:
          self.manager.indel_map[int(pos)] = len(self.manager.indel_list)
          self.manager.indel_list.append( variation.IndelVariation( pos, ref, alt ) )
      count += 1
    self.log( 'VCF.load: %i lines loaded' % count )

  def snp( self, pos, ref, alt, coverage=None, confidence=0.5 ):
    '''
      adds a snp
    '''
    self.snp_map[int(pos)] = len(self.snp_list)
    self.snp_list.append( { 'pos': int(pos), 'ref': ref, 'alt': alt, 'conf': confidence } )
    if self.writer is not None:
      self.writer.snp( pos, ref, alt, coverage=coverage, confidence=confidence )

  def indel( self, pos, before, after, coverage=None, confidence=0.5 ):
    '''
      adds an indel
    '''
    self.manager.indel( pos, before, after )
    if self.writer is not None:
      self.writer.indel( pos, before, after, coverage=coverage, confidence=confidence )

  def write_all( self, writer ):
    '''
      writes a .vcf to the writer handle
    '''
    vcf_writer = VCFWriter(writer)
    for snp in self.snp_list:
      vcf_writer.snp( snp['pos'], snp['ref'], snp['alt'] )
    for indel in self.manager.indel_list:
      vcf_writer.indel( indel.pos, indel.before, indel.after )

  def net_insertions_to_reference_position(self, position):
    '''
      returns the net change in position given a reference position
    '''
    net = 0
    for indel in self.manager.indel_list: # TODO performance
      if indel.pos < position:
        net += len(indel.after) - len(indel.before)
      else:
        break # assume indel_list is ordered
    return net

  def net_insertions_to_candidate_position(self, candidate_position):
    '''
      returns the net change in position given a candidate position
    '''
    net_insertions = 0
    for indel in self.manager.indel_list: # TODO performance
      current_reference_position = indel.pos
      current_candidate_position = indel.pos + net_insertions
      if current_candidate_position < candidate_position:
        net_insertions += len(indel.after) - len(indel.before)
      else:
        break # assume indel_list is ordered
    return net_insertions

  def candidate_position_to_reference_position(self, candidate_position):
    '''
      returns the corresponding reference position, plus any offset back to the reference position
      e.g. ref = GA; candidate = GTTA; indel should be { pos 2 before GA after GTTA }
      0 -> (0, 0)
      1 -> (0, 1)
      2 -> (0, 2)
      3 -> (1, 0)
    '''
    return self.manager.candidate_position_to_reference_position( candidate_position )

  def find_indel( self, pos ):
    '''
      any indel on this position?
    '''
    # TODO need a better search
    for indel in self.manager.indel_list:
      if indel.pos <= pos <= ( indel.pos + len(indel.after) - len(indel.before) ):
        return indel
    return None
    #if pos in self.indel_map:
    #  return self.indel_list[self.indel_map[pos]]
    #else:
    #  return None

  def bisect(self, pos_list, value):
    '''returns an index that is less than or equal to value'''
    start = 0 # inclusive
    end = len(pos_list) # exclusive
    while end - start > 1:
      mid = start + ( end - start ) / 2
      if pos_list[mid]['pos'] == value:
        return mid
      elif pos_list[mid]['pos'] < value:
        start = mid
      else: # pos > value
        end = mid
      #print "bisect: start: %i end: %i" % ( start, end )
        
    return start

  def variations(self, start, end, offset=0 ):
    '''
      returns a string representation of all variations within the specified range.
      start and end are both positions on the *reference*
    '''
    #print "variations for %i %i" % ( start, end )
    result = set()
    start_pos = self.bisect( self.snp_list, start ) # snps must be sorted
    if start_pos < len(self.snp_list):
      for snp_idx in xrange(start_pos, len(self.snp_list)): #self.snp_list[start_pos:]:
        snp = self.snp_list[snp_idx]
        if snp['pos'] > end:
          break
        if start <= snp['pos'] < end:
          result.add( 'S%i-0' % ( snp['pos'] - start ) )
          #print snp
    #start_pos = self.bisect( self.indel_list, start ) # indels must be sorted
    #if start_pos < len(self.indel_list):
    #  for indel_idx in xrange( start_pos, len(self.indel_list) ): #self.indel_list[start_pos]:
    #    indel = self.indel_list[indel_idx]
    #    if indel['pos'] > end:
    #      break
    #    if start < indel['pos'] < end or include_start and start == indel['pos']:
    #      net = len(indel['after']) - len(indel['before'])
    #      if net < 0: # deletion
    #        result.add( 'D%i-%i' % ( indel['pos'] - start, -net ) )
    #      elif net > 0: # insertion
    #        indel_pos = indel['pos']
    #        result.add( 'I%i-%i' % ( indel_pos - start, net ) )
    if offset > 0: # offset means we are in an insertion, don't include again
      start_search = start + 2
      #include_start = False
    else:
      start_search = start + 1 # no offset means we are past any insertion, don't include again

    #for indel_idx in xrange( 0, len(self.manager.indel_list) ):
    indels_overlapping = self.manager.find_indels_overlapping( start_search, end )
    for indel in indels_overlapping:
      net = len(indel.after) - len(indel.before)
      #if indel.pos >= start_search and indel.pos < end: # the indel starts inside the read
      if net < 0: # deletion
        result.add( 'D%i-%i' % ( indel.pos - start, -net ) )
      elif net > 0: # insertion
        indel_pos = indel.pos
        result.add( 'I%i-%i[%i,%i,%i]' % ( indel_pos - start, net, indel_pos, offset, end ) )
          #print "added full indel pos %i rel %i net %i" % ( indel_pos, indel_pos - start, net )
      #if indel['pos'] < start and indel['pos'] + net >= start: # the read starts inside an indel
      #  result.add( 'I0-%i' % ( net - ( start - indel['pos'] ) ) ) # partial insert
        #print "added partial indel 0 %i" % ( start - indel['pos'] )
      #print result, " for", start, "", end
    return result

  def __str__( self ):
    return 'snps: %s indels: %s' % ( self.snp_list, self.manager )

class VCFWriter(object):
  '''
    write a vcf file
  '''

  def __init__( self, writer ):
    '''
      @writer: file like object
    '''
    self.writer = writer
    self.writer.write( '''##fileformat=VCFv4.1
##fileDate=%s
##source=SuperniftyMutationSimulator
##chr\tpos\tid\tref\talt\tqual\tfilter\tinfo\tformat\tna00001
''' % datetime.datetime.now().strftime("%Y%m%d") )

  def snp( self, pos, ref, alt, confidence=0.5, coverage=None, zygosity=None ):
    '''
      note that pos in a vcf is 1-based
      @coverage: depth of coverage on this variation
      @zygosity: for diploid, "0/1" = heterozygous; "1/1" = alt homozygous
    '''
    if confidence < 1.:
      quality = '%.0f' % ( -10 * math.log(1.-confidence, 10 ))
    else:
      quality = '.'

    self.writer.write( '%s\t%i\t%s\t%s\t%s\t%s\t%s\tDP=%i\tGT\t%s\n' % (
        # chrom, pos, id, ref, alt, qual, filter, info
        '.', (pos + 1), '.', ref, alt, quality, 'PASS', coverage or DEFAULT_VCF_COVERAGE, zygosity or DEFAULT_VCF_ZYGOSITY
      ) )

  def indel( self, pos, before, after, coverage=None, zygosity=None ):
    self.writer.write( '%s\t%i\t%s\t%s\t%s\t%s\t%s\tDP=%i\tGT\t%s\n' % (
        # chrom, pos, id, ref, alt, qual, filter, info
        '.', (pos + 1), '.', before, after, '.', 'PASS', coverage or DEFAULT_VCF_COVERAGE, zygosity or DEFAULT_VCF_ZYGOSITY
      ) )

class VCFDiff(object):
  '''
    compare two vcfs
  '''

  def __init__( self, vcf_correct, vcf_candidate, log=bio.log_stderr, generate_positions=False ):
    '''
      generates self.stats and self.buckets, containing tp, fp, and fn
      @vcf_correct: the correct VCF()
      @vcf_candidate: the candidate VCF()
    '''
    self.stats = { 'tp': 0, 'fp': 0, 'fn': 0 }
    self.generate_positions = generate_positions
    self.positions = { 'tp': [], 'fp': [], 'fn': [] }
    self.buckets = [ { 'tp': 0, 'fp': 0, 'fn': 0 } for _ in xrange(0, 20) ]
    max_snp_pos = 0
    if len( vcf_correct.snp_list ) > 0:
      max_snp_pos = vcf_correct.snp_list[-1]['pos']
    if len( vcf_candidate.snp_list ) > 0 and vcf_candidate.snp_list[-1]['pos'] > max_snp_pos:
      max_snp_pos = vcf_candidate.snp_list[-1]['pos']
    max_indel_pos = 0
    if len( vcf_correct.manager.indel_list ) > 0:
      max_snp_pos = vcf_correct.manager.indel_list[-1].pos
    if len( vcf_candidate.manager.indel_list ) > 0 and vcf_candidate.manager.indel_list[-1].pos > max_indel_pos:
      max_indel_pos = vcf_candidate.manager.indel_list[-1].pos
    bucket_size = ( 1 + max( max_snp_pos, max_indel_pos ) ) / 20.
    log( '%i correct snps, %i candidate snps; %i correct snp map; %i candidate snp map' % ( len(vcf_correct.snp_list), len(vcf_candidate.snp_list), len(vcf_correct.snp_map), len(vcf_candidate.snp_map) ) )
    log( '%i correct indels, %i candidate indels; %i correct indels map; %i candidate indels map' % ( len(vcf_correct.manager.indel_list), len(vcf_candidate.manager.indel_list), len(vcf_correct.manager.indel_map), len(vcf_candidate.manager.indel_map) ) )
    # snps
    for true_snp in vcf_correct.snp_list:
      bucket = int( math.floor( true_snp['pos'] / bucket_size ) )
      if true_snp['pos'] in vcf_candidate.snp_map: # there's a similarly placed snp
        candidate_snp = vcf_candidate.snp_list[vcf_candidate.snp_map[true_snp['pos']]]
        if candidate_snp['ref'] == true_snp['ref'] and candidate_snp['alt'] == true_snp['alt']:
          self.stats['tp'] += 1
          self.buckets[bucket]['tp'] += 1
          if self.generate_positions:
            self.positions['tp'].append( true_snp['pos'] )
        else:
          self.stats['fn'] += 1
          self.buckets[bucket]['fn'] += 1
          #log( 'fn(near): %i' % true_snp['pos'] )
          if self.generate_positions:
            self.positions['fn'].append( true_snp['pos'] )
      else:
        self.stats['fn'] += 1
        self.buckets[bucket]['fn'] += 1
        if self.generate_positions:
          self.positions['fn'].append( true_snp['pos'] )
        log( 'fn: %i' % true_snp['pos'] )
      
    for candidate_snp in vcf_candidate.snp_list:
      bucket = int( math.floor( candidate_snp['pos'] / bucket_size ) )
      if candidate_snp['pos'] in vcf_correct.snp_map: 
        true_snp = vcf_correct.snp_list[vcf_correct.snp_map[candidate_snp['pos']]]
        if candidate_snp['ref'] != true_snp['ref'] or candidate_snp['alt'] != true_snp['alt']:
          self.stats['fp'] += 1
          self.buckets[bucket]['fp'] += 1
          if self.generate_positions:
            self.positions['fp'].append( candidate_snp['pos'] )
          #log( 'fp(near): %i' % candidate_snp['pos'] )
        else:
          pass # tp already found
      else: # false positive
        self.stats['fp'] += 1
        self.buckets[bucket]['fp'] += 1
        if self.generate_positions:
          self.positions['fp'].append( candidate_snp['pos'] )
        log( 'fp: %i' % ( candidate_snp['pos'] ) ) # debug

    # indels - true positives
    for true_indel in vcf_correct.manager.indel_list:
      bucket = int( math.floor( true_indel.pos / bucket_size ) )
      if vcf_candidate.manager.find_indel_match( true_indel ) is None:
        self.stats['fn'] += 1
        self.buckets[bucket]['fn'] += 1
        log( 'fn: %s' % true_indel )
      else:
        self.stats['tp'] += 1
        self.buckets[bucket]['tp'] += 1

    # indels - false positives
    for candidate_indel in vcf_candidate.manager.indel_list:
      bucket = int( math.floor( candidate_indel.pos / bucket_size ) )
      if vcf_correct.manager.find_indel_match( candidate_indel ) is None:
        self.stats['fp'] += 1
        self.buckets[bucket]['fp'] += 1
        log( 'fp: %s' % candidate_indel )
      else:
        pass # tp already found

class VCFFastaDiff(object):
  '''
    given a vcf from the reference and a probabilisticfasta from the candidate, 
    calculate found/not found variations
    i.e. compare variations, calculate confusion matrix
  '''
  def __init__( self, vcf, fasta, log ):
    '''
      @vcf: VCF object
      @fasta: ProbabilisticFasta object
    '''
    self.snp_stats = { 'tp': 0, 'fn': 0, 'other': 0 }
    for true_snp in vcf.snp_list: # vcf.snp has pos, ref, alt
      pos = true_snp['pos']
      candidate_value = fasta.consensus_at(pos)[1]
      if candidate_value == true_snp['ref']: # missed snp
        self.snp_stats['fn'] += 1
      elif candidate_value == true_snp['alt']: # found snp
        self.snp_stats['tp'] += 1
      else: # neither snp ref or alt
        log( 'consensus at %i: %s %s %s ref %s alt %s' % ( 
          pos, 
          fasta.consensus_at(pos-1)[1], 
          candidate_value, 
          fasta.consensus_at(pos+1)[1], 
          true_snp['ref'], true_snp['alt'] ) )
        self.snp_stats['other'] += 1

    self.ins_stats = { 'tp': 0, 'fn': 0, 'fp': 0, 'other': 0 }
    self.del_stats = { 'tp': 0, 'fn': 0, 'fp': 0, 'other': 0 }
    # recall
    for true_indel in vcf.manager.indel_list: # vcf.indel has pos, before, after
      pos = true_indel.pos + 1
      if len(true_indel.before) > len(true_indel.after): # deletion
        if pos in fasta.deleted:
          self.del_stats['tp'] += 1
        else: # not found
          self.del_stats['fn'] += 1
      else: # insertion
        if pos in fasta.inserted:
          self.ins_stats['tp'] += 1
        else: # not found
          self.ins_stats['fn'] += 1

    indel_list_positions = set()
    for indel in vcf.manager.indel_list:
      indel_list_positions.add( indel.pos )

    # precision
    for predicted_insertion in fasta.inserted:
      if predicted_insertion not in indel_list_positions:
        self.ins_stats['fp'] += 1
    for predicted_deletion in fasta.deleted:
      if predicted_deletion not in indel_list_positions:
        self.del_stats['fp'] += 1

