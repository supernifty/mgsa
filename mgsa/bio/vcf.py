
import datetime
import random
import re
import sys

import bio

class VCF(object):
  '''
    manage a VCF set
  '''

  def __init__( self, reader=None, writer=None, log=bio.log_stderr ):
    '''
      @reader: file handle
      @writer: VCFWriter object
    '''
    self.log = log
    self.snp_list = []
    self.snp_map = {} # maps pos snp
    self.indel_list = []
    self.indel_map = {} # maps pos to indel
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
      if len(fields) > 5:
        pos = int( fields[1] ) - 1 # vcf file is 1-based; bio uses 0-based
        ref = fields[3]
        alt = fields[4]
        if len(ref) == 1 and len(alt) == 1:
          self.snp_map[pos] = len(self.snp_list)
          self.snp_list.append( { 'pos': pos, 'ref': ref, 'alt': alt } )
        else:
          self.indel_map[int(pos)] = len(self.indel_list)
          self.indel_list.append( { 'pos': pos, 'before': ref, 'after': alt } )
      count += 1
    self.log( 'VCF.load: %i lines loaded' % count )

  def snp( self, pos, ref, alt ):
    '''
      adds a snp
    '''
    self.snp_map[int(pos)] = len(self.snp_list)
    self.snp_list.append( { 'pos': int(pos), 'ref': ref, 'alt': alt } )
    if self.writer is not None:
      self.writer.snp( pos, ref, alt )

  def indel( self, pos, before, after ):
    '''
      adds an indel
    '''
    self.indel_map[int(pos)] = len(self.indel_list)
    self.indel_list.append( { 'pos': int(pos), 'before': before, 'after': after } )
    if self.writer is not None:
      self.writer.indel( pos, before, after )

  def write_all( self, writer ):
    '''
      writes a .vcf to the writer handle
    '''
    vcf_writer = VCFWriter(writer)
    for snp in self.snp_list:
      vcf_writer.snp( snp['pos'], snp['ref'], snp['alt'] )
    for indel in self.indel_list:
      vcf_writer.indel( indel['pos'], indel['before'], indel['after'] )

  def net_insertions(self, position):
    '''
      returns the net change in position
    '''
    net = 0
    for indel in self.indel_list: # TODO performance
      if indel['pos'] < position:
        net += len(indel['after']) - len(indel['before'])
      else:
        break # assume indel_list is ordered
    return net

  def find_indel( self, pos ):
    '''
      any indel on this position?
    '''
    # TODO need a better search
    for indel in self.indel_list:
      if indel['pos'] <= pos <= ( indel['pos'] + len(indel['after']) - len(indel['before']) ):
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

  def variations(self, start, end, include_start=False ):
    '''
      returns a string representation of all variations within the specified range
    '''
    result = set()
    start_pos = self.bisect( self.snp_list, start ) # snps must be sorted
    if start_pos < len(self.snp_list):
      for snp_idx in xrange(start_pos, len(self.snp_list)): #self.snp_list[start_pos:]:
        snp = self.snp_list[snp_idx]
        if snp['pos'] > end:
          break
        if start < snp['pos'] < end or include_start and start == snp['pos']:
          result.add( 'S%i-0' % ( snp['pos'] - start ) )
          #print snp
    start_pos = self.bisect( self.indel_list, start ) # indels must be sorted
    if start_pos < len(self.indel_list):
      for indel_idx in xrange( start_pos, len(self.indel_list) ): #self.indel_list[start_pos]:
        indel = self.indel_list[indel_idx]
        if indel['pos'] > end:
          break
        if start < indel['pos'] < end or include_start and start == indel['pos']:
          net = len(indel['after']) - len(indel['before'])
          if net < 0: # deletion
            result.add( 'D%i-%i' % ( indel['pos'] - start, -net ) )
          elif net > 0: # insertion
            indel_pos = indel['pos']
            result.add( 'I%i-%i' % ( indel_pos - start, net ) )
      #print result, " for", start, "", end
    return result

class VCFWriter(object):
  '''
    write a vcf file
  '''

  def __init__( self, writer ):
    self.writer = writer
    self.writer.write( '''##fileformat=VCFv4.1
##fileDate=%s
##source=SuperniftyMutationSimulator
##chr\tpos\tid\tref\talt\tqual\tfilter\tinfo
''' % datetime.datetime.now().strftime("%Y%m%d") )

  def snp( self, pos, ref, alt ):
    '''
      note that pos in a vcf is 1-based
    '''
    self.writer.write( '%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        # chrom, pos, id, ref, alt, qual, filter, info
        '.', (pos + 1), '.', ref, alt, '.', 'PASS', 'DP=100'
      ) )

  def indel( self, pos, before, after ):
    self.writer.write( '%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        # chrom, pos, id, ref, alt, qual, filter, info
        '.', (pos + 1), '.', before, after, '.', 'PASS', 'DP=100'
      ) )

class VCFDiff(object):
  '''
    compare two vcfs
  '''

  def __init__( self, vcf_correct, vcf_candidate, log ):
    '''
      @vcf_correct: the correct VCF
      @vcf_candidate: the candidate VCF
    '''
    self.stats = { 'tp': 0, 'fp': 0, 'fn': 0 }
    #log( '%i correct snps, %i candidate snps; %i correct snp map; %i candidate snp map' % ( len(vcf_correct.snp_list), len(vcf_candidate.snp_list), len(vcf_correct.snp_map), len(vcf_candidate.snp_map) ) )
    for true_snp in vcf_correct.snp_list:
      if true_snp['pos'] in vcf_candidate.snp_map: # there's a similarly placed snp
        candidate_snp = vcf_candidate.snp_list[vcf_candidate.snp_map[true_snp['pos']]]
        if candidate_snp['ref'] == true_snp['ref'] and candidate_snp['alt'] == true_snp['alt']:
          self.stats['tp'] += 1
        else:
          self.stats['fn'] += 1
          log( 'fn(near): %i' % true_snp['pos'] )
      else:
        self.stats['fn'] += 1
        log( 'fn: %i' % true_snp['pos'] )
      
    for candidate_snp in vcf_candidate.snp_list:
      if candidate_snp['pos'] in vcf_correct.snp_map: 
        true_snp = vcf_correct.snp_list[vcf_correct.snp_map[candidate_snp['pos']]]
        if candidate_snp['ref'] != true_snp['ref'] or candidate_snp['alt'] != true_snp['alt']:
          self.stats['fp'] += 1
          log( 'fp(near): %i' % candidate_snp['pos'] )
        else:
          pass # tp already found
      else: # false positive
        self.stats['fp'] += 1
        log( 'fp: %i' % candidate_snp['pos'] )

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
    for true_indel in vcf.indel_list: # vcf.indel has pos, before, after
      pos = true_indel['pos'] + 1
      if len(true_indel['before']) > len(true_indel['after']): # deletion
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
    for indel in vcf.indel_list:
      indel_list_positions.add( indel['pos'] )

    # precision
    for predicted_insertion in fasta.inserted:
      if predicted_insertion not in indel_list_positions:
        self.ins_stats['fp'] += 1
    for predicted_deletion in fasta.deleted:
      if predicted_deletion not in indel_list_positions:
        self.del_stats['fp'] += 1

