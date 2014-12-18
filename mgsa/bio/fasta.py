
import collections
import datetime
import random
import re
import sys

import bio
import vcf

DELETE_BASE = '-'

class ProbabilisticFasta(object):
  '''
    generate fasta with probabilities attached from mapped fragments
  '''
  
  def __init__( self, log=None ):
    self.genome = {} # counts over position { 'A': [ 1, 0, 2, ... ], 'G': [ 2, 2, 0, ... ] } roughly means [AGG][GG][AA]...
    self.insertions = {} # { 0: { 'AAA': 1, 'BB': 2 } } means at position 0, insertion of AAA with 1 confidence; BB with 2 confidence
    self.inserted = {}
    self.deleted = set()
    self.prior = {} # overall counts
    self.total = 0 # total # of nucleotides
    self.length = 0 # length of genome
    self.log = log

  def add( self, fragment, start, confidence=1.0 ):
    '''
      process a directly mapped fragment of dna
      @fragment: the piece of dna
      @fragment: where in the genome it starts
      @confidence: how confident in this mapping
    '''
    for i in xrange(0, len(fragment) ):
      value = fragment[i]
      position = start + i
      if value not in self.genome:
        self.genome[value] = []
        self.prior[value] = 0
      if len(self.genome[value]) < position + 1:
        self.genome[value].extend( [0] * ( position - len(self.genome[value]) + 1 ) )
        if len(self.genome[value]) > self.length:
          self.length = len(self.genome[value])
      #self.log( "adding confidence %f to %s at %i; start %i" % ( confidence, value, position, start ) )
      #self.log( "capacity %i" % ( len(self.genome[value] ) ) )
      if position >= 0: # mappers can clip off the end
        self.genome[value][position] += confidence
        self.prior[value] += 1
        self.total += 1

  def insert( self, fragment, start, confidence=1.0 ):
    '''
      @fragment: the piece of dna to insert
      @start: where the insertion starts on the reference
      @confidence: how confident in this insertion
    '''
    if start not in self.insertions:
      self.insertions[start] = {}
    if fragment not in self.insertions[start]:
      self.insertions[start] = { fragment: 0 }
    self.insertions[start][fragment] += confidence

  def delete( self, start, count, confidence=1.0 ):
    fragment = DELETE_BASE * count
    self.add( fragment, start, confidence )

  def count( self, position ):
    '''
      counts of nucleotides seen at a position
      e.g. { 'A': 10, 'G': 4 }
    '''
    result = collections.defaultdict(int)
    for value in self.genome:
      if len(self.genome[value]) > position:
        result[value] = self.genome[value][position]
    return result

  def confidence( self, position ):
    '''
      returns probabilities of nucleotides
      TODO priors
    '''
    counts = self.count( position ) # e.g. { 'A': 4, 'C': 1 }
    denominator = sum( counts[k] for k in counts ) + len(self.prior)
    probabilities = {}
    for count in counts:
      probabilities[count] = ( counts[count] + 1.0 ) / denominator
    return probabilities

  def consensus_at( self, i=0 ):
    '''
      returns the majority base at a given position, or N if no coverage (move, best, confidence, coverage)
    '''
    best = 'N' # default if no coverage
    best_value = 0
    coverage = 0
    result = ''
    move = 0
    # consensus at this base
    for value in self.genome:
      if i < len(self.genome[value]):
        coverage += self.genome[value][i]
        if self.genome[value][i] > best_value:
          best_value = self.genome[value][i]
          best = value
    # find consensus insertion
    if i in self.insertions:
      # pick highest
      best_insertion = None
      best_insertion_value = 0
      insertion_coverage = 0
      for insertion in self.insertions[i]:
        insertion_coverage += self.insertions[i][insertion] 
        if best_insertion is None or self.insertions[i][insertion] > best_insertion_value:
          best_insertion = insertion
          best_insertion_value = self.insertions[i][insertion] 
      if best_insertion is not None and best_insertion_value > ( coverage - insertion_coverage ):
        result = best_insertion
        self.inserted[i] = best_insertion
        #if self.log:
        #  self.log( 'included insertion %s at ref %i' % (best_insertion, i) )
        move = len(best_insertion)
      else:
        pass
        if self.log:
          pass #self.log( 'skipped insertion at %i with val %f with noinsert %f' % ( i, best_insertion_value, coverage - insertion_coverage ) )
    # was the consensus to delete?
    if best == DELETE_BASE:
      self.deleted.add(i)
      #if self.log:
      #  self.log( 'deletion at %i' % (i) )
    else:
      result += best
      move += 1
     
    return (move, result, best, coverage)
  
  def consensus_count( self, start=0, count=1 ):
    '''
      get the next count characters from start
      can return more than count if insertions at the end of the segment
      can return less than count if at end of string
    '''
    result = ''
    total_move = 0
    reference_move = 0
    while len(result) <= count and start + reference_move < self.length:
      move, s, confidence, coverage = self.consensus_at( start + reference_move )
      result += s 
      reference_move += 1
    #self.log( 'consensus for %i count %i is %s (%i)' % (start, count, result, len(result)) )
    return (reference_move, result)

  def consensus( self, start=0, end=-1 ):
    '''
      return the corresponding string for the segment specified.
      returned segment can be longer or shorter than end-start - segment applies to the reference.
    '''
    if end == -1 or end > self.length:
      end = self.length
    result = ''
    self.insertion_count = 0
    for i in xrange(start, end):
      move, add, confidence, coverage = self.consensus_at( i )
      result += add
    #self.log( 'consensus for %i to %i is %s (%i)' % (start, end, result, len(result)) )
    return result 

class FastaMutate(object):
  '''
    change a reference fasta
  '''
  probabilities = 'AAACCTTGGG'

  def __init__( self, reader, log=bio.log_stderr, vcf_file=None, snp_prob=0.01, insert_prob=0.01, delete_prob=0.01, min_insert_len=1, max_insert_len=1, min_delete_len=1, max_delete_len=1, min_variation_dist=0, probabilistic=True, insert_source='random', allow_end_mutate=False ):
    '''
      @reader: FastaReader
      @vcf_file: write mutations to vcf
    '''
    self.reader = reader
    self.snp_prob = snp_prob
    self.insert_prob = insert_prob
    self.delete_prob = delete_prob
    self.min_insert_len = min_insert_len
    self.max_insert_len = max_insert_len
    self.min_delete_len = min_delete_len
    self.max_delete_len = max_delete_len
    self.min_variation_dist = min_variation_dist
    self.deletion_remain = 0
    self.mutations = 0
    self.vcf_file = vcf_file
    self.probabilistic = probabilistic
    self.insert_source = insert_source
    if vcf_file is not None:
      self.vcf = vcf.VCF( writer=vcf.VCFWriter(vcf_file) )
    else:
      self.vcf = vcf.VCF()
    self.pos = 0
    self.last_variation_pos = None # this is the end position of the last variation
    self.allow_end_mutate = allow_end_mutate

    seed = random.randint(0, sys.maxint)
    random.seed(seed)
    if log is not None:
      log( 'seed: %i' % seed )
  
  def items(self):
    while True:
      fragment = self.reader.next_item()
      if fragment is None: # no more fragments
        if self.deletion_remain > 0:
          self.end_deletion()
        break
      # apply mutations
      if self.allow_end_mutate or self.reader.has_next_item(): # don't mutate last fragment
        fragment = self.mutate( fragment )
      yield fragment

  def add_snp(self, c):
    '''
      generates a single snp, adds to vcf, returns new base
    '''
    new_c = self.probabilities[random.randint(0, len(self.probabilities)-1)]
    while new_c == c:
      new_c = self.probabilities[random.randint(0, len(self.probabilities)-1)]
    self.mutations += 1
    if self.vcf is not None:
      self.vcf.snp( self.pos, c, new_c )
    return new_c

  def add_insertion(self, c):
    '''
      generates a new insertion and returns it
    '''
    insert_len = random.randint(self.min_insert_len, self.max_insert_len) # decide insertion len
    # generate actual insertion
    if self.insert_source == 'random':
      new_c = ''
      while insert_len > 0:
        insert_len -= 1
        new_c += self.probabilities[random.randint(0, len(self.probabilities)-1)]
    elif self.insert_source == 'repeated':
      pass
    elif self.insert_source == 'novel':
      pass
    self.mutations += 1
    # add to vcf
    if self.vcf is not None:
      self.vcf.indel( self.pos, self.previous + c, self.previous + new_c + c )
    return new_c

  def add_deletion(self, c):
    self.mutations += 1
    self.deletion_remain = random.randint(self.min_delete_len, self.max_delete_len)
    self.deleted = ''
    self.deletion_start = self.pos
    self.deletion_previous = self.previous
    self.continue_deletion( c )

  def continue_deletion(self, c):
    self.deletion_remain -= 1
    self.deleted += c
    if self.deletion_remain == 0:
      self.end_deletion()

  def end_deletion(self):
      if self.vcf is not None:
        self.vcf.indel( self.deletion_start - 1, self.deletion_previous + self.deleted, self.deletion_previous )

  def mutate(self, fragment):
    result = ''
    for c in fragment: # iterate over each base in fragment
      if self.deletion_remain > 0:
        self.continue_deletion( c )
      elif self.probabilistic:
        # snp
        if random.uniform(0, 1) < self.snp_prob and ( self.last_variation_pos is None or self.last_variation_pos + self.min_variation_dist <= self.pos ):
          new_c = self.add_snp( c )
          result += new_c
          self.last_variation_pos = self.pos
        # insert
        elif random.uniform(0, 1) < self.insert_prob and self.pos > self.max_insert_len and ( self.last_variation_pos is None or self.last_variation_pos + self.min_variation_dist <= self.pos ): # TODO reads can get -ve reference
          new_c = self.add_insertion( c )
          result += new_c + c # insertion gets placed before current base
          self.last_variation_pos = self.pos - 1 # -1 because insertion is placed before current
        # delete
        elif self.pos > 0 and random.uniform(0, 1) < self.delete_prob and ( self.last_variation_pos is None or self.last_variation_pos + self.min_variation_dist <= self.pos ): 
          self.add_deletion( c )
          self.last_variation_pos = self.pos + self.deletion_remain
        # no mutation
        else: 
          result += c
      else: # deterministic
        pass
  
      self.pos += 1
      self.previous = c
    return result

class FastaDiff(object):
  '''
    compare fastas
  '''

  def __init__( self, reader, candidate, log, vcf=None ):
    '''
      @reader: FastaReader of known donor
      @candidate: ProbabilisticFasta of imputed genome
      @vcf: use to keep fastas properly aligned
    '''
    self.candidate_pos = 1
    self.last_log_pos = 0
    self.last_error_total_logged = 0
    self.errors = {}
    self.error_positions = set() # set of incorrect positions
    self.error_total = 0
    self.candidate = candidate
    self.log = log
    self.vcf = vcf
    candidate_remainder = ''
    for item in reader.items():
      # gets the corresponding segment from the candidate
      move, candidate_item = candidate.consensus_count( self.candidate_pos, len(item) )
      candidate_item = candidate_remainder + candidate_item
      candidate_remainder = candidate_item[len(item):]
      #self.log( 'truth: %s pred:  %s' % ( item, candidate_item ) ) # debug
      if candidate_item != item:
        # find differences
        self.find_differences( item, candidate_item, self.candidate_pos )
      self.candidate_pos += move #len(item)
      if self.candidate_pos < 10000 and self.candidate_pos - self.last_log_pos >= 1000 or self.candidate_pos - self.last_log_pos >= 10000:
        self.log_error()

  def log_error( self, i=0 ):
    self.log( 'processed cand %i items - %i errors' % ( self.candidate_pos + i, self.error_total ) )
    self.last_log_pos = self.candidate_pos
    self.last_error_total_logged = self.error_total

  def find_differences( self, item, candidate, start ):
    for i in xrange(0, min(len(item), len(candidate))):
      if item[i] != candidate[i]:
        self.error_total += 1
        if self.error_total < 100 and self.error_total != self.last_error_total_logged:
          self.log_error(i)
        key = '%s->%s' % ( item[i], candidate[i] )
        if self.error_total < 100 or self.error_total % 1000 == 0:
          self.log( 'error: actual->predicted %s at cand %i + %i (%i): counts %s\ntruth: %s\npred:  %s' % ( key, start, i, start + i, self.candidate.count( start + i ), item, candidate ) )
        if key not in self.errors:
          self.errors[key] = 0
        self.errors[key] += 1
        self.error_positions.add( start + i )

    if len(item) != len(candidate):
      key = 'length'
      if key not in self.errors:
        self.errors[key] = 0
      self.errors[key] += 1
    
class Fasta(object):
  '''
    keeps the whole fasta in memory for random access
  '''
  def __init__(self, reader):
    self.reader = reader
    self._length = None
    self.fasta = ''

  @property
  def length(self):
    if self._length is None:
      # read the whole thing
      #print "getlength"
      for item in self.reader.items():
         #print "lgot", item
         self.fasta += item
      self._length = len(self.fasta)
    #print "getlength out"
    return self._length

  def base_at(self, i):
    if i >= len(self.fasta):
      for item in self.reader.items():
        #print "got", item
        self.fasta += item
        if i < len(self.fasta): # success
          break
    if i < len(self.fasta):
      return self.fasta[i]
    else:
      return None

  def fragment(self, start=0, end=-1 ):
    if end == -1:
      end = self.length
    result = ''
    for i in xrange(start, end):
      result += self.base_at(i)
    return result

class FastaReader(object):
  '''
    yields fragments from a fasta file object
>>> import StringIO
>>> g = StringIO.StringIO( '>\nabc\ndef' )
>>> p = bio.FastaReader( g )
>>> [ f for f in p.items() ]
['abc', 'def']
  '''
  def __init__(self, genome):
    self.genome = genome
    self.has_next = True
    self.future_fragment = self._next_item()

  def items(self):
    while True:
      current_fragment = self.future_fragment
      if current_fragment is None:
        break
      self.future_fragment = self._next_item()
      self.has_next = self.future_fragment is not None
      #print "current %s future %s" % ( current_fragment, self.future_fragment )
      yield current_fragment

  #def items(self):
  #  while True:
  #    current_fragment = self.next_item()
  #    if current_fragment is None:
  #      break
  #    yield current_fragment

  def next_item(self):
    current_fragment = self.future_fragment
    self.future_fragment = self._next_item()
    self.has_next = self.future_fragment is not None
    return current_fragment    
 
  def _next_item(self):
    for line in self.genome:
      if not line.startswith( '>' ):
        return line.strip()
    return None

  def has_next_item(self):
    return self.has_next

class FastaStats(object):
  '''calculate some overall stats for a fasta file'''
  def __init__( self, fasta, log=bio.log_stderr ):
    self.stats = { 'count': 0 }
    for line in fasta:
      line = line.strip()
      if not line.startswith( '>' ):
        self.stats['count'] += len(line)
    log( self.stats )

class ErrorGenerator(object):
  def __init__( self, error_profile ):
    self.error_profile = error_profile

  def apply_errors( self, dna ):
    result = []
    for x in dna:
      result.append( self.error_profile( x ) )
    return ''.join( result )

  @staticmethod
  def create_uniform_error_profile( error_prob ):
    transitions = { 'A': 'TGC', 'T': 'GCA', 'G': 'ACT', 'C': 'AGT', 'N': 'N' }
    def uniform_error_profile( bp ):
      if random.random() < error_prob:
        return transitions[bp][random.randint(0, len(transitions[bp]) - 1)]
      else:
        return bp
    return uniform_error_profile

if __name__ == "__main__":
  #s = SamToFasta( sys.stdin, log )
  #log( 'size is %i' % s.fasta.size )

  #log( s.fasta.prior )
  #for i in xrange(1, 10):
  #  print i, s.fasta.count(i), s.fasta.confidence(i)

  #pos = 1
  #while pos < s.fasta.size:
  #  print s.fasta.consensus( pos, pos + 70 )
  #  pos += 70

  #d = FastaDiff( FastaReader( open( '../data/test/test1k.fasta', 'r' ) ), s.fasta, log )
  #log( d.error_total )
  #log( d.errors )

  m = FastaMutate( FastaReader( open( '../data/test/test1k.fasta', 'r' ) ), log, open( '../data/test/test1ksnp.vcf', 'w' ), 0.1 )
  o = open( '../data/test/test1ksnp.fasta', 'w' )
  o.write( '>mutated-snp-test1k\n')
  for line in m.items():
    o.write( '%s\n' % line )
  o.close()
  log( m.mutations )
