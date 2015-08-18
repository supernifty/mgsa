
import collections
import datetime
import random
import re
import StringIO
import sys

import bio
import features
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

  def add( self, fragment, start, confidence=1.0, debug=None ):
    '''
      process a directly mapped fragment of dna
      @fragment: the piece of dna
      @start: where in the reference genome it starts
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
        #if position == 39415: # debugging
        #  self.log( 'start %i: added %f to %s, total %f, line %s' % ( start, confidence, value, self.genome[value][position], debug ) )

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
    denominator = sum( counts[k] for k in counts ) + len(self.prior) # laplacian smoothing
    probabilities = {}
    for count in counts:
      probabilities[count] = ( counts[count] + 1.0 ) / denominator
    return probabilities

  def consensus_at( self, i=0 ):
    '''
      returns the majority base at given position i, or N if no coverage (move, best, confidence, coverage)
      @return move, result, best, coverage
       - move: how many bases to move candidate (1 is normal)
       - result: best variation
       - best: confidence of this variation
       - coverage: total # reads covering this base
    '''
    best = 'N' # default if no coverage
    best_value = 0
    coverage = 0
    result = ''
    move = 0
    # consensus at this base
    for value in self.genome: # i.e. A, C, G, T
      if i < len(self.genome[value]): # check in range
        coverage += self.genome[value][i]
        if self.genome[value][i] > best_value:
          best_value = self.genome[value][i] # best count
          best = value # best base
    #self.log( 'consensus at %i: best %s: %i of %i' % ( i, best, best_value, coverage ) )
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
     
    return (move, result, best_value, coverage) # TODO only returning confidence for variation, not indel
  
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
      inclusive of start, exclusive of end
      
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

class MultiFastaMutate(object):
  def __init__( self, multi_reader, log=bio.log_stderr, vcf_file=None, snp_prob=0.01, insert_prob=0.01, delete_prob=0.01, min_insert_len=1, max_insert_len=1, min_delete_len=1, max_delete_len=1, min_variation_dist=0, min_variation_start=0, probabilistic=True, insert_source='random', allow_end_mutate=False, tandem_count=1 ):
    for reader in multi_reader.items():
      FastaMutate( reader, log, vcf_file, snp_prob, insert_prob, delete_prob, min_insert_len, max_insert_len, min_delete_len, max_delete_len, min_variation_dist, min_variation_start, probabilistic, insert_source, allow_end_mutate, tandem_count=1 )

class FastaMutate(object):
  '''
    change a reference fasta
  '''

  def __init__( self, reader, log=bio.log_stderr, vcf_file=None, snp_prob=0.01, insert_prob=0.01, delete_prob=0.01, min_insert_len=1, max_insert_len=1, min_delete_len=1, max_delete_len=1, min_variation_dist=0, min_variation_start=0, probabilistic=True, insert_source='random', allow_end_mutate=False, probabilities='AACCCTTGGG', tandem_count=1 ):
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
    self.min_variation_start = min_variation_start
    self.deletion_remain = 0
    self.mutations = 0
    self.vcf_file = vcf_file
    self.probabilistic = probabilistic
    self.probabilities = probabilities
    self.insert_source = insert_source
    self.tandem_count = tandem_count
    self.insert_source_data = ''
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
      log( 'seed: %i; insert_source %s' % (seed, insert_source) )
    self.log = log

  def items(self):
    while True:
      fragment = self.reader.next_item()
      if fragment is None: # no more fragments
        if self.deletion_remain > 0:
          self.end_deletion()
        break
      if self.insert_source == 'repeat':
        self.insert_source_data += str(fragment)
      elif self.insert_source == 'tandem':
        self.insert_source_data += str(fragment)
        if len(self.insert_source_data) > len(str(fragment)) + self.max_insert_len: # keep this fragment + length
          self.insert_source_data = self.insert_source_data[-len(str(fragment)) - self.max_insert_len:]
        #self.log( 'insert_source_data: %s' % self.insert_source_data )
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

  def add_insertion(self, c, fragment_pos, fragment_len ):
    '''
      generates a new insertion and returns it
      insertion gets placed before current base
    '''
    insert_len = random.randint(self.min_insert_len, self.max_insert_len) # decide insertion len
    # generate actual insertion
    if self.insert_source == 'repeat' and len(self.insert_source_data) >= insert_len:
      #self.log( 'choosing repeated fragment up to %i' % len(self.insert_source_data) )
      fragment_start = random.randint(0, len(self.insert_source_data) - insert_len)
      new_c = self.insert_source_data[fragment_start:fragment_start + insert_len]
      #self.log( "repeat %s from %i:%i" % ( new_c, fragment, position ) )
    elif self.insert_source == 'tandem' and len(self.insert_source_data) >= insert_len + fragment_len:
      #self.log( 'choosing tandem fragment up to %i from %s' % ( len(self.insert_source_data), self.insert_source_data ) )
      # want the insert_len chars directly before fragment_pos, where data is len fragment_len + insert_len
      dist_from_end = fragment_len - fragment_pos
      new_c = self.insert_source_data[ -insert_len -dist_from_end: -dist_from_end ]
    elif self.insert_source == 'random':
      new_c = ''
      while insert_len > 0:
        insert_len -= 1
        new_c += self.probabilities[random.randint(0, len(self.probabilities)-1)]
      #self.log( "random %s" % ( new_c ) )
    else: #if self.insert_source == 'simple':
      possibles = 'ACGT'
      new_c = possibles[random.randint(0, len(possibles)-1)] * insert_len
    #elif self.insert_source == 'novel':
    #  pass
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
    for fragment_pos, c in enumerate(fragment): # iterate over each base in fragment
      #self.log( 'evaluating %s' % c )
      if self.deletion_remain > 0:
        self.continue_deletion( c )
      elif self.probabilistic:
        # snp
        if random.uniform(0, 1) < self.snp_prob and self.pos >= self.min_variation_start and ( self.last_variation_pos is None or self.last_variation_pos + self.min_variation_dist <= self.pos ):
          new_c = self.add_snp( c )
          result += new_c
          self.last_variation_pos = self.pos
        # insert
        elif random.uniform(0, 1) < self.insert_prob and self.pos >= self.min_variation_start and self.pos > self.max_insert_len and ( self.last_variation_pos is None or self.last_variation_pos + self.min_variation_dist <= self.pos ): # TODO reads can get -ve reference
          new_c = self.add_insertion( c, fragment_pos, len(fragment) )
          self.log( 'added insertion at %i' % self.pos )
          result += new_c + c # insertion gets placed before current base
          self.last_variation_pos = self.pos - 1 # -1 because insertion is placed before current
        # delete
        elif self.pos > 0 and random.uniform(0, 1) < self.delete_prob and self.pos >= self.min_variation_start and ( self.last_variation_pos is None or self.last_variation_pos + self.min_variation_dist <= self.pos ): 
          self.add_deletion( c )
          self.log( 'added deletion at %i' % self.pos )
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
    
class FastaLength(object):
  def __init__(self, reader):
    self.reader = reader
    self.length = 0
    for item in self.reader.items():
       self.length += len(item)

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
      for item in self.reader.items():
         self.fasta += item
      self._length = len(self.fasta)
    return self._length

  def base_at(self, i):
    if i >= len(self.fasta):
      for item in self.reader.items():
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

  def read_to( self, pos ):
    '''
      read up to pos; return false if unable
    '''
    if pos >= len(self.fasta):
      for item in self.reader.items():
        self.fasta += item
        if pos < len(self.fasta): # success
          return True
      return False # failed
    else:
      return True

class FastaReaderFromVCF(object):
  def __init__(self, fasta, vcf):
    self.reader = fasta
    self.vcf = vcf
    self.pos = 0

  def items(self):
    while True:
      next_fragment = self.next_item()
      if next_fragment is None:
        break
      else:
        yield next_fragment

  def next_item(self):
    fragment = self.reader.next_item()
    if fragment is None:
      return None
    else:
      fragment = list(fragment)
      new_pos = self.pos + len(fragment)
      for i in xrange(self.pos, new_pos):
        if i in self.vcf.snp_map:
          fragment[i - self.pos] = self.vcf.snp_list[self.vcf.snp_map[i]]['alt']
      self.pos = new_pos
      return ''.join(fragment)

class FastaReader(object):
  '''
    yields fragments from a fasta file object
>>> import StringIO
>>> g = StringIO.StringIO( '>\nabc\ndef' )
>>> p = bio.FastaReader( g )
>>> [ f for f in p.items() ]
['abc', 'def']
  '''
  def __init__(self, genome, include_headers=False):
    self.genome = genome
    self.name = None
    self.has_next = True
    self.include_headers = include_headers
    self.future_fragment = self._next_item()

  def items(self):
    while True:
      current_fragment = self.future_fragment
      if current_fragment is None:
        break
      self.future_fragment = self._next_item()
      self.has_next = self.future_fragment is not None
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
      if line.startswith( '>' ):
        self.name = line[1:].strip()
      if line.startswith( '>' ) and not self.include_headers:
          pass
      else:
        return line.strip()
    return None

  def has_next_item(self):
    return self.has_next

  def __repr__(self):
    return "name: %s" % self.name

class FastaStats(object):
  '''calculate some overall stats for a fasta file'''
  def __init__( self, fasta, read_length=100, log=bio.log_stderr ):
    '''
      @fasta: file handle
    '''
    self.stats = { 'count': 0, 'gc': [], 'entropy': [] }
    current = ''
    lines = 0
    for line in fasta:
      line = line.strip()
      if not line.startswith( '>' ):
        self.stats['count'] += len(line) # total fasta length
        current += line
        if len(current) >= read_length:
          add = len(current) - read_length
          for idx in xrange(0, add):
            feature = features.ReadFeature( current[idx:idx+read_length] )
            self.stats['gc'].append( feature.gc() )
            self.stats['entropy'].append( feature.entropy() )
          current = current[add:]
      lines += 1
      if lines % 1000 == 0:
        log( '%i lines processed' % lines )
    #log( self.stats )

class MultiFastaReaderContainer(object):
  def __init__(self, genome):
    reader = MultiFastaReader( genome )
    self.fastas = {}
    for item in reader.items():
      self.fastas[item.name] = item

  def find_chromosome( self, name ):
    if name in self.fastas:
      return self.fastas[name]
    for key in self.fastas:
      if key.split( ' ' )[0] == name:
        return self.fastas[key]
    return None

class MultiFastaReader(object):
  '''
    use items to iterate over a list of FastaReader objects
  '''
  def __init__(self, genome):
    '''
      genome is a file like object
    '''
    self.genome = genome

  def items( self ):
    '''
      iterate over fasta sequences
    '''
    current = []
    found_header = False
    for line in self.genome:
      if line.startswith( '>' ):
        if found_header:
          yield FastaReader( genome=StringIO.StringIO( ''.join(current) ) ) # TODO inefficient
          current = [ line ]
        else:
          found_header = True
          current.append( line )
      else:
        current.append( line )
    yield FastaReader( genome = StringIO.StringIO( ''.join(current) ) )

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
    '''snv errors with probability error_prob'''
    transitions = { 'A': 'TGC', 'T': 'GCA', 'G': 'ACT', 'C': 'AGT', 'N': 'N' }
    def uniform_error_profile( bp ):
      if random.random() < error_prob:
        return transitions[bp][random.randint(0, len(transitions[bp]) - 1)]
      else:
        return bp
    return uniform_error_profile

  @staticmethod
  def create_homopolymer_error_profile( error_prob, error_length ):
    possibles = 'ACGT' 
    def homopolymer_error_profile( bp ):
      if random.random() < error_prob:
        return ''.join( ( bp, possibles[random.randint(0, len(possibles)-1)] * error_length ) )
      else:
        return bp
    return homopolymer_error_profile

class RepeatedMultiFastaGenerator( object ):
  def __init__( self, multi_reader, out_fh, multiplier, cfg ):
    for reader in multi_reader.items():
      RepeatedFastaGenerator( reader, out_fh, multiplier, cfg )

class RepeatedFastaGenerator( object ):
  def __init__( self, reader, out_fh, multiplier, cfg, log=bio.log_quiet ):
    '''
      this just makes repeats of an entire fasta with a specified snp mutation prob
      @reader: FastaReader
      @cfg: config dictionary
    '''
    # reader = bio.FastaReader( open( sys.argv[2], 'r' ) )
    out_fh.write( '>generated fasta %ix from %s\n' % ( multiplier, reader.name ) )
    sequence = StringIO.StringIO()
    for item in reader.items():
      if len(item) > 0:
        out_fh.write( '%s\n' % item ) # write unmodified
        sequence.write( '%s\n' % item ) # remember

    for i in xrange(0, multiplier - 1):
      sequence.seek(0)
      reader = FastaReader( sequence )
      mutator = FastaMutate( reader, snp_prob=float(cfg['mult_snp_prob']), insert_prob=0, delete_prob=0, log=log )
      for item in mutator.items():
        if len(item) > 0:
          out_fh.write( '%s\n' % item ) # write unmodified
       #sys.stdout.write( sequence )

BASES = 'ACGT' 

class SequenceGenerator( object ):
  def __init__( self, length, probs=(0.25, 0.25, 0.25, 0.25) ):
    result = []
    for _ in xrange( length ):
      x = random.random()
      for prob_idx in xrange(len(probs)):
        if x < sum(probs[:prob_idx+1]):
          result.append(BASES[prob_idx])
          break
    self.sequence = ''.join(result)  

  def mutate( self, position ):
    current_base_index = BASES.index(self.sequence[position])
    next_base_index = ( current_base_index + 1 ) % len(BASES)
    new_base = BASES[ next_base_index ]
    self.sequence = self.sequence[:position] + new_base + self.sequence[position+1:]

if __name__ == "__main__":
  pass
