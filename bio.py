
import datetime
import random
import sys

class ProbabilisticFasta(object):
  '''
    generate fasta with probabilities attached from mapped fragments
  '''
  
  def __init__( self ):
    self.genome = {} # counts over position
    self.prior = {} # overall counts
    self.total = 0 # total # of nucleotides
    self.size = 0 # length of genome

  def add( self, fragment, start ):
    for i in xrange(0, len(fragment) ):
      value = fragment[i]
      position = start + i
      if value not in self.genome:
        self.genome[value] = []
        self.prior[value] = 0
      if len(self.genome[value]) < position + 1:
        self.genome[value].extend( [0] * ( position - len(self.genome[value]) + 1 ) )
        if len(self.genome[value]) > self.size:
          self.size = len(self.genome[value])
      self.genome[value][position] += 1
      self.prior[value] += 1
      self.total += 1

  def count( self, position ):
    '''
      counts of nucleotides seen at a position
    '''
    result = {}
    for value in self.genome:
      if len(self.genome[value]) > position:
        result[value] = self.genome[value][position]
    return result

  def confidence( self, position ):
    '''
      returns probabilities of nucleotides
      TODO priors
    '''
    counts = self.count( position )
    denominator = sum( counts[k] for k in counts ) + len(self.prior)
    probabilities = {}
    for count in counts:
      probabilities[count] = ( counts[count] + 1.0 ) / denominator
    return probabilities

  def consensus( self, start=1, end=-1 ):
    if end == -1 or end > self.size:
      end = self.size
    result = ''
    for i in xrange(start, end):
      best = '0'
      best_value = 0
      for value in self.genome:
        if len(self.genome[value]) > i and self.genome[value][i] > best_value:
          best_value = self.genome[value][i]
          best = value
      result += best
    return result 

class FastaMutate(object):
  '''
    change a reference fasta
  '''
  probabilities = 'AAACCTTGGG'

  def __init__( self, reader, log, snp_probability=0.01 ):
    self.reader = reader
    self.snp_prob = snp_probability
    self.mutations = 0
  
  def items(self):
    while True:
      fragment = self.reader.next_item()
      if fragment is None:
        break
      # apply mutations
      fragment = self.mutate( fragment )
      yield fragment

  def mutate(self, fragment):
    result = ''
    for c in fragment:
      if random.uniform(0, 1) < self.snp_prob:
        new_c = self.probabilities[random.randint(0, len(self.probabilities)-1)]
        if new_c != c:
          self.mutations += 1
        result += new_c
      else:
        result += c
    return result

class FastaDiff(object):
  '''
    compare fastas
  '''

  def __init__( self, reader, candidate, log ):
    pos = 1
    last_log_pos = 0
    self.errors = {}
    self.error_total = 0
    self.candidate = candidate
    for item in reader.items():
      candidate_item = candidate.consensus( pos, pos + len(item) )
      if candidate_item != item:
        # find differences
        self.find_differences( item, candidate_item, pos )
      pos += len(item)
      if pos < 1000 and pos - last_log_pos >= 100 or pos - last_log_pos >= 10000:
        log( 'processed %i items - %i errors' % ( pos, self.error_total ) )
        last_log_pos = pos

  def find_differences( self, item, candidate, start ):
    for i in xrange(0, min(len(item), len(candidate))):
      if item[i] != candidate[i]:
        self.error_total += 1
        key = '%s->%s' % ( item[i], candidate[i] )
        log( 'error: %s at %i: counts %s' % ( key, start + i, self.candidate.count( start + i ) ) )
        if key not in self.errors:
          self.errors[key] = 0
        self.errors[key] += 1

    if len(item) != len(candidate):
      key = 'length'
      if key not in self.errors:
        self.errors[key] = 0
      self.errors[key] += 1
    
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

  def items(self):
    while True:
      fragment = self.next_item()
      if fragment is None:
        break
      yield fragment

  def next_item(self):
    for line in self.genome:
      if not line.startswith( '>' ):
        return line.strip()
    return None

class SamToFasta(object):
  '''
    read and evaluate assembler data
  '''  

  def __init__( self, sam, log ):
    self.log = log
    self.fasta = ProbabilisticFasta()
    self.stats = { 'mapped': 0, 'unmapped': 0, 'unknown_mapping': 0, 'lines': 0 }
    for line in sam:
      self.parse_line( line.strip() )
      self.stats['lines'] += 1
      if self.stats['lines'] < 10 or self.stats['lines'] < 1000 and self.stats['lines'] % 100 == 0 or self.stats['lines'] % 10000 == 0:
        self.log( '%i lines processed' % self.stats['lines'] )

  def parse_line( self, line ):
    #0 SimSeq_1101
    #1 163
    #2 gi|205355313|ref|NZ_ABGQ01000001.1|
    #3 1
    #4 44
    #5 100M
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
      if flag & 0x04 != 0:
        self.stats['unmapped'] += 1
      elif flag & 0x02 != 0:
        self.stats['mapped'] += 1
        pos = int(fields[3])
        self.fasta.add( fields[9], pos )
      else:
        self.stats['unknown_mapping'] += 1

def log(msg):
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  sys.stderr.write( '%s: %s\n' % ( when, msg ) )

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

  m = FastaMutate( FastaReader( open( '../data/test/test1k.fasta', 'r' ) ), log )
  o = open( '../data/test/test1ksnp.fasta', 'w' )
  o.write( '>mutated-snp-test1k\n')
  for line in m.items():
    o.write( '%s\n' % line )
  o.close()
  log( m.mutations )
