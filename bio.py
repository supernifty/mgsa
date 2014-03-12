
import datetime
import random
import re
import sys

DELETE_BASE = '-'

class ProbabilisticFasta(object):
  '''
    generate fasta with probabilities attached from mapped fragments
  '''
  
  def __init__( self, log=None ):
    self.genome = {} # counts over position { 'A': [ 1, 0, 2, ... ], 'G': [ ... ] }
    self.insertions = {} # { 0: { 'AAA': 1, 'BB': 2 } }
    self.inserted = {}
    self.deleted = set()
    self.prior = {} # overall counts
    self.total = 0 # total # of nucleotides
    self.size = 0 # length of genome
    self.log = log

  def add( self, fragment, start, confidence=1.0 ):
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
      self.genome[value][position] += confidence
      self.prior[value] += 1
      self.total += 1

  def insert( self, fragment, start, confidence=1.0 ):
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

  def consensus_at( self, i=1 ):
    best = '0'
    best_value = 0
    coverage = 0
    result = ''
    move = 0
    for value in self.genome:
      if len(self.genome[value]) > i:
        coverage += self.genome[value][i]
        if self.genome[value][i] > best_value:
          best_value = self.genome[value][i]
          best = value
    # check for insert 
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
        result += best_insertion
        self.inserted[i] = best_insertion
        self.log( 'included insertion %s at %i' % (best_insertion, i) )
        move += len(best_insertion)
      else:
        pass#self.log( 'skipped insertion at %i with val %f with noinsert %f' % ( i, best_insertion_value, coverage - insertion_coverage ) )
    # include base if not deleted
    if best == DELETE_BASE:
      self.deleted.add(i)
      self.log( 'deletion at %i' % (i) )
    else:
      result += best
      move += 1
     
    return (move, result)
  
  def consensus_count( self, start=1, count=1 ):
    '''
      get the next count characters from start
    '''
    result = ''
    total_move = 0
    reference_move = 0
    while len(result) < count and start + reference_move < self.size:
      move, s = self.consensus_at( start + reference_move )
      result += s 
      reference_move += 1
    #self.log( 'consensus for %i count %i is %s (%i)' % (start, count, result, len(result)) )
    return (reference_move, result)

  def consensus( self, start=1, end=-1 ):
    '''
      return the corresponding string for the segment specified.
      returned segment can be longer or shorter than end-start - segment applies to the reference.
    '''
    if end == -1 or end > self.size:
      end = self.size
    result = ''
    self.insertion_count = 0
    for i in xrange(start, end):
      move, add = self.consensus_at( i )
      result += add
    #self.log( 'consensus for %i to %i is %s (%i)' % (start, end, result, len(result)) )
    return result 

class FastaMutate(object):
  '''
    change a reference fasta
  '''
  probabilities = 'AAACCTTGGG'

  def __init__( self, reader, log, vcf_file=None, snp_prob=0.01, insert_prob=0.01, delete_prob=0.01 ):
    self.reader = reader
    self.snp_prob = snp_prob
    self.insert_prob = insert_prob
    self.delete_prob = delete_prob
    self.mutations = 0
    self.vcf_file = vcf_file
    if vcf_file is not None:
      self.vcf = VCF( writer=VCFWriter(vcf_file) )
    else:
      self.vcf = VCF()
    self.pos = 0

    seed = random.randint(0, sys.maxint)
    random.seed(seed)
    if log is not None:
      log( 'seed: %i' % seed )
  
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
      # snp
      if random.uniform(0, 1) < self.snp_prob:
        new_c = self.probabilities[random.randint(0, len(self.probabilities)-1)]
        if new_c != c:
          self.mutations += 1
          if self.vcf is not None:
            self.vcf.snp( self.pos, c, new_c )
        result += new_c
      # insert
      elif random.uniform(0, 1) < self.insert_prob:
        new_c = self.probabilities[random.randint(0, len(self.probabilities)-1)]
        self.mutations += 1
        if self.vcf is not None:
          self.vcf.indel( self.pos, c, new_c + c )
        result += new_c + c
      
      elif self.pos > 0 and random.uniform(0, 1) < self.delete_prob: # delete
        self.mutations += 1
        if self.vcf is not None:
          self.vcf.indel( self.pos - 1, self.previous + c, c )

      else: # no mutation
        result += c

      self.pos += 1
      self.previous = c
    return result

class VCF(object):
  '''
    manage a VCF set
  '''

  def __init__( self, reader=None, writer=None ):
    '''
      @writer: VCFWriter object
    '''
    self.snp_list = []
    self.indel_list = []
    self.writer = writer
    if reader is not None:
      self.load( reader )

  def load( self, reader ):
    for line in reader:
      fields = line.split()
      if len(fields) > 5:
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        if len(ref) == 1 and len(alt) == 1:
          self.snp_list.append( { 'pos': int(pos), 'ref': ref, 'alt': alt } )
        else:
          self.indel_list.append( { 'pos': int(pos), 'before': ref, 'after': alt } )

  def snp( self, pos, ref, alt ):
    self.snp_list.append( { 'pos': pos, 'ref': ref, 'alt': alt } )
    if self.writer is not None:
      self.writer.snp( pos, ref, alt )

  def indel( self, pos, before, after ):
    self.indel_list.append( { 'pos': pos, 'before': before, 'after': after } )
    if self.writer is not None:
      self.writer.indel( pos, before, after )

  def write_all( self, writer ):
    vcf_writer = VCFWriter(writer)
    for snp in self.snp_list:
      vcf_writer.snp( snp['pos'], snp['ref'], snp['alt'] )
    for indel in self.indel_list:
      vcf_writer.indel( indel['pos'], indel['before'], indel['after'] )

class VCFWriter(object):
  '''
    write a vcf file
  '''

  def __init__( self, writer ):
    self.writer = writer
    self.writer.write( '''##fileformat=VCFv4.1
##fileDate=%s
##source=SuperniftyMutationSimulator
''' % datetime.datetime.now().strftime("%Y%m%d") )

  def snp( self, pos, ref, alt ):
    self.writer.write( '%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        # chrom, pos, id, ref, alt, qual, filter, info
        '.', pos, '.', ref, alt, '.', 'PASS', 'DP=100'
      ) )

  def indel( self, pos, before, after ):
    self.writer.write( '%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        # chrom, pos, id, ref, alt, qual, filter, info
        '.', pos, '.', before, after, '.', 'PASS', 'DP=100'
      ) )

class VCFDiff(object):
  '''
    compare variations, calculate confusion matrix
  '''
  def __init__( self, vcf, fasta, log ):
    '''
      @vcf: VCF object
      @fasta: ProbabilisticFasta object
    '''
    self.snp_stats = { 'tp': 0, 'fn': 0, 'other': 0 }
    for true_snp in vcf.snp_list: # vcf.snp has pos, ref, alt
      pos = true_snp['pos'] + 1
      value = fasta.consensus_at(pos)[1]
      if value == true_snp['ref']:
        self.snp_stats['fn'] += 1
      elif value == true_snp['alt']:
        self.snp_stats['tp'] += 1
      else: # neither snp ref or alt
        log( 'consensus at %i: %s %s %s ref %s alt %s' % ( 
          pos, 
          fasta.consensus_at(pos-1)[1], 
          value, 
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

class FastaDiff(object):
  '''
    compare fastas
  '''

  def __init__( self, reader, candidate, log ):
    pos = 1
    last_log_pos = 0
    self.errors = {}
    self.error_positions = set() # set of incorrect positions
    self.error_total = 0
    self.candidate = candidate
    self.log = log
    for item in reader.items():
      move, candidate_item = candidate.consensus_count( pos, len(item) )
      if candidate_item != item:
        # find differences
        self.find_differences( item, candidate_item, pos )
      pos += move #len(item)
      if pos < 10000 and pos - last_log_pos >= 1000 or pos - last_log_pos >= 10000:
        self.log( 'processed %i items - %i errors' % ( pos, self.error_total ) )
        last_log_pos = pos

  def find_differences( self, item, candidate, start ):
    for i in xrange(0, min(len(item), len(candidate))):
      if item[i] != candidate[i]:
        self.error_total += 1
        key = '%s->%s' % ( item[i], candidate[i] )
        if self.error_total < 10 or self.error_total % 1000 == 0:
          self.log( 'error: actual->predicted %s at %i: counts %s' % ( key, start + i, self.candidate.count( start + i ) ) )
        if key not in self.errors:
          self.errors[key] = 0
        self.errors[key] += 1
        self.error_positions.add( start + i )

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

  def __init__( self, sam, log, allow_indels=True ):
    self.allow_indels = allow_indels
    self.log = log
    self.fasta = ProbabilisticFasta(log)
    self.stats = { 'mapped': 0, 'unmapped': 0, 'unknown_mapping': 0, 'lines': 0 }
    for line in sam:
      self.parse_line( line.strip() )
      self.stats['lines'] += 1
      if self.stats['lines'] < 5000 and self.stats['lines'] % 1000 == 0 or self.stats['lines'] % 10000 == 0:
        self.log( '%i lines processed' % self.stats['lines'] )
    self.log( self.stats )

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
      if flag & 0x04 != 0:
        self.stats['unmapped'] += 1
      elif flag & 0x02 != 0:
        self.stats['mapped'] += 1
        pos = int(fields[3]) # pos in genome
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
            self.fasta.add( fields[9][fragment_pos:fragment_pos+cigar_len], pos + fragment_pos, confidence=0.5 )
            genome_pos += cigar_len
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
      else:
        self.stats['unknown_mapping'] += 1

#def log(msg):
#  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
#  sys.stderr.write( '%s: %s\n' % ( when, msg ) )

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
