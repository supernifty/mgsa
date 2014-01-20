import threading
import datetime
import time

class MultipleGenomeSequenceAlignerThread(threading.Thread):
  '''
    Threaded wrapper to MultipleGenomeSequenceAligner
  '''

  def __init__(self, genomes, sequences):
    threading.Thread.__init__(self)
    self.worker = MultipleGenomeSequenceAligner(genomes, sequences, self)
    self.status = []
    self.running = True
    self.start()

  def run(self):
    self.worker.align()
    self.running = False

  def stop_align(self):
    self.running = False

  def update( self, completion, status ):
    self.status.append( { 'when': datetime.datetime.now().strftime("%Y%m%d-%H%M%S"), 'what': status, 'complete': completion } )
    return self.running

  def clear_status( self ):
    result = self.status[:]
    self.status = []
    return result

class InterruptedError(Exception):
  pass

class MultipleGenomeSequenceAligner(object):
  '''
    Main class that accepts multiple genomes and attempts to align a sequence to them
  '''

  def __init__(self, genomes, sequences, status):
    self.genomes = [ FastaReader(genome) for genome in genomes ]
    self.sequences = PairedFastqReader(sequences)
    self.status = status

  def _update_status(self, completed, what):
    if self.status is not None:
      result = self.status.update( completed, what )
      if not result:
        raise InterruptedError

  def align(self):
    try:
      self._update_status( 0, "starting..." ) 
      time.sleep(5)
      self._update_status( 100, "finished" )
    except InterruptedError:
      pass

class Aligner(object):
  def __init__(self, genomes):
    self.genomes = genomes
    self.mapped_fragments = []
  
  def map_fragments( self, fragments ):
    '''
      input is an array of paired fragments
    '''
    for fragment in fragments:
      # TODO map this fragment!
      fragment.mapped = True

  def write_sam( self, target ):
    writer = SamGenerator(target, alignment)
    writer.write()

class Fragment(object):
  def __init__(self, sequence):
    self.sequence = sequence
    self.mapped = False

  def __repr__(self):
    return self.sequence

class SamGenerator(object):
  '''
    writes out an alignment in SAM format
  '''
  def __init__(self, target, alignment):
    self.target = target
    self.alignment = alignment

  def _write_header(self):
    pass

  def write(self):
    self._write_header()
    for fragment in self.alignment.mapped_fragments:
      target.write( '%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\n' % ( 
        '*', # qname SimSeq_2ab57
        0, # flag 99  
        '*', # rname gi|556503834|ref|NC_000913.3|
        0, # pos 14 
        255, # mapq 44 
        '*', # cigar  100M 
        '*', # rnext  = 
        0, # pnext 245 
        0, # tlen 331 
        '*', # seq GACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGTTTACCTGCCGTGAGTAAATTAAACTTTTATT 
        '*' ) # qual GDAGGFGG?GEGGAEG=F>FGGFDBG8;GAGFEGGFFG?A:DEEDAGDCEGDFCD@=AE?BB:EFA=C>,CDBECEEDF#CD=F=DF.DFF#GC=AB#:>
        # optional AS:i:189  XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:67G23A8  YS:i:198  YT:Z:CP
      )

class PairedFastqReader(object):
  '''
    yields fragments from a fastq file object
>>> import StringIO
>>> fq = StringIO.StringIO( '@\nabc\n@\ndef' )
>>> fq2 = StringIO.StringIO( '@\nab2\n@\nde2' )
>>> p = mgsa.PairedFastqReader( [ fq, fq2 ] )
>>> [ f for f in p.items() ]
[['abc', 'ab2'], ['def', 'de2']]
  '''
  
  def __init__(self, sequences):
    self.sequences = sequences

  def items(self):
    '''
      gets the next fragment from each sequence file and returns as an array of fragments
    '''
    while True:
      fragments = self.next_item()
      if fragments[0] is None:
        break
      yield fragments

  def next_item(self):
    fragments = []
    for sequence in self.sequences:
      fragments.append( self._next(sequence) )
    return fragments

  def _next(self, sequence):
    for line in sequence.readline():
      if line.startswith( '@' ):
        fragment = Fragment( sequence.readline().strip() )
        return fragment
    return None

class FastaReader(object):
  '''
    yields fragments from a fasta file object
>>> import StringIO
>>> g = StringIO.StringIO( '>\nabc\ndef' )
>>> p = mgsa.FastaReader( g )
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


if __name__ == "__main__":
  #import doctest
  #doctest.testmod()
  pass
