
import unittest
import StringIO

import bio

class TestFastaDiff( unittest.TestCase ):

  def setUp(self):
    pass

  def _make_reader( self, l ):
    return bio.FastaReader( StringIO.StringIO( '\n'.join( l ) ) )

  def _make_fasta( self, l, log=None ):
    if log is None:
      log = self.quiet_log
    result = bio.ProbabilisticFasta(log=log)
    for fragment in l:
      result.add( fragment, 1 )
    return result

  def test_perfect(self):
    diff = bio.FastaDiff( 
      self._make_reader( ( 'gg', 'cc' ) ), 
      self._make_fasta( ( 'ggcc', ) ), 
      self.quiet_log )
    self.assertEqual( 0, diff.error_total )
    self.assertEqual( 5, diff.candidate_pos )

  def test_snp_error(self):
    diff = bio.FastaDiff( 
      self._make_reader( ( 'gg', 'cc' ) ), 
      self._make_fasta( ( 'gtcc', ) ), 
      self.quiet_log )
    self.assertEqual( 1, diff.error_total )
    self.assertEqual( 5, diff.candidate_pos )

  def test_single_insert_ok(self):
    candidate = self._make_fasta( ( 'ggccaa', ) ) 
    candidate.insert( 't', 3 )
    diff = bio.FastaDiff( 
      self._make_reader( ( 'gg', 'tcc', 'aa' ) ), 
      candidate,
      self.quiet_log )
    #print candidate.insertions
    self.assertEqual( 0, diff.error_total )
    self.assertEqual( 7, diff.candidate_pos )

  def test_consecutive_insert_ok(self):
    candidate = self._make_fasta( ( 'gggatgcg', ), log=self.quiet_log ) 
    candidate.insert( 't', 4 )
    candidate.insert( 'a', 5 )
    # correct is gggtaatggg
    diff = bio.FastaDiff( 
      self._make_reader( ( 'ggg', 'taat', 'gcg' ) ), 
      candidate,
      self.quiet_log )
    #print candidate.insertions
    self.assertEqual( 0, diff.error_total )
    self.assertEqual( 9, diff.candidate_pos )

  def test_contig_insert_ok(self):
    candidate = self._make_fasta( ( 'gggatgcg', ), log=self.quiet_log ) 
    candidate.insert( 'ta', 4 )
    # correct is gggtaatggg
    diff = bio.FastaDiff( 
      self._make_reader( ( 'ggg', 'taat', 'gcg' ) ), 
      candidate,
      self.quiet_log )
    #print candidate.insertions
    self.assertEqual( 0, diff.error_total )
    self.assertEqual( 9, diff.candidate_pos )

  def test_ambiguous_insert_ok(self):
    candidate = self._make_fasta( ( 'cagt', ), log=self.quiet_log ) 
    candidate.insert( 'a', 2 )
    # correct is caagt
    diff = bio.FastaDiff( 
      self._make_reader( ( 'ca', 'agt' ) ), 
      candidate,
      self.quiet_log )
    #print candidate.insertions
    self.assertEqual( 0, diff.error_total )
    self.assertEqual( 5, diff.candidate_pos )

  def test_wrong_insert(self):
    #print "enter test_insert_miss"
    candidate = self._make_fasta( ( 'cagt', ), log=self.quiet_log ) 
    candidate.insert( 'g', 2 )
    # correct is caagt
    reference = self._make_reader( ( 'ca', 'agt' ) )
    diff = bio.FastaDiff( 
      reference,
      candidate,
      self.quiet_log )
    #print "insertions:", candidate.insertions
    self.assertEqual( 1, diff.error_total )
    self.assertEqual( 5, diff.candidate_pos )
    #print "exit test_insert_miss"

  def log(self, msg):
    print msg

  def quiet_log(self, msg):
    pass

class TestProbabilisticFasta( unittest.TestCase ):
  def setUp(self):
    pass

  def test_confidence( self ):
    p = bio.ProbabilisticFasta()
    p.add( 'AA', 0, 0.5 )
    p.add( 'CC', 1, 0.5 )
    self.assertEqual( 0.6, p.confidence(0)['A'] )
    self.assertEqual( 0.4, p.confidence(0)['C'] )
    self.assertEqual( 0.5, p.confidence(1)['A'] )
    self.assertEqual( 0.5, p.confidence(1)['C'] )
    self.assertEqual( 0.6, p.confidence(2)['C'] )

if __name__ == '__main__':
  unittest.main()
