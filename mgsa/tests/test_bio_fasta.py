
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

  def test_consensus_snp( self ):
    p = bio.ProbabilisticFasta()
    p.add( 'AA', 0, 1. ) # AA-
    p.add( 'CC', 1, 1. ) # -CC
    p.add( 'CG', 1, 1. ) # -CG
    c = p.consensus_at( 0 )
    self.assertEqual( 1, c[0] )
    self.assertEqual( 'A', c[1] )
    self.assertEqual( 1, c[2] )
    self.assertEqual( 1, c[3] )
    c = p.consensus_at( 1 )
    self.assertEqual( 'C', c[1] )
    self.assertEqual( 2, c[2] )
    self.assertEqual( 3, c[3] )
    c = p.consensus_at( 2 )
    self.assertEqual( 'C', c[1] )
    self.assertEqual( 1, c[2] )
    self.assertEqual( 2, c[3] )

  def test_simple_count( self ):
    p = bio.ProbabilisticFasta( log = bio.log_quiet )
    p.add( 'TTATT', 0, 1. )
    self.assertEqual( 1.0, p.count(0)['T'] ) # 
    self.assertEqual( 0.0, p.count(0)['A'] )

  def test_simple( self ):
    p = bio.ProbabilisticFasta( log = bio.log_quiet )
    p.add( 'TTATT', 0, 1. )
    # num = base count + 1 => 2
    # denom = 1 (total count) + 2 (#base types) => 3
    self.assertEqual( 2./3, p.confidence(0)['T'] ) 
    self.assertEqual( 2./3, p.confidence(1)['T'] )
    self.assertEqual( 2./3, p.confidence(2)['A'] )
    self.assertEqual( 2./3, p.confidence(3)['T'] )
    self.assertEqual( 2./3, p.confidence(4)['T'] )

  def test_insert( self ):
    p = bio.ProbabilisticFasta( log = bio.log_quiet )
    # ref = TTATT
    # read = TTAATT
    # cigar = 3M1I2M
    p.add( 'TTA', 0, 1. ) # 2M
    p.insert( 2, 1 ) # 1I
    p.add( 'TT', 3 )
    self.assertEqual( [1., 1., 0, 1., 1.], p.genome['T'] )
    self.assertEqual( [0., 0., 1.], p.genome['A'] )

  def test_delete( self ):
    p = bio.ProbabilisticFasta( log = bio.log_quiet )
    # ref = TTAATT
    # read = TTATT
    # cigar = 3M1D2M
    p.add( 'TTA', 0, 1. ) # 2M
    p.delete( 3, 1 ) # 1D
    p.add( 'TT', 4 )
    self.assertEqual( [1., 1., 0, 0., 1., 1.], p.genome['T'] )
    self.assertEqual( [0., 0., 1.], p.genome['A'] )
    self.assertEqual( [0., 0., 0., 1.], p.genome['-'] )

class TestFastaMutate( unittest.TestCase ):
  def setUp(self):
    pass

  def _make_reader( self, l ):
    return bio.FastaReader( StringIO.StringIO( '\n'.join( l ) ) )

  def test_no_mutate(self):
    #print "test_no_mutate: start"
    reader = self._make_reader( ( 'gg', 'cc' ) )
    mutator = bio.FastaMutate( reader, log=self.quiet_log, snp_prob=0, insert_prob=0, delete_prob=0 )
    items = []
    for item in mutator.items():
      #print "appending", item
      items.append(item)
    #print "items", items
    self.assertTrue( 2, len(items) )
    self.assertEqual( 'gg', items[0] )
    self.assertEqual( 'cc', items[1] )
    #print "test_no_mutate: end"

  def test_snp(self):
    #print "test_snp: start"
    reader = self._make_reader( ( 'gg', 'cc' ) )
    mutator = bio.FastaMutate( reader, log=self.quiet_log, snp_prob=1, insert_prob=0, delete_prob=0, allow_end_mutate=True )
    items = []
    for item in mutator.items():
      #print "appending", item
      items.append(item)
    #print "items", items
    self.assertTrue( 2, len(items) )
    self.assertTrue( 2, len(items[0]) )
    self.assertTrue( 2, len(items[1]) )
    self.assertNotEqual( 'g', items[0][0] )
    self.assertNotEqual( 'g', items[0][1] )
    self.assertNotEqual( 'c', items[1][0] )
    self.assertNotEqual( 'c', items[1][1] )

  def test_insert(self):
    #print "test_snp: start"
    reader = self._make_reader( ( 'gg', 'cc' ) )
    mutator = bio.FastaMutate( reader, log=self.quiet_log, snp_prob=0, insert_prob=1, delete_prob=0, allow_end_mutate=True )
    items = []
    for item in mutator.items():
      items.append(item)
    self.assertTrue( 2, len(items) )
    self.assertTrue( 2, len(items[0]) )
    self.assertTrue( 2, len(items[1]) )
    # 1st two don't have insertions
    self.assertEqual( 'gg', items[0] )
    # 2nd two do have insertions
    self.assertEqual( 4, len(items[1]) )
    self.assertEqual( 'c', items[1][1] )
    self.assertEqual( 'c', items[1][3] )

  def quiet_log(self, msg):
    pass

class TestFastaReader( unittest.TestCase ):
  def setUp(self):
    pass

  def test_items( self ):
    g = StringIO.StringIO( '>\nabc\ndef' )
    p = bio.FastaReader( g )
    items = [ f for f in p.items() ]
    self.assertEqual( ['abc', 'def'], items )

  def test_items( self ):
    g = StringIO.StringIO( '>\nabc\ndef' )
    p = bio.FastaReader( g )
    count = 0
    for item in p.items():
      if count == 0:
        self.assertEqual( 'abc', item )
        self.assertTrue( p.has_next_item() )
      if count == 1:
        self.assertEqual( 'def', item )
        self.assertFalse( p.has_next_item() )
      count += 1
    self.assertEqual( 2, count )

class TestMultiFastaReader( unittest.TestCase ):
  def test_simple( self ):
    mg = StringIO.StringIO( '>s1\npqr\n>s2\nstu\nvwx' )
    mfg = bio.MultiFastaReader( mg )
    count = 0
    for fr in mfg.items():
      if count == 0:
        self.assertEqual( 's1', fr.name )
        items = [ f for f in fr.items() ]
        self.assertEqual( [ 'pqr' ], items )
      if count == 1:
        self.assertEqual( 's2', fr.name )
        items = [ f for f in fr.items() ]
        self.assertEqual( [ 'stu', 'vwx' ], items )
      count += 1
    self.assertEqual( 2, count )

class TestRepeatedFastaGenerator( unittest.TestCase ):
  def test_multi( self ):
    mg = StringIO.StringIO( '>s1\npqr\n>s2\nstu\nvwx' )
    mfg = bio.MultiFastaReader( mg )
    og = StringIO.StringIO()
    cfg = { 'mult_snp_prob': 0. }
    bio.RepeatedMultiFastaGenerator( mfg, og, 2, cfg )
    self.assertEqual( '>generated fasta 2x from s1\npqr\npqr\n>generated fasta 2x from s2\nstu\nvwx\nstu\nvwx\n', og.getvalue() )

if __name__ == '__main__':
  unittest.main()
