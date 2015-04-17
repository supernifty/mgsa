
import unittest
import StringIO

import numpy as np

import bio

class TestHelpers( unittest.TestCase ):

  def setUp(self):
    pass

  def test_fasta_stats(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nTTTTTTTTTT' )
    stats = bio.fasta_stats( fasta, log=self.quiet_log )
    self.assertEqual( 20, stats.stats['count'] )

  def test_fasta_set_width(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nTTTTTTTTTT' )
    target = StringIO.StringIO()
    bio.fasta_set_width( fasta, target, 6, log=self.quiet_log )
    lines = target.getvalue().split()
    self.assertEqual( '>comment', lines[0] )
    self.assertEqual( 'AAAAAA', lines[1] )
    self.assertEqual( 'AAAATT', lines[2] )
    self.assertEqual( 'TTTTTT', lines[3] )
    self.assertEqual( 'TT', lines[4] )
    self.assertEqual( 5, len(lines) )

  def test_fasta_edit_insert(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nTTTTTTTTTT' )
    target = StringIO.StringIO()
    bio.fasta_edit( fasta, target, 'i.5.CC' )
    lines = target.getvalue().split()
    self.assertEqual( '>comment', lines[0] )
    self.assertEqual( 'AAAACCAAAAAATTTTTTTTTT', lines[1] )
    self.assertEqual( 2, len(lines) )

  def test_fasta_edit_delete(self):
    fasta = StringIO.StringIO( '>comment\nAATG\nCCTG' )
    target = StringIO.StringIO()
    bio.fasta_edit( fasta, target, 'd.3.4' )
    lines = target.getvalue().split()
    self.assertEqual( '>comment', lines[0] )
    self.assertEqual( 'AATG', lines[1] )
    self.assertEqual( 2, len(lines) )

  def test_fasta_filter_start(self):
    fasta = StringIO.StringIO( '>comment\nAATG\nCCTG>another\nTTTT\n' )
    target = StringIO.StringIO()
    bio.fasta_filter( fasta, target, 'comment:1-4' )
    lines = target.getvalue().split()
    self.assertEqual( ['AATG'], lines )

  def test_fasta_filter_cross(self):
    fasta = StringIO.StringIO( '>comment\nAATG\nCCTG>another\nTTTT\n' )
    target = StringIO.StringIO()
    bio.fasta_filter( fasta, target, 'comment:3-6' )
    lines = target.getvalue().split()
    self.assertEqual( ['TGCC'], lines )

  def test_fasta_filter_next_end(self):
    fasta = StringIO.StringIO( '>comment\nAATG\nCCTG\n>another\nGGTT\n' )
    target = StringIO.StringIO()
    bio.fasta_filter( fasta, target, 'another:3-4' )
    lines = target.getvalue().split()
    self.assertEqual( ['TT'], lines )

  def test_fasta_find( self ):
    fasta = StringIO.StringIO( '>comment\nAATG\nCCTG\n>another\nGGTT\n' )
    self.assertEqual( [2, 6], bio.fasta_find( fasta, 'TG' ) )

  def test_choose( self ):
    self.assertEqual( 1,  bio.choose( 10, 0 ) )
    self.assertEqual( 45, bio.choose( 10, 2 ) )
    self.assertEqual( 45, bio.choose( 10, 8 ) )
    self.assertEqual( 1, bio.choose( 10, 10 ) )

  def test_binomial_prob( self ):
    self.assertAlmostEqual( 0.6057704364907279, bio.binomial_prob( k=0, m=100, n=200 ) )
    self.assertAlmostEqual( 0.30440725451795375, bio.binomial_prob( k=1, m=100, n=200 ) )
    self.assertAlmostEqual( 0.3660323412732292, bio.binomial_prob( k=0, m=100, n=100 ) )
    self.assertAlmostEqual( 0.36972963764972644, bio.binomial_prob( k=1, m=100, n=100 ) )
    self.assertAlmostEqual( 0.367695424771, bio.binomial_prob( k=0, m=1000, n=1000 ) )
    self.assertAlmostEqual( 0.368063488259, bio.binomial_prob( k=1, m=1000, n=1000 ) )
    self.assertAlmostEqual( 0.36786104643297046, bio.binomial_prob( k=0, m=10000, n=10000 ) )
    self.assertAlmostEqual( 0.3678978362165921, bio.binomial_prob( k=1, m=10000, n=10000 ) )
    self.assertAlmostEqual( 0.3678776017682465, bio.binomial_prob( k=0, m=100000, n=100000 ) )
    self.assertAlmostEqual( 0.36788128058105235, bio.binomial_prob( k=1, m=100000, n=100000 ) )
    self.assertAlmostEqual( 0.18394064029052615, bio.binomial_prob( k=2, m=100000, n=100000 ) )

  def test_poisson( self ):
    x = []
    for y in xrange(0, 10000):
      x.append( bio.poisson( 10 ) )
    self.assertEqual( 10, round( np.mean(x), 0 ) )
    self.assertEqual( 10, round( np.var(x), 0 ) )

  def test_logistic( self ):
    self.assertEqual( 0.5, bio.logistic( 0. ) )

  def test_trim_edges( self ):
    self.assertEqual( [ 1, 2, 3 ], bio.trim_edges( [ 0, 0, 1, 2, 3, 0, 0 ] ) )

  def test_run_length( self ):
    self.assertEqual( (3, 4), bio.longest_run( [ 1, 1, 2, 2, 1, 1, 1 ], 1 ) )
    self.assertEqual( (2, 2), bio.longest_run( [ 1, 1, 2, 2, 1, 1, 1 ], 2 ) )

  def test_all_runs( self ):
    self.assertEqual( { 2: 2, 3: 1 }, bio.all_runs( [ 1, 1, 2, 2, 1, 1, 1, 2, 1, 1 ], 1 ) )

  
  def log(self, msg):
    print msg

  def quiet_log(self, msg):
    pass


if __name__ == '__main__':
  unittest.main()
