
import unittest
import StringIO

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

  def log(self, msg):
    print msg

  def quiet_log(self, msg):
    pass


if __name__ == '__main__':
  unittest.main()
