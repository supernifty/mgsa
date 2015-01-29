
import unittest
import StringIO

import bio

class TestFastqGenerator( unittest.TestCase ):

  def setUp(self):
    pass

  def test_simple(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nTTTTTTTTTT' )
    cfg = bio.Config().read_config_line( 'coverage 2, read_length 8' )
    target = StringIO.StringIO()
    fq = bio.FastqGenerator( fasta, cfg, vcf=None, target_fh=target, log=self.quiet_log )
    self.assertEqual( '@mgsa_seq_3~0~0\nAAAAAAAT\n+\n~~~~~~~~\n@mgsa_seq_7~0~0\nAAATTTTT\n+\n~~~~~~~~\n@mgsa_seq_11~0~0\nTTTTTTTT\n+\n~~~~~~~~\n', target.getvalue() )

  def log(self, msg):
    print msg

  def quiet_log(self, msg):
    pass

class TestFastqPosGenerator( unittest.TestCase ):

  def test_simple(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nGTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=10, read_length=5 ) # pos 10 is the first T
    self.assertEqual( 5, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 21, len(lines) )
    self.assertEqual( 'AAAAG', lines[1] )
    self.assertEqual( 'AAAGT', lines[5] )
    self.assertEqual( 'AAGTT', lines[9] )
    self.assertEqual( 'AGTTT', lines[13] )
    self.assertEqual( 'GTTTT', lines[17] )
    self.assertEqual( 20, len(fq) )

  def test_end(self):
    fasta = StringIO.StringIO( '>comment\nAAAAA\nTTTTG' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=9, read_length=5 ) # pos 10 is the first T
    self.assertEqual( 1, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 5, len(lines) )
    self.assertEqual( 'TTTTG', lines[1] )
    self.assertEqual( 10, len(fq) )
 
  def test_snp(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nTTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=10, read_length=5, variation='snp' ) # pos 10 is the first T
    self.assertEqual( 5, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 21, len(lines) )
    self.assertEqual( 'AAAAA', lines[1] )
    self.assertEqual( 'AAAAT', lines[5] )
    self.assertEqual( 'AAATT', lines[9] )
    self.assertEqual( 'AATTT', lines[13] )
    self.assertEqual( 'ATTTT', lines[17] )
    self.assertEqual( 20, len(fq) )

if __name__ == '__main__':
  unittest.main()
