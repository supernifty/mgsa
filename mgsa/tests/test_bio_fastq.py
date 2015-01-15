
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


if __name__ == '__main__':
  unittest.main()
