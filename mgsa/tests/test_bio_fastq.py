
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

  def test_insert(self):
    fasta = StringIO.StringIO( '>comment\nCCCCCCCCCC\nGTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=10, read_length=5, variation='insert 1' ) # pos 10 is the first T
    self.assertEqual( 5, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 21, len(lines) )
    self.assertEqual( 'CCCAG', lines[1] )
    self.assertEqual( 'CCAGT', lines[5] )
    self.assertEqual( 'CAGTT', lines[9] )
    self.assertEqual( 'AGTTT', lines[13] )
    self.assertEqual( 'GTTTT', lines[17] )
    self.assertEqual( 20, len(fq) )

  def test_insert_start(self):
    fasta = StringIO.StringIO( '>comment\nGCCCCCCCCC\nTTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=0, read_length=5, variation='insert 1' ) # pos 10 is the first T
    self.assertEqual( 2, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 9, len(lines) )
    self.assertEqual( 'AGCCC', lines[1] )
    self.assertEqual( 'GCCCC', lines[5] )
    self.assertEqual( 20, len(fq) )

  def test_delete(self):
    fasta = StringIO.StringIO( '>comment\nCCCCCCCACA\nGATTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=9, read_length=5, variation='delete 1' ) # pos 10 is the G
    lines = target.getvalue().split( '\n' )
    #print lines
    self.assertEqual( 4, written )
    self.assertEqual( 17, len(lines) )
    self.assertEqual( 'CCACG', lines[1] )
    self.assertEqual( 'CACGA', lines[5] )
    self.assertEqual( 'ACGAT', lines[9] )
    self.assertEqual( 'CGATT', lines[13] )
    self.assertEqual( 20, len(fq) )

  def test_delete_3(self):
    fasta = StringIO.StringIO( '>comment\nCCACTCCACA\nGATTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=9, read_length=5, variation='delete 3' ) # pos 9 is the A; expect to delete ACA
    lines = target.getvalue().split( '\n' )
    #print lines
    self.assertEqual( 4, written )
    self.assertEqual( 17, len(lines) )
    self.assertEqual( 'CTCCG', lines[1] )
    self.assertEqual( 'TCCGA', lines[5] )
    self.assertEqual( 'CCGAT', lines[9] )
    self.assertEqual( 'CGATT', lines[13] )
    self.assertEqual( 20, len(fq) )

  def test_delete_long(self):
    fasta = StringIO.StringIO( '>comment\nABCDEFGHIJKLMNOPQRSTUVWXYZ' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=15, read_length=15, variation='delete 10' ) # 
    lines = target.getvalue().split( '\n' )
    # ['@mgsa_seq_0~0~0_variation_D5-10', 'ABCDEFQRSTUVWXY', '+', '~~~~~~~~~~~~~~~', '@mgsa_seq_0~1~0_variation_D4-10', 'BCDEFQRSTUVWXYZ', '+', '~~~~~~~~~~~~~~~', '']
    self.assertEqual( 2, written )
    self.assertEqual( 9, len(lines) )
    self.assertEqual( 'ABCDEFQRSTUVWXY', lines[1] )
    self.assertEqual( 'BCDEFQRSTUVWXYZ', lines[5] )
    self.assertEqual( 26, len(fq) )

  def test_snp_error(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nTTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=10, read_length=5, variation='snp', error='snp 1' ) # pos 10 is the first T
    self.assertEqual( 5, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 21, len(lines) )
    self.assertEqual( 1, bio.Distance( 'AAAAA', lines[1] ).hamming() )
    self.assertEqual( 1, bio.Distance( 'AAAAT', lines[5] ).hamming() )
    self.assertEqual( 1, bio.Distance( 'AAATT', lines[9] ).hamming() )
    self.assertEqual( 1, bio.Distance( 'AATTT', lines[13] ).hamming() )
    self.assertEqual( 1, bio.Distance( 'ATTTT', lines[17] ).hamming() )
    self.assertEqual( 20, len(fq) )

  def test_ins_error(self):
    fasta = StringIO.StringIO( '>comment\nCCCCCCCCCC\nTTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=10, read_length=5, variation=None, error='insert 1 1' ) # pos 10 is the first T
    self.assertEqual( 4, written )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 17, len(lines) )
    self.assertEqual( 1, bio.Distance( 'CCCT', lines[1] ).simple_indel() )
    self.assertEqual( 1, bio.Distance( 'CCTT', lines[5] ).simple_indel() )
    self.assertEqual( 1, bio.Distance( 'CTTT', lines[9] ).simple_indel() )
    self.assertEqual( 1, bio.Distance( 'TTTT', lines[13] ).simple_indel() )
    self.assertEqual( 20, len(fq) )

  def test_del_error(self):
    fasta = StringIO.StringIO( '>comment\nCCCCCCCCAC\nGTATTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    target = StringIO.StringIO()
    written = fq.write( target, pos=10, read_length=5, variation=None, error='delete 1 1' ) # pos 10 is the G
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 6, written )
    self.assertEqual( 25, len(lines) )
    self.assertEqual( -1, bio.Distance( 'CCCACG', lines[1] ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'CCACGT', lines[5] ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'CACGTA', lines[9] ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'ACGTAT', lines[13] ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'CGTATT', lines[17] ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'GTATTT', lines[21] ).simple_indel() )
    self.assertEqual( 20, len(fq) )


if __name__ == '__main__':
  unittest.main()
