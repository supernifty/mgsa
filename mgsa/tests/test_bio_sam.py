
import unittest
import StringIO

import bio

class TestSamToVCF( unittest.TestCase ):

  def setUp(self):
    pass

  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AAAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, self.quiet_log )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )

  def test_snp(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, self.quiet_log )
    self.assertEqual( 1, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )

  def test_insert(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      2M2I15M       *       0       0       AACCAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, self.quiet_log )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )

  def test_delete(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_0~0~0  0       generated       1      60      2M2D18M       *       0       0       AAAAAAAAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AACCAAAAAAAA\nTTTTTTTTTT' ) # 10A, 10T
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, self.quiet_log )
    #print target_vcf
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )
    self.assertEqual( 'ACC', target_vcf.manager.indel_list[0].before )
    self.assertEqual( 'A', target_vcf.manager.indel_list[0].after )
    
  def log(self, msg):
    print msg

  def quiet_log(self, msg):
    pass

class TestBamReader( unittest.TestCase ):
  def setUp(self):
    pass

