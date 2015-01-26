
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
    s = bio.SamToVCF( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )

  def test_snp(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 1, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )

  def test_insert(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      2M2I15M       *       0       0       AACCAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )

  def test_delete(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_0~0~0  0       generated       1      60      2M2D18M       *       0       0       AAAAAAAAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AACCAAAAAAAA\nTTTTTTTTTT' ) # 10A, 10T
    target_vcf = bio.VCF()
    s = bio.SamToVCF( sam, reference, target_vcf, bio.log_quiet )
    #print target_vcf
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )
    self.assertEqual( 'ACC', target_vcf.manager.indel_list[0].before )
    self.assertEqual( 'A', target_vcf.manager.indel_list[0].after )
    
  def log(self, msg):
    print msg

class TestSamToFasta( unittest.TestCase ):

  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    s = bio.SamToFasta( sam, log=bio.log_quiet )
    # stats
    self.assertEqual( 1, s.stats['unknown_mapping'] )
    self.assertEqual( 2, s.stats['lines'] )
    # consensus (note sequence starts at 6)
    # move, result, best, coverage
    self.assertEqual( (1, 'A', 1.0, 1.0), s.fasta.consensus_at( 5 ) )
    self.assertEqual( (1, 'A', 1.0, 1.0), s.fasta.consensus_at( 6 ) )
    self.assertEqual( (1, 'C', 1.0, 1.0), s.fasta.consensus_at( 7 ) )

class TestFastaMutate( unittest.TestCase ):

  def test_simple(self):
    l = ( 'AAAA', 'CCCC')
    fr = bio.FastaReader( StringIO.StringIO( '\n'.join( l ) ) )
    mutator = bio.FastaMutate( fr, snp_prob=1., log=bio.log_quiet )
    result = []
    for fragment in mutator.items():
      result.append( fragment )
    self.assertEqual( 2, len(result) )
    self.assertEqual( 'CCCC', result[1] )
    self.assertEqual( 4, len(result[0]) )
    self.assertTrue( result[0][0] in ('C', 'T', 'G') )
    self.assertTrue( result[0][1] in ('C', 'T', 'G') )
    self.assertTrue( result[0][2] in ('C', 'T', 'G') )
    self.assertTrue( result[0][3] in ('C', 'T', 'G') )
    self.assertEqual( 4, len(mutator.vcf.snp_map) )

class TestSamAccuracyEvaluator( unittest.TestCase ):
  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    acc = bio.SamAccuracyEvaluator( sam )
    self.assertEqual( 0, acc.stats['hard_clipping'] )
    self.assertEqual( 60.0, acc.stats['correct_mapq'] )
    self.assertEqual( 15, acc.stats['matched'] )
    self.assertEqual( 0, acc.stats['incorrect'] )
    self.assertEqual( 0, acc.stats['soft_clipping'] )
    self.assertEqual( 2, acc.stats['lines'] )
    self.assertEqual( 1, acc.stats['unknown_mapping'] )
    self.assertEqual( 1, acc.stats['correct'] )
    self.assertEqual( 0, acc.stats['mapped'] )
    self.assertEqual( 0, acc.stats['unmapped'] )
    self.assertEqual( 0.0, acc.stats['incorrect_mapq'] )

class TestBamReader( unittest.TestCase ):
  def setUp(self):
    pass

