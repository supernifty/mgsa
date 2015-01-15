
import unittest
import StringIO

import bio

class TestBiasReport( unittest.TestCase ):

  def setUp(self):
    pass

  def test_simple(self):
    sam_file = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       1      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  0       generated       1      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference_fasta = StringIO.StringIO( 'AAAAA\nTTTTTTTTTT' )
    vcf = bio.VCF( log=self.quiet_log )
    vcf.snp( 2, 'A', 'C' ) # correct
    vcf.snp( 10, 'T', 'G' ) # incorrect

    report = bio.BiasReport(
      candidate_fasta=bio.SamToFasta( sam=sam_file, log=self.quiet_log ).fasta, # what the aligner has built
      reference_fasta=reference_fasta, # the reference
      donor_vcf=vcf, # variations to donor
      log=self.quiet_log,
      buckets_count=3 )
    
    self.assertEqual( 13, report.stats['mapped_no_variation'] ) # mapped ok without any variation
    self.assertEqual( 2, report.stats['total'] ) # total variations
    self.assertEqual( 0.5, report.stats['mean_reference'] ) # 1 good, 1 bad
    self.assertEqual( 2, report.stats['mapped'] ) # variation
    self.assertEqual( 0, report.stats['unmapped'] ) # variation
    self.assertEqual( 1.0, report.stats['reference_percent_total'] ) # 1 totally bad
    self.assertEqual( 1.0, report.stats['incorrect_percent_total'] ) # 1 also totally incorrect
    self.assertEqual( 0.0, report.stats['error_percent_total'] ) # but not due to an error
    self.assertEqual( 0.0, report.stats['mean_error'] ) # no error
    self.assertEqual( 0, report.stats['unmapped_no_variation'] ) # no unmapped

    self.assertEqual( [0.5, 0.0, 0.0, 0.5], report.reference_histogram ) # 50% has no reference bias, 50% has high reference bias
    self.assertEqual( [1.0, 0.0, 0.0, 0.0], report.error_histogram ) # 100% has no error

  def log(self, msg):
    print msg

  def quiet_log(self, msg):
    pass


if __name__ == '__main__':
  unittest.main()
