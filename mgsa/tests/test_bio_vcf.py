
import unittest
import StringIO

import bio

class TestVCF( unittest.TestCase ):

  def setUp(self):
    pass

  def test_indel_position(self):
    vcf = bio.VCF()
    vcf.indel( 2, 'GG', 'GCCG' )
    vcf.indel( 4, 'GG', 'GAACCG' )
    # 01II23IIII45
    # 0123456789
    
    self.assertEqual( (0,0,0), vcf.candidate_position_to_reference_position( 0 ) ) 
    self.assertEqual( (1,0,0), vcf.candidate_position_to_reference_position( 1 ) ) 
    self.assertEqual( (1,1,2), vcf.candidate_position_to_reference_position( 2 ) ) 
    self.assertEqual( (1,2,2), vcf.candidate_position_to_reference_position( 3 ) ) 
    self.assertEqual( (2,0,0), vcf.candidate_position_to_reference_position( 4 ) ) 
    self.assertEqual( (3,1,4), vcf.candidate_position_to_reference_position( 6 ) ) 
    self.assertEqual( (3,2,4), vcf.candidate_position_to_reference_position( 7 ) ) 
    self.assertEqual( (3,3,4), vcf.candidate_position_to_reference_position( 8 ) ) 
    self.assertEqual( (3,4,4), vcf.candidate_position_to_reference_position( 9 ) ) 
    self.assertEqual( (4,0,0), vcf.candidate_position_to_reference_position( 10 ) ) 