
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

  def test_variations(self):
    vcf = bio.VCF()
    vcf.snp( 5, 'G', 'T' )
    vcf.snp( 12, 'C', 'T' )
    vcf.indel( 7, 'AA', 'ATTA' )
    vcf.indel( 10, 'AA', 'AGGA' )
    self.assertEqual( 2, len(vcf.variations( 0, 10 )) )
    self.assertTrue( 'I7-2[7,0,10]' in vcf.variations( 0, 10 ) )
    self.assertTrue( 'S5-0' in vcf.variations( 0, 10 ) )
    self.assertEqual( 1, len(vcf.variations( 10, 20 )) )
    self.assertTrue( 'S2-0' in vcf.variations( 10, 20 ) )
