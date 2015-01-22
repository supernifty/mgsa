
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

  def test_short_indel(self):
    vcf = bio.VCF()
    vcf.indel( 5, 'AA', 'ATA' )
    self.assertEqual( 1, len(vcf.variations( 0, 10 )) )
    self.assertTrue( 'I5-1[5,0,10]' in vcf.variations( 0, 10 ) )

class TestVCFWriter( unittest.TestCase ):
  def setUp(self):
    pass

  def test_write_snp_simple( self ):
    target = StringIO.StringIO()
    vcf = bio.VCFWriter(target)
    vcf.snp( pos=1, ref='A', alt='G' )
    lines = target.getvalue().split('\n')
    self.assertEqual( 6, len(lines) )
    self.assertEqual( '.\t2\t.\tA\tG\t.\tPASS\tDP=10\tGT\t1', lines[4] )

  def test_write_snp_confidence( self ):
    target = StringIO.StringIO()
    vcf = bio.VCFWriter(target)
    vcf.snp( pos=1, ref='A', alt='G', confidence=0.5 )
    lines = target.getvalue().split('\n')
    self.assertEqual( 6, len(lines) )
    self.assertEqual( '.\t2\t.\tA\tG\t3\tPASS\tDP=10\tGT\t1', lines[4] )

class TestVCFDiff( unittest.TestCase ):
  def setUp(self):
    pass

  def test_simple( self ):
    correct = bio.VCF()
    candidate = bio.VCF()
    correct.snp( 2, 'G', 'T' )
    correct.snp( 5, 'G', 'T' )
    candidate.snp( 5, 'G', 'T' )
    candidate.snp( 8, 'G', 'T' )
    diff = bio.VCFDiff( correct, candidate, log = bio.log_quiet, generate_positions = True )
    self.assertEqual( 1, diff.stats['tp'] )
    self.assertEqual( 1, diff.stats['fp'] )
    self.assertEqual( 1, diff.stats['fn'] )
    self.assertEqual( 5, diff.positions['tp'][0] )
    self.assertEqual( 2, diff.positions['fn'][0] )
    self.assertEqual( 8, diff.positions['fp'][0] )
