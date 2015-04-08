
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

  def test_load(self):
    src = StringIO.StringIO( '.\t3\t.\tA\tG\t3.0\tPASS\tDP=10\tGT\t1' )
    vcf = bio.VCF( reader=src, log=bio.log_quiet )
    self.assertEqual( 1, len(vcf.snp_list) )
    snp = vcf.snp_list[0]
    self.assertEqual( 'A', snp['ref'] )
    self.assertEqual( 'G', snp['alt'] )
    self.assertEqual( 2, snp['pos'] )

  def test_breakpoints(self):
    vcf = bio.VCF()
    vcf.indel( 5, 'AA', 'ATA' ) # insertion
    self.assertEqual( [5, 6], vcf.breakpoints( all_affected=False ) )
    self.assertEqual( [5, 6],  vcf.breakpoints( all_affected=True ) )

    vcf = bio.VCF()
    vcf.indel( 5, 'ATA', 'AA' ) # deletion
    self.assertEqual( [5, 7], vcf.breakpoints( all_affected=False ) )
    self.assertEqual( [5, 6, 7],  vcf.breakpoints( all_affected=True ) )

class TestMultiChromosomeVCF( unittest.TestCase ):
  def test_read( self ):
    target = StringIO.StringIO()
    vcf = bio.VCFWriter(target)
    vcf.snp( pos=1, ref='A', alt='G', chromosome='X' )
    vcf.snp( pos=1, ref='C', alt='T', chromosome='Y' )
    lines = target.getvalue().split('\n')
    self.assertEqual( 'X\t2\t.\tA\tG\t3.0\tPASS\tDP=10\tGT\t1', lines[4] )
    self.assertEqual( 'Y\t2\t.\tC\tT\t3.0\tPASS\tDP=10\tGT\t1', lines[5] )
    target.seek(0)
    multi_vcf = bio.MultiChromosomeVCF( target, log=bio.log_quiet )
    self.assertEqual( 2, len( multi_vcf.vcfs ) )
    self.assertEqual( 'A', multi_vcf.vcfs['X'].snp_list[0]['ref'] )
    self.assertEqual( 'G', multi_vcf.vcfs['X'].snp_list[0]['alt'] )
    self.assertEqual( 'C', multi_vcf.vcfs['Y'].snp_list[0]['ref'] )
    self.assertEqual( 'T', multi_vcf.vcfs['Y'].snp_list[0]['alt'] )

class TestVCFWriter( unittest.TestCase ):
  def setUp(self):
    pass

  def test_write_snp_simple( self ):
    target = StringIO.StringIO()
    vcf = bio.VCFWriter(target)
    vcf.snp( pos=1, ref='A', alt='G' )
    lines = target.getvalue().split('\n')
    self.assertEqual( 6, len(lines) )
    self.assertEqual( '.\t2\t.\tA\tG\t3.0\tPASS\tDP=10\tGT\t1', lines[4] ) # 3 is default quality

  def test_write_snp_confidence( self ):
    target = StringIO.StringIO()
    vcf = bio.VCFWriter(target)
    vcf.snp( pos=1, ref='A', alt='G', confidence=0.5 )
    lines = target.getvalue().split('\n')
    self.assertEqual( 6, len(lines) )
    self.assertEqual( '.\t2\t.\tA\tG\t3.0\tPASS\tDP=10\tGT\t1', lines[4] )

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

  def test_confidence_map( self ):
    correct = bio.VCF()
    candidate = bio.VCF()
    correct.snp( 2, 'G', 'T' )
    correct.snp( 5, 'G', 'T' )
    candidate.snp( 5, 'G', 'T' )
    candidate.snp( 8, 'G', 'T' )
    diff = bio.VCFDiff( correct, candidate, log = bio.log_quiet, generate_positions = True )
    x, y = diff.confidence_map( 'tp' )
    self.assertEqual( 1, len(x) )
    self.assertEqual( 5, x[0] )
    self.assertAlmostEqual( 3.0102999566398116, y[0] ) # 0.5 conf
    x, y = diff.confidence_map( 'fp' )
    self.assertEqual( 1, len(x) )
    self.assertEqual( 8, x[0] )
    self.assertAlmostEqual( 3.0102999566398116, y[0] ) # 0.5 conf
    x, y = diff.confidence_map( 'fn' )
    self.assertEqual( 1, len(x) )
    self.assertEqual( 2, x[0] )
    self.assertAlmostEqual( 3.0102999566398116, y[0] ) # 0.5 conf
 
class TestVCFFastaDiff( unittest.TestCase ):

  def test_snp_fn( self ):
    fa = bio.ProbabilisticFasta(log=bio.log_quiet)
    fa.add( fragment='gg', start=1 )
    fa.add( fragment='cc', start=3 )
    vcf = bio.VCF()
    vcf.snp( pos=2, ref='g', alt='t' )
    diff = bio.VCFFastaDiff( vcf, fa, log=bio.log_quiet )
    self.assertEqual( 1, diff.snp_stats['fn'] )
    self.assertEqual( 0, diff.snp_stats['tp'] )
    self.assertEqual( 0, diff.snp_stats['other'] )

  def test_snp_tp( self ):
    fa = bio.ProbabilisticFasta(log=bio.log_quiet)
    fa.add( fragment='gt', start=1 )
    fa.add( fragment='cc', start=3 )
    vcf = bio.VCF()
    vcf.snp( pos=2, ref='g', alt='t' )
    diff = bio.VCFFastaDiff( vcf, fa, log=bio.log_quiet )
    self.assertEqual( 0, diff.snp_stats['fn'] )
    self.assertEqual( 1, diff.snp_stats['tp'] )
    self.assertEqual( 0, diff.snp_stats['other'] )

  def test_snp_other( self ):
    fa = bio.ProbabilisticFasta(log=bio.log_quiet)
    fa.add( fragment='ga', start=1 )
    fa.add( fragment='cc', start=3 )
    vcf = bio.VCF()
    vcf.snp( pos=2, ref='g', alt='t' )
    diff = bio.VCFFastaDiff( vcf, fa, log=bio.log_quiet )
    self.assertEqual( 0, diff.snp_stats['fn'] )
    self.assertEqual( 0, diff.snp_stats['tp'] )
    self.assertEqual( 1, diff.snp_stats['other'] )

