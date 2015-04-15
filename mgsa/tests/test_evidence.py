
import unittest
import StringIO

import bio

class MockFasta(object):
  def __init__( self, depths ):
    '''
      @depths: array of depths
    '''
    self.depths = depths
    self.length = len(depths)

  def count( self, pos ):
    return { 'A': self.depths[pos] }

  def consensus( self, start, end ):
    return 'N' * ( end - start )

class MockReference(object):
  def __init__( self ):
    pass

  def fragment( self, start, end ):
    return 'N' * ( end - start )

  def read_to( self, pos ):
    return True
    
class TestEvidence( unittest.TestCase ):

  def setUp(self):
    pass

  def test_linear_evidence_good(self):
    e = bio.LinearShapeDifference( 0, 1, 3, 4 ) # 1234
    good = bio.ProbabilisticFasta()
    good.add( 'AACGTT', 1 )
    good.add( 'AACGTT', 2 )
    good.add( 'AACGTT', 3 )
    good.add( 'AACGTT', 4 )
    self.assertEqual( 1., e.evidence( good, 1 ) )

  def test_linear_evidence_good_2(self):
    e = bio.LinearShapeDifference( 0, 4, 3, 1 ) # 4321
    good = bio.ProbabilisticFasta()
    good.add( 'AACGTT', 1 )
    good.add( 'AACGT', 1 )
    good.add( 'AACG', 1 )
    good.add( 'AAC', 1 )
    self.assertEqual( 1., e.evidence( good, 3 ) )

  def test_linear_evidence_ok(self):
    e = bio.LinearShapeDifference( 0, 4, 3, 1, log=bio.log_quiet ) # 4321
    good = bio.ProbabilisticFasta() # 4[4432]1
    good.add( 'AACGTT', 1 )
    good.add( 'AACGT', 1 )
    good.add( 'AACG', 1 )
    good.add( 'AAC', 1 )
    #self.assertAlmostEqual( 0.5921640105587143, e.evidence( good, 2 ) )
    #self.assertAlmostEqual( 0.3291796067500631, e.evidence( good, 2 ) )
    self.assertAlmostEqual( 0.9256423851565307, e.evidence( good, 2, mean_depth=1. ) )

  def test_linear_evidence_bad(self):
    e = bio.LinearShapeDifference( 0, 4, 3, 1 ) # 4321
    good = bio.ProbabilisticFasta()
    good.add( 'AACGTT', 1 )
    good.add( 'AACGTT', 2 )
    good.add( 'AACGTT', 3 )
    good.add( 'AACGTT', 4 )
    #self.assertAlmostEqual( 0.19311601249364552, e.evidence( good, 1 ) )
    self.assertAlmostEqual( 0.1417505640449146, e.evidence( good, 1, mean_depth=1. ) )

  def test_multi_shape_good(self):
    s1 = bio.LinearShapeDifference( -3, 1, -1, 3 ) # 123
    s2 = bio.LinearShapeDifference( 0, 4, 3, 1 ) # 4321
    e = bio.MultipleLinearShapeDifference( [ s1, s2 ] )
    good = bio.ProbabilisticFasta()
    good.add( 'AACG', 1 )
    good.add(  'AACG', 2 )
    good.add(   'AACG', 3 )
    good.add(    'AACG', 4 )
    # 1234321
    self.assertAlmostEqual( 1, e.evidence( good, 4 ) )
    
  def test_multi_shape_bad(self):
    s1 = bio.LinearShapeDifference( -3, 4, 0, 1 ) # 4321
    s2 = bio.LinearShapeDifference( 1, 2, 3, 4 ) # 234
    e = bio.MultipleLinearShapeDifference( [ s1, s2 ] )
    bad = bio.ProbabilisticFasta()
    bad.add( 'AACG', 1 )
    bad.add(  'AACG', 2 )
    bad.add(   'AACG', 3 )
    bad.add(    'AACG', 4 )
    # 1234321
    self.assertAlmostEqual( 0, e.evidence( bad, 4 ) )
    
  def test_deletion_detector_good(self):
    e = bio.LongDeletionDetector( 100 ) #, log=bio.log_stderr )
    depths = [ 10, 50, 30, 80, 0 ]
    depths.extend( range( 100, 69, -1 ) )
    depths.extend( [ 0 ] * 100 )
    depths.extend( range( 50, 100 ) )
    good = MockFasta( depths )
    self.assertAlmostEqual( 1, e.evidence( good, 35, mean_depth=100. ) )

  def test_deletion_detector_good_vcf(self):
    v = bio.VCF()
    e = bio.LongDeletionDetector( 100, vcf_writer=v, threshold=0.5, reference_fasta=MockReference())#, log=bio.log_stderr )
    depths = [ 10, 50, 30, 80, 0 ] # leading depths
    depths.extend( range( 100, 69, -1 ) ) # breakpoint
    depths.extend( [ 0 ] * 100 ) # deletion
    depths.extend( range( 50, 100 ) ) # breakpoint
    good = MockFasta( depths )
    e.evidence( good, 30, mean_depth=100. )
    e.evidence( good, 31, mean_depth=100. )
    e.evidence( good, 32, mean_depth=100. )
    e.evidence( good, 33, mean_depth=100. )
    e.evidence( good, 34, mean_depth=100. )
    self.assertAlmostEqual( 1, e.evidence( good, 35, mean_depth=100. ) )
    e.evidence( good, 36, mean_depth=100. )
    e.evidence( good, 37, mean_depth=100. )
    e.evidence( good, 38, mean_depth=100. )
    e.evidence( good, 39, mean_depth=100. )
    self.assertEqual( 1, len( v.manager.indel_list ) )
    self.assertEqual( 35, v.manager.indel_list[0].pos )
    self.assertEqual( 100, len( v.manager.indel_list[0].before ) - len( v.manager.indel_list[0].after )  )
