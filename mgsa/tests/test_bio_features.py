
import unittest

import bio

class TestFeatures(unittest.TestCase):
  
  def test_entropy(self):
    self.assertEqual( 0.0, bio.ReadFeature( 'AAAAAAAAAA' ).entropy() )
    self.assertAlmostEqual( 0.40922323189505994, bio.ReadFeature( 'ACAAAAAAAA' ).entropy() )
    self.assertAlmostEqual( 1.0, bio.ReadFeature( 'ACGTTAAGAT' ).entropy() )
    self.assertEqual( 0.0, bio.ReadFeature( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT' ).entropy() )

  def test_depth(self):
    p = bio.ProbabilisticFasta()
    p.add( 'AAAAA', 0, 1 )
    p.add( 'CCCCC', 1, 1 )
    p.add( 'GGGGG', 2, 1 ) # i.e. A AC ACG ACG ACG CG G
    
    mf = bio.MapFeature( p )
    self.assertEqual( 1, mf.depth( 0 ) )
    self.assertEqual( 2, mf.depth( 1 ) )
    self.assertEqual( 3, mf.depth( 2 ) )
    self.assertEqual( 3, mf.depth( 3 ) )
    self.assertEqual( 3, mf.depth( 4 ) )
    self.assertEqual( 2, mf.depth( 5 ) )
    self.assertEqual( 1, mf.depth( 6 ) )
 
