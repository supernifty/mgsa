
import unittest

import bio

class TestFeatures(unittest.TestCase):
  
  def test_entropy(self):
    self.assertEqual( 0.0, bio.ReadFeature( 'AAAAAAAAAA' ).entropy() )
    self.assertAlmostEqual( 0.40922323189505994, bio.ReadFeature( 'ACAAAAAAAA' ).entropy() )
    self.assertAlmostEqual( 1.0, bio.ReadFeature( 'ACGTTAAGAT' ).entropy() )
