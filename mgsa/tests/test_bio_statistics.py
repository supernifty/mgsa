
import unittest

import bio

class TestStatistics(unittest.TestCase):
  
  def test_bucket(self):
    self.assertEqual( [1, 0], bio.bucket( [10], [0,11] ) )
    self.assertEqual( [0, 1], bio.bucket( [10], [0,10] ) )
    self.assertEqual( [1, 1], bio.bucket( [0,10], [0,10] ) )

  def test_sign_bucket(self):
    b = bio.sign_bucket( [ -1, 0, 1, 2, 3 ] )
    self.assertEqual( 0.6, b['positive'] )
    self.assertEqual( 0.2, b['negative'] )
    self.assertEqual( 0.2, b['zero'] )
    self.assertAlmostEqual( 0.75, b['positive_proportion'] )
