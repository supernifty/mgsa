
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

  def test_roc(self):
    hist_n = [ 3, 2, 1 ]
    hist_p = [ 1, 2, 3 ]
    # positives are at lower values
    rocs = bio.generate_roc( hist_n, hist_p )

    self.assertAlmostEqual( 0, rocs['fpr'][0] )
    self.assertAlmostEqual( 0.5, rocs['fpr'][1] )
    self.assertAlmostEqual( 5./6, rocs['fpr'][2] )

    self.assertAlmostEqual( 0, rocs['tpr'][0] )
    self.assertAlmostEqual( 1./6, rocs['tpr'][1] )
    self.assertAlmostEqual( 0.5, rocs['tpr'][2] )
