
import unittest
import StringIO

import bio

class TestDistance( unittest.TestCase ):

  def test_hamming( self ):
    self.assertEqual( 0, bio.Distance( 'aa', 'aa' ).hamming() )
    self.assertEqual( 1, bio.Distance( 'ba', 'aa' ).hamming() )
    self.assertEqual( 1, bio.Distance( 'ba', 'bb' ).hamming() )
    self.assertEqual( 2, bio.Distance( 'ba', 'ab' ).hamming() )
    self.assertRaises(Exception, lambda: bio.Distance( 'b', 'ab' ).hamming() )
