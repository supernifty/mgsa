
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

  def test_simple_indel( self ):
    self.assertEqual( 0, bio.Distance( 'aa', 'aa' ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'aba', 'aa' ).simple_indel() )
    self.assertEqual( 1, bio.Distance( 'ba', 'bab' ).simple_indel() )
    self.assertEqual( -2, bio.Distance( 'baa', 'b' ).simple_indel() )
    self.assertEqual( 1, bio.Distance( 'b', 'ab' ).simple_indel() )
    self.assertEqual( 2, bio.Distance( 'b', 'aab' ).simple_indel() )
    self.assertEqual( -1, bio.Distance( 'cca', 'ca' ).simple_indel() )
