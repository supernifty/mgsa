
import unittest
import StringIO

import bio

class TestConfig(unittest.TestCase):
  
  def test_read_config_file( self ):
    s = StringIO.StringIO( '#comment\ncommand value' )
    c = bio.Config()
    cfg = c.read_config_file( s )
    self.assertEqual( 'value', cfg['command'] )
