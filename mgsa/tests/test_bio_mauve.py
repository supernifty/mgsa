
import StringIO
import unittest

import bio

class TestMauveMap(unittest.TestCase):
  
  def test_simple(self):
    src = StringIO.StringIO( '> 1:1-4 +\nA--BC-D\n> 2:1-7 +\nAxyBCzD\n=' )
    m = bio.MauveMap( src, log=bio.log_quiet )
    self.assertEqual( 1, m.coverage[1] )
    self.assertEqual( 4, m.coverage[2] )
    self.assertEqual( 5, m.coverage[3] )
    self.assertEqual( 7, m.coverage[4] )
    self.assertEqual( 4, len(m.coverage))
  
  def test_simple_2(self):
    src = StringIO.StringIO( '> 1:1-4 +\nA--BC-D\n> 2:1-7 +\nAxyBCzD\n=' )
    m = bio.MauveMap( src, src_strand=2, target_strand=1, log=bio.log_quiet )
    self.assertEqual( 1, m.coverage[1] )
    self.assertEqual( 2, m.coverage[4] )
    self.assertEqual( 3, m.coverage[5] )
    self.assertEqual( 4, m.coverage[7] )
    self.assertEqual( 4, len(m.coverage))

  def test_remap_ok( self ):
    src = StringIO.StringIO( '> 1:1-10 +\n1234567890\n> 2:3-13 +\n3456789012\n=' )
    target = StringIO.StringIO()
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       3      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    m = bio.MauveMap( src, src_strand=1, target_strand=2, log=bio.log_quiet )
    m.remap( sam1, target )
    self.assertEqual( 1,  m.stats['mapped'] )
    self.assertEqual( '5', target.getvalue().split( '\n' )[1].split('\t')[3] )
    
  def test_remap_partial( self ):
    src = StringIO.StringIO( '> 1:1-10 +\n1234567890\n> 2:3-13 +\n345--89012\n=' )
    target = StringIO.StringIO()
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       3      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    m = bio.MauveMap( src, src_strand=1, target_strand=2, log=bio.log_quiet )
    m.remap( sam1, target )
    self.assertEqual( 1,  m.stats['reads_partial'] )
    self.assertEqual( '0', target.getvalue().split( '\n' )[1].split('\t')[3] )
    
  def test_remap_unmapped( self ):
    src = StringIO.StringIO( '> 1:1-10 +\n1234567890\n> 2:3-13 +\n34------12\n=' )
    target = StringIO.StringIO()
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       3      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    m = bio.MauveMap( src, src_strand=1, target_strand=2, log=bio.log_quiet )
    m.remap( sam1, target )
    self.assertEqual( 1,  m.stats['reads_notcovered'] )
    self.assertEqual( '0', target.getvalue().split( '\n' )[1].split('\t')[3] )

  def test_remap_multiple( self ):
    src = StringIO.StringIO( '> 1:1-10 +\n1234567890\n> 2:3-13 +\n3456789012\n=\n> 1:50-60 +\n1234567890\n> 2:80-90 +\n1234567890\n=' )
    target = StringIO.StringIO()
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       3      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  0       generated       53      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    m = bio.MauveMap( src, src_strand=1, target_strand=2, log=bio.log_quiet )
    m.remap( sam1, target )
    self.assertEqual( 2,  m.stats['mapped'] )
    self.assertEqual( 20,  m.genome_stats['count'] )
    self.assertEqual( 1,  m.genome_stats['xmin'] )
    self.assertEqual( 59,  m.genome_stats['xmax'] )
    self.assertEqual( 3,  m.genome_stats['ymin'] )
    self.assertEqual( 89,  m.genome_stats['ymax'] )
    self.assertEqual( '5', target.getvalue().split( '\n' )[1].split('\t')[3] )
    self.assertEqual( '83', target.getvalue().split( '\n' )[2].split('\t')[3] )
    
  def test_remap_reverse( self ):
    src = StringIO.StringIO( '> 1:1-10 +\n1234567890\n> 2:3-13 +\n3456789012\n=\n> 1:50-60 -\n0987654321\n> 2:80-90 +\n1234567890\n=' )
    target = StringIO.StringIO()
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       3      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  0       generated       53      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    m = bio.MauveMap( src, src_strand=1, target_strand=2, log=bio.log_quiet )
    m.remap( sam1, target )
    self.assertEqual( '5', target.getvalue().split( '\n' )[1].split('\t')[3] )
    self.assertEqual( '87', target.getvalue().split( '\n' )[2].split('\t')[3] )
    
  def test_remap_reference( self ):
    src = StringIO.StringIO( '> 1:1-10 +\n1234567890\n> 2:3-13 +\n3456789012\n=\n> 1:50-60 -\n0987654321\n> 2:80-90 +\n1234567890\n=' )
    target = StringIO.StringIO()
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       3      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  0       generated       53      60      5M       *       0       0       AACAA    ~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    m = bio.MauveMap( src, src_strand=1, target_strand=2, new_reference='hello', log=bio.log_quiet )
    m.remap( sam1, target )
    self.assertEqual( '@SQ     SN:hello    LN:4023', target.getvalue().split( '\n' )[0] )
    self.assertEqual( 'hello', target.getvalue().split( '\n' )[1].split('\t')[2] )
    self.assertEqual( 'hello', target.getvalue().split( '\n' )[2].split('\t')[2] )
            
