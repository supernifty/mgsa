
import StringIO
import unittest

import bio

class TestChart( unittest.TestCase ):

  def test_multi_chromosome_comparison( self ):
    target = StringIO.StringIO()
    writer = bio.VCFWriter(target)
    writer.snp( pos=1, ref='A', alt='G', chromosome='X' )
    writer.snp( pos=2, ref='G', alt='T', chromosome='X' )
    writer.snp( pos=1, ref='C', alt='T', chromosome='Y' )
    target.seek(0)
    parent = bio.MultiChromosomeVCF( target )

    target = StringIO.StringIO()
    writer = bio.VCFWriter(target)
    writer.snp( pos=2, ref='G', alt='T', chromosome='X' )
    writer.snp( pos=1, ref='C', alt='T', chromosome='Y' )
    writer.snp( pos=2, ref='C', alt='T', chromosome='Y' )
    target.seek(0)
    child = bio.MultiChromosomeVCF( target )

    target = StringIO.StringIO()
    writer = bio.VCFWriter(target)
    writer.snp( pos=1, ref='A', alt='G', chromosome='X' )
    writer.snp( pos=1, ref='C', alt='T', chromosome='Y' )
    writer.snp( pos=2, ref='C', alt='T', chromosome='Y' )
    target.seek(0)
    truth = bio.MultiChromosomeVCF( target )

    result = bio.multi_chromosome_comparison( parent, child, truth )
    self.assertEqual( [(0, 1, 'X'), (3, 0, 'Y')], result['chromosome_names'] )
    self.assertEqual( [3, 6], result['chromosome_offsets'] )
    self.assertEqual( [2, 4], result['x_tp'] ) # parent, child match
    self.assertEqual( [5], result['x_fp'] ) # child only (Y-2)
    self.assertEqual( [1], result['x_fn'] ) # parent only (X-1)
    self.assertEqual( [5], result['x_t'] ) # truth and child match (Y-2) but parent doesn't
