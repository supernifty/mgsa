
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
    parent = bio.MultiChromosomeVCF( target, log=bio.log_quiet ) # parent: X 1, 2; Y 1

    target = StringIO.StringIO()
    writer = bio.VCFWriter(target)
    writer.snp( pos=2, ref='G', alt='T', chromosome='X' )
    writer.snp( pos=1, ref='C', alt='T', chromosome='Y' )
    writer.snp( pos=2, ref='C', alt='T', chromosome='Y' )
    target.seek(0)
    child = bio.MultiChromosomeVCF( target, log=bio.log_quiet ) # child: X 2; Y 1 2

    target = StringIO.StringIO()
    writer = bio.VCFWriter(target)
    writer.snp( pos=1, ref='A', alt='G', chromosome='X' )
    writer.snp( pos=1, ref='C', alt='T', chromosome='Y' )
    writer.snp( pos=2, ref='C', alt='T', chromosome='Y' )
    target.seek(0)
    truth = bio.MultiChromosomeVCF( target, log=bio.log_quiet ) # truth: X 1; Y 1 2

    result = bio.multi_chromosome_comparison( parent, child, truth )
    self.assertEqual( [(0, 1, 'Y'), (3, 0, 'X')], result['chromosome_names'] )
    self.assertEqual( [3, 6], result['chromosome_offsets'] )
    self.assertEqual( [1, 5], result['x_tp'] ) # parent, child match
    self.assertEqual( [2], result['x_fp'] ) # child only (Y-2)
    self.assertEqual( [4], result['x_fn'] ) # parent only (X-1)
    self.assertEqual( [2], result['x_t'] ) # truth and child match (Y-2) but parent doesn't
