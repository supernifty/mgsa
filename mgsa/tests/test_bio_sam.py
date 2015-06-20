
import unittest
import StringIO

import bio

class TestSamToVCF( unittest.TestCase ):

  def setUp(self):
    pass

  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AAAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )

  def test_snp(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 1, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )

  def test_insert(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      2M2I15M       *       0       0       AACCAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AAAAAAAAAA\nTTTTTTTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )

  def test_delete(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_0~0~0  0       generated       1      60      2M2D18M       *       0       0       AAAAAAAAAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'AACCAAAAAAAA\nTTTTTTTTTT' ) # 10A, 10T
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    #print target_vcf
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )
    self.assertEqual( 'ACC', target_vcf.manager.indel_list[0].before )
    self.assertEqual( 'A', target_vcf.manager.indel_list[0].after )
    
  def test_consensus_ins(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_0~0~0  0       generated       1      60      3M1I2M       *       0       0       TTAATT    ~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'TTATT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    genome = s.sam_to_fasta.fasta.genome
    self.assertEqual( [1., 1., 0, 1., 1.], genome['T'] )
    self.assertEqual( [0., 0., 1.], genome['A'] )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )
    self.assertEqual( 'AT', target_vcf.manager.indel_list[0].before )
    self.assertEqual( 'AAT', target_vcf.manager.indel_list[0].after )

  def test_consensus_del(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_0~0~0  0       generated       1      60      3M1D2M       *       0       0       TTATT    ~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'TTAATT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    genome = s.sam_to_fasta.fasta.genome
    self.assertEqual( [1., 1., 0, 0., 1., 1.], genome['T'] )
    self.assertEqual( [0., 0., 1.], genome['A'] )
    self.assertEqual( [0., 0., 0., 1.], genome['-'] )
    self.assertEqual( 0, len(target_vcf.snp_list) )
    self.assertEqual( 1, len(target_vcf.manager.indel_list) )
    self.assertEqual( 'AA', target_vcf.manager.indel_list[0].before )
    self.assertEqual( 'A', target_vcf.manager.indel_list[0].after )

  def test_soft_clip(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_0~0~0  0       generated       3      60      2S4M       *       0       0       GGTTTT    ~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    reference = StringIO.StringIO( 'CCTTTT' )
    target_vcf = bio.VCF()
    s = bio.SamToVCF.instance( sam, reference, target_vcf, bio.log_quiet )
    genome = s.sam_to_fasta.fasta.genome
    self.assertEqual( [0, 0, 1., 1., 1., 1.], genome['T'] )
    self.assertEqual( [0, 0], genome['G'] )
  
  def log(self, msg):
    print msg

class TestSamToMultiChromosomeVCF( unittest.TestCase ):
  def test_snp(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       chr1       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  0       chr2       6      60      15M       *       0       0       GGGCATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0',)
    reference = StringIO.StringIO( '>chr1\nNNNNNAAAAA\nTTTTTTTTTT\n>chr2\nNNNNNGGGGA\nTTTTTTTTTT\n' )
    target_vcf = bio.VCF()
    s = bio.SamToMultiChromosomeVCF( sam, reference, target_vcf, bio.log_quiet )
    self.assertEqual( 2, len(target_vcf.snp_list) )
    self.assertEqual( 0, len(target_vcf.manager.indel_list) )
    if target_vcf.snp_list[0]['chr'] == 'chr1':
      chr1pos = 0
    else:
      chr1pos = 1
    self.assertEqual( 'A', target_vcf.snp_list[chr1pos]['ref'] )
    self.assertEqual( 'G', target_vcf.snp_list[1-chr1pos]['ref'] )
    self.assertEqual( 7, target_vcf.snp_list[chr1pos]['pos'] )
    self.assertEqual( 8, target_vcf.snp_list[1-chr1pos]['pos'] )

class TestSamToFasta( unittest.TestCase ):

  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    s = bio.SamToFasta( sam, log=bio.log_quiet )
    # stats
    self.assertEqual( 1, s.stats['unknown_mapping'] )
    self.assertEqual( 2, s.stats['lines'] )
    # consensus (note sequence starts at 6)
    # move, result, best, coverage
    self.assertEqual( (1, 'A', 1.0, 1.0), s.fasta.consensus_at( 5 ) )
    self.assertEqual( (1, 'A', 1.0, 1.0), s.fasta.consensus_at( 6 ) )
    self.assertEqual( (1, 'C', 1.0, 1.0), s.fasta.consensus_at( 7 ) )

  def test_write(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    s = bio.SamToFasta( sam, log=bio.log_quiet )
    target = StringIO.StringIO()
    s.write( target )
    lines = target.getvalue().split( '\n' )
    self.assertEqual( 3, len(lines) )
    self.assertTrue( lines[0].startswith('>') )
    self.assertEqual( 'NNNNAACAATTTTTTTTTT', lines[1] )
 
class TestSamToMultiChromosomeFasta( unittest.TestCase ):

  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       chr1       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  0       chr2       6      60      15M       *       0       0       GGTAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0')
    s = bio.SamToMultiChromosomeFasta( sam, log=bio.log_quiet )
    # stats
    self.assertEqual( 2, s.stats['unknown_mapping'] )
    self.assertEqual( 3, s.stats['total_lines'] )
    self.assertEqual( 2, s.stats['lines'] )
    self.assertEqual( 1, s.fastas['chr1'].stats['unknown_mapping'] )
    self.assertEqual( 1, s.fastas['chr1'].stats['lines'] )
    self.assertEqual( 1, s.fastas['chr2'].stats['unknown_mapping'] )
    self.assertEqual( 1, s.fastas['chr2'].stats['lines'] )
    # consensus (note sequence starts at 6)
    # move, result, best, coverage
    self.assertEqual( (1, 'A', 1.0, 1.0), s.fastas['chr1'].fasta.consensus_at( 5 ) )
    self.assertEqual( (1, 'A', 1.0, 1.0), s.fastas['chr1'].fasta.consensus_at( 6 ) )
    self.assertEqual( (1, 'C', 1.0, 1.0), s.fastas['chr1'].fasta.consensus_at( 7 ) )
    self.assertEqual( (1, 'G', 1.0, 1.0), s.fastas['chr2'].fasta.consensus_at( 5 ) )
    self.assertEqual( (1, 'G', 1.0, 1.0), s.fastas['chr2'].fasta.consensus_at( 6 ) )
    self.assertEqual( (1, 'T', 1.0, 1.0), s.fastas['chr2'].fasta.consensus_at( 7 ) )

class TestFastaMutate( unittest.TestCase ):

  def test_simple(self):
    l = ( 'AAAA', 'CCCC')
    fr = bio.FastaReader( StringIO.StringIO( '\n'.join( l ) ) )
    mutator = bio.FastaMutate( fr, snp_prob=1., log=bio.log_quiet )
    result = []
    for fragment in mutator.items():
      result.append( fragment )
    self.assertEqual( 2, len(result) )
    self.assertEqual( 'CCCC', result[1] )
    self.assertEqual( 4, len(result[0]) )
    self.assertTrue( result[0][0] in ('C', 'T', 'G') )
    self.assertTrue( result[0][1] in ('C', 'T', 'G') )
    self.assertTrue( result[0][2] in ('C', 'T', 'G') )
    self.assertTrue( result[0][3] in ('C', 'T', 'G') )
    self.assertEqual( 4, len(mutator.vcf.snp_map) )

class TestSamDiff( unittest.TestCase ):
  def test_simple(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet )
    # diff.totals
    self.assertEqual( 1, len(diff.totals.keys()))
    self.assertEqual( 1, diff.totals[3] )

  def test_simple_compare_pos(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet, compare_position=True )
    # diff.totals
    self.assertEqual( 1, len(diff.totals.keys()))
    self.assertEqual( 1, diff.totals[3] ) # 11 -> 1
    self.assertEqual( 2, len(diff.position_totals.keys()))
    self.assertEqual( 1, diff.position_totals[1] ) # 1 -> 1
    self.assertEqual( 1, diff.position_totals[2] ) # 2 -> 1

  def test_simple_subset_quality(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      40      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet, compare_position=True, subset_detail=True )
    self.assertEqual( [40], diff.mapq_totals[1] ) # 01
    self.assertEqual( [60], diff.mapq_totals[2] ) # 10
    self.assertEqual( [], diff.mapq_totals[3] ) # 10
    self.assertEqual( 40.0, diff.mapq_subset_stats[1]['mean'] )
    self.assertEqual( 60.0, diff.mapq_subset_stats[2]['mean'] )

  def test_simple_subset_quality_agree(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      40      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet, compare_position=True, subset_detail=True )
    self.assertEqual( [40, 60], diff.mapq_totals[3] ) # 11
    self.assertEqual( 50.0, diff.mapq_subset_stats[3]['mean'] )
    self.assertEqual( 1, diff.mapq_subset_stats[3]['hist'][8] )
    self.assertEqual( 1, diff.mapq_subset_stats[3]['hist'][12] )

  def test_split(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  4       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet )
    # diff.totals
    self.assertEqual( 1, len(diff.totals.keys()))
    self.assertEqual( 1, diff.totals[1] )

  def test_unmapped(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  4       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  4       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet )
    # diff.totals
    self.assertEqual( 1, len(diff.totals.keys()))
    self.assertEqual( 1, diff.totals[0] )

  def test_mismatch_tracking_no_mismatch(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      40      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet, compare_position=True, subset_detail=True, mismatch_detail=1 )
    self.assertEqual( 0, len(diff.mismatch_stats) )

  def test_mismatch_tracking_with_mismatch(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      40      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet, compare_position=True, subset_detail=True, mismatch_detail=1 )
    self.assertEqual( 1, len(diff.mismatch_stats) )
    self.assertEqual( 8, diff.mismatch_stats['mgsa_seq_5~0~0']['pos'] )
    self.assertEqual( 6, diff.mismatch_stats['mgsa_seq_5~0~0']['alt'] )

  def test_mismatch_tracking_with_mismatch(self):
    sam1 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      40      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    sam2 = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       8      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    diff = bio.SamDiff( [ sam1, sam2], log=bio.log_quiet, compare_position=True, subset_detail=True, mismatch_detail=0 )
    self.assertEqual( 1, len(diff.mismatch_stats) )
    self.assertEqual( 8, diff.mismatch_stats['mgsa_seq_5~0~0']['alt'] )
    self.assertEqual( 6, diff.mismatch_stats['mgsa_seq_5~0~0']['pos'] )

class TestSamAccuracyEvaluator( unittest.TestCase ):
  def test_simple(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTTTTTTT    ~~~~~~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    acc = bio.SamAccuracyEvaluator( sam )
    self.assertEqual( 0, acc.stats['hard_clipping'] )
    self.assertEqual( 60.0, acc.stats['correct_mapq'] )
    self.assertEqual( 15, acc.stats['matched'] )
    self.assertEqual( 0, acc.stats['incorrect'] )
    self.assertEqual( 0, acc.stats['soft_clipping'] )
    self.assertEqual( 2, acc.stats['lines'] )
    self.assertEqual( 1, acc.stats['unknown_mapping'] )
    self.assertEqual( 1, acc.stats['correct'] )
    self.assertEqual( 0, acc.stats['mapped'] )
    self.assertEqual( 0, acc.stats['unmapped'] )
    self.assertEqual( 0.0, acc.stats['incorrect_mapq'] )

class TestSamStats( unittest.TestCase ):
  def test_gc(self):
    sam = ( '@SQ     SN:generated    LN:4023', 'mgsa_seq_5~0~0  0       generated       6      60      15M       *       0       0       AACAATTTTT    ~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', 'mgsa_seq_5~0~0  4       generated       6      60      15M       *       0       0       AAGCATTTTT    ~~~~~~~~~~~~~~~    NM:i:10 AS:i:84 XS:i:0', )
    stats = bio.SamStats( sam )
    self.assertEqual( [0.1], stats.mapped['gc'] )
    self.assertEqual( [0.2], stats.unmapped['gc'] )
 
class TestBamReaderExternal( unittest.TestCase ):
  def test_simple(self):
    ext = bio.BamReaderExternal( "echo %s", "'hello'" )
    lines = [ x for x in ext ]
    self.assertEqual( 1, len(lines) )
    self.assertEqual( 'hello\n', lines[0] )

  def test_multiple(self):
    ext = bio.BamReaderExternal( "echo 'hello\n%s'", "goodbye" )
    lines = [ x for x in ext ]
    self.assertEqual( 2, len(lines) )
    self.assertEqual( 'hello\n', lines[0] )
    self.assertEqual( 'goodbye\n', lines[1] )

class TestBamReader( unittest.TestCase ):
  def setUp(self):
    pass

