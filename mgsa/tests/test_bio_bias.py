
import unittest
import StringIO

import bio

class TestMappability( unittest.TestCase ):
  def test_simple(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nGTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    mapper = bio.Mappability( fq, self._matching_mapper_function, 5, min_pos=3, max_pos=5, log=bio.log_quiet )
    self.assertEqual( 0, mapper.summary['max_incorrect'] )
    self.assertEqual( 5, mapper.summary['read_length'] )
    self.assertEqual( 2, mapper.summary['genome_length'] )
    self.assertEqual( 0, mapper.summary['max_unmapped'] )
    self.assertEqual( 1.0, mapper.summary['max_accuracy'] ) 
    self.assertEqual( 1.0, mapper.summary['min_accuracy'] )
    self.assertEqual( [1.0, 1.0], mapper.accuracy )

  def test_incorrect(self):
    fasta = StringIO.StringIO( '>comment\nAAAAAAAAAA\nGTTTTTTTTT' )
    fq = bio.FastqPosGenerator( fasta, log=bio.log_quiet )
    mapper = bio.Mappability( fq, self._incorrect_mapper_function, 5, min_pos=3, max_pos=5, log=bio.log_quiet )
    self.assertEqual( 1., mapper.summary['max_incorrect'] )
    self.assertEqual( 5, mapper.summary['read_length'] )
    self.assertEqual( 2, mapper.summary['genome_length'] )
    self.assertEqual( 0, mapper.summary['max_unmapped'] )
    self.assertEqual( 0.0, mapper.summary['max_accuracy'] ) 
    self.assertEqual( 0.0, mapper.summary['min_accuracy'] )
    self.assertEqual( [0.0, 0.0], mapper.accuracy )

  def _matching_mapper_function( self, fq, sam, index=True ):
    # write a bogus matching fastq file to sam
    with open( sam, 'w' ) as sam_fh:
      sam_fh.write( '@HD     VN:1.0  SO:unsorted\n@SQ     SN:generated    LN:4023\n@PG     ID:dummy      PN:dummy      VN:2.1.0\nmgsa_seq_5~0~0       0       generated       6    44      5M     *       0       0       AAAAG      ~~~~~      AS:i:100        XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YT:Z:UU' )

  def _incorrect_mapper_function( self, fq, sam, index=True ):
    # write a bogus matching fastq file to sam
    with open( sam, 'w' ) as sam_fh:
      sam_fh.write( '@HD     VN:1.0  SO:unsorted\n@SQ     SN:generated    LN:4023\n@PG     ID:dummy      PN:dummy      VN:2.1.0\nmgsa_seq_8~0~0       0       generated       6    44      5M     *       0       0       AAAAG      ~~~~~      AS:i:100        XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YT:Z:UU' )

