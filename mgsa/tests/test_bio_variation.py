
import unittest
import StringIO

import bio

class TestIndelVariation( unittest.TestCase ):

  def setUp(self):
    pass

  def test_all_same_base(self):
    x = bio.IndelVariation( 0, '', '' )
    self.assertTrue( x.all_same_base( 'GGG' ) )
    self.assertFalse( x.all_same_base( 'GCG' ) )

  def test_exact_insertion_ok(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
  
  def test_inexact_insertion_ok(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 20, 'GG', 'GCCTGAAG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )

  def test_insertion_overlap(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 21, 'CG', 'CCTTAAGG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
  
  def test_insertion_no_overlap(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 22, 'CG', 'CCTTAAGG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )
  
  def test_insertion_no_overlap_content(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 21, 'CG', 'CCATAAGG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )
  
  def test_exact_deletion_wrong_overlap(self):
    v1 = bio.IndelVariation( 22, 'GCCTTAAG', 'GG' )
    v2 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )

  # allow this
  #def test_exact_deletion_too_far(self):
  #  v1 = bio.IndelVariation( 30, 'GCCTTAAG', 'GG' )
  #  v2 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
  #  self.assertFalse( v1.matches( v2 ) )
  #  self.assertFalse( v2.matches( v1 ) )
  
  def test_inexact_deletion_ok(self):
    v1 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    v2 = bio.IndelVariation( 20, 'GCCTGAAG', 'GG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )

  def test_exact_deletion_ok(self):
    v1 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    v2 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
  
  def test_inexact_deletion_ok(self):
    v1 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    v2 = bio.IndelVariation( 20, 'GCCTGAAG', 'GG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )

  def test_deletion_overlap_ok(self):
    v1 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    v2 = bio.IndelVariation( 21, 'CCTTAAGG', 'CG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
 
  def test_deletion_overlap_ok_2(self):
    v1 = bio.IndelVariation( 20, 'GGAAC', 'GC' )
    v2 = bio.IndelVariation( 21, 'GAACC', 'GC' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
 
  def test_short_overlap(self):
    v1 = bio.IndelVariation( 5, 'CA', 'CAA' )
    v2 = bio.IndelVariation( 6, 'AC', 'AAC' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
 
  def test_short_overlap_2(self):
    # ok if .....CAAAC
    v1 = bio.IndelVariation( 5, 'CA', 'CAA' )
    v2 = bio.IndelVariation( 7, 'AC', 'AAC' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
 
  def test_short_no_overlap_1(self):
    # ok if .....CAAAC
    v1 = bio.IndelVariation( 5, 'CA', 'CAA' )
    v2 = bio.IndelVariation( 7, 'AC', 'AGC' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )

  def test_repr_simple(self):
    v1 = bio.IndelVariation( 5, 'CA', 'CAA' )
    self.assertEqual( 'pos=5 before=CA after=CAA', v1.__repr__() )

  def test_repr_long(self):
    before = 'A' * 200
    after = 'G' * 200
    v1 = bio.IndelVariation( 5, 'CA', 'CAA' )
    self.assertTrue( len( v1.__repr__() ) < 300 )

