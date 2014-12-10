
import unittest
import StringIO

import bio

class TestIndelVariation( unittest.TestCase ):

  def setUp(self):
    pass

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

  def test_exact_deletion_too_far(self):
    v1 = bio.IndelVariation( 30, 'GCCTTAAG', 'GG' )
    v2 = bio.IndelVariation( 20, 'GCCTTAAG', 'GG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )
  
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
 
