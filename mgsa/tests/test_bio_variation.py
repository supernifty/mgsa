
import unittest
import StringIO

import bio

class TestIndelVariation( unittest.TestCase ):

  def setUp(self):
    pass

  def test_exact_ok(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
  
  def test_inexact_ok(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 20, 'GG', 'GCCTGAAG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )

  def test_overlap(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 21, 'CG', 'CCTTAAGG' )
    self.assertTrue( v1.matches( v2 ) )
    self.assertTrue( v2.matches( v1 ) )
  
  def test_no_overlap(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 22, 'CG', 'CCTTAAGG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )
  
  def test_no_overlap_content(self):
    v1 = bio.IndelVariation( 20, 'GG', 'GCCTTAAG' )
    v2 = bio.IndelVariation( 21, 'CG', 'CCATAAGG' )
    self.assertFalse( v1.matches( v2 ) )
    self.assertFalse( v2.matches( v1 ) )
  
