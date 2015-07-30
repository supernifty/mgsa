
import collections
import math

import os

WORD_SIZE = 3 # i.e. triplets

class ReadFeature (object):
  '''features of a read'''

  def __init__( self, content, word_size = WORD_SIZE ):
    self.content = content
    self.word_size = word_size
    self.alphabet = min( 4 ** word_size, len(self.content) - word_size)

  def entropy( self ):
    if self.alphabet <= 1:
      return 0

    counts = collections.defaultdict(int)
    words = len(self.content) - self.word_size
    for i in xrange( words ):
      counts[ self.content[i:i+self.word_size] ] += 1 # count how many of each (overlapping) word is in content

    result = 0.
    max_count = max( self.alphabet, words )
    for word in counts:
      p = 1. * counts[word] / max_count
      result -= p * math.log( p, self.alphabet ) 
      #print "word", word, "result", result, "counts", counts[word], "alphabet", self.alphabet

    return result

  def complexity( self ):
    '''based on dust TODO'''
    pass

  def gc( self ):
    denom = self.content.count('G') + self.content.count('C') + self.content.count('A') + self.content.count('T') 
    if denom == 0:
      return None
    else:
      return 1. * ( self.content.count('G') + self.content.count('C') ) / denom

class MapFeature (object):
  '''features of a mapped location'''

  def __init__( self, fasta ):
    '''
      @fasta: ProbabilisticFasta
    '''
    self.fasta = fasta

  def depth( self, position ):
    '''
      number of reads covering the position
      @position: int
    '''
    counts =  self.fasta.count( position )
    return sum( counts.values() )
