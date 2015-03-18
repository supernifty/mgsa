
import collections
import math

WORD_SIZE = 3

class ReadFeature (object):
  def __init__( self, content, word_size = WORD_SIZE ):
    self.content = content
    self.word_size = word_size
    self.alphabet = min( 4 ** word_size, len(self.content) - word_size)

  def entropy( self ):
    counts = collections.defaultdict(int)
    for i in xrange( len(self.content) - self.word_size ):
      counts[ self.content[i:i+self.word_size] ] += 1

    result = 0.
    for word in counts:
      result -= 1. * counts[word] / self.alphabet * math.log( 1. * counts[word] / self.alphabet, self.alphabet ) 

    return result

  def complexity( self ):
    '''based on dust TODO'''
    pass
