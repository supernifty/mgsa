
class Distance(object):
  def __init__( self, a, b ):
    self.a = a
    self.b = b

  def hamming( self ):
    if len( self.a ) != len( self.b ):
      raise Exception( "Lengths differ" )
    result = 0
    for i in xrange( len(self.a) ):
      if self.a[i] != self.b[i]:
        result += 1
    return result
