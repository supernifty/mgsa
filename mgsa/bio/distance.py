
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

  def simple_indel( self ):
    '''
      find a single indel in a and b and return its length
    '''
    if len( self.a ) > len( self.b ):
      longer = self.a
      shorter = self.b
      sign = -1
    elif len( self.a ) < len( self.b ):
      longer = self.b
      shorter = self.a
      sign = 1
    else:
      return 0
    start_break_found = False
    for start_break in xrange(len(longer)):
      if start_break >= len(shorter): # insertion on end of longer
        size = len(longer) - len(shorter)
        return size * sign
      if longer[start_break] != shorter[start_break]:
        start_break_found = True
        break
    if not start_break_found:
      return 0

    for end_break in xrange(-1, -len(longer) + start_break - 1, -1):
      if -end_break > len(shorter): # insertion at start of shorter
        size = len(longer) - len(shorter)
        return size * sign
      if -end_break > len(shorter) - start_break or longer[end_break] != shorter[end_break]: # insertion at start of difference, or mismatch found
        longer_end_pos = len(longer) + end_break
        size = longer_end_pos - start_break + 1
        return size * sign
