#
# manage different types of variation
#

class SNPVariation(object):
  pass

class IndelVariation(object):
  def __init__( self, pos, before, after ):
    self.pos = pos
    self.before = before
    self.after = after

  def matches( self, other ):
    '''
      exact match of variation
    '''
    if self.pos == other.pos:
      return self.before == other.before and self.after == other.after
    
    # same length, not far apart
    elif self.might_match( other ):
      return self.overlap_matches( other )

    return False

  def might_match( self, other ):
    '''
      does it overlap with other
    '''
    return len(self.after) - len(self.before) == len(other.after) - len(other.before) and abs(other.pos - self.pos) < len(self.after) - len(self.before);

  def overlap_matches( self, other ):
    '''
      overlap between variation matches - so it is "probably" a match
      assumes other has same length and not too far away
    '''
    pos_diff = other.pos - self.pos
    if pos_diff > 0: # self is first
      self_overlap = self.after[:-pos_diff]
      other_overlap = other.after[pos_diff:]
    else: # other is first; pos_diff is -ve
      other_overlap = self.after[:pos_diff]
      self_overlap = other.after[-pos_diff:]
