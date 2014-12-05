#
# manage different types of variation
#

class VariationManager(object):
  def __init__(self):
    self.snp_list = []
    self.snp_map = {} # maps pos snp
    self.indel_list = []
    self.indel_map = {} # maps pos to indel

  def indel( self, pos, before, after ):
    self.indel_map[int(pos)] = len(self.indel_list)
    self.indel_list.append( IndelVariation( int(pos), before, after ) )

  def find_indel_match( self, indel ):
    '''
      look for a matching indel in this list of indels
    '''
    floor = self.bisect( self.indel_list, indel.pos )
    if indel.matches( self.indel_list[floor] ):
      return self.indel_list[floor]
    if floor + 1 < len(self.indel_list) and indel.matches( self.indel_list[floor + 1] ):
      return self.indel_list[floor+1]
    return None

  def bisect(self, pos_list, value):
    '''returns an index that is less than or equal to value'''
    start = 0 # inclusive
    end = len(pos_list) # exclusive
    while end - start > 1:
      mid = start + ( end - start ) / 2
      if pos_list[mid].pos == value:
        return mid
      elif pos_list[mid].pos < value:
        start = mid
      else: # pos > value
        end = mid
      #print "bisect: start: %i end: %i" % ( start, end )
    return start

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
      self_overlap = self.after[pos_diff:]
      other_overlap = other.after[:-pos_diff]
    else: # other is first; pos_diff is -ve
      other_overlap = other.after[-pos_diff:]
      self_overlap = self.after[:pos_diff]

    #print "overlap_matches: %s %s: overlaps %s %s; result %s" % ( self, other, self_overlap, other_overlap, self_overlap == other_overlap )
    return self_overlap == other_overlap

  def __str__( self ):
    return 'pos=%i before=%s after=%s' % ( self.pos, self.before, self.after )
