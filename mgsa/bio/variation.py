#
# manage different types of variation
#

class VariationManager(object):
  '''
    manages a collection of variations
  '''
  def __init__(self):
    self.snp_list = []
    self.snp_map = {} # maps pos snp
    self.indel_list = []
    self.indel_map = {} # maps pos to indel
    self.net_diff = 0
    self.reference_to_candidate = []

  def indel( self, pos, before, after ):
    '''
      add an indel to the collection of variations
    '''
    self.indel_map[int(pos)] = len(self.indel_list)
    self.indel_list.append( IndelVariation( int(pos), before, after ) )
    self.reference_to_candidate.append( pos + self.net_diff ) # don't include most recent indel
    diff = len(after) - len(before)
    self.net_diff += diff

  def find_indel_match( self, indel ):
    '''
      look for a matching indel in this list of indels
    '''
    floor = self.bisect_indel( self.indel_list, indel.pos )
    if indel.matches( self.indel_list[floor] ):
      return self.indel_list[floor]
    if floor + 1 < len(self.indel_list) and indel.matches( self.indel_list[floor + 1] ):
      return self.indel_list[floor+1]
    return None

  def bisect_indel(self, pos_list, value):
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

  def bisect_candidate_position(self, candidate_pos):
    '''returns an index that is less than or equal to value'''
    start = 0 # inclusive
    end = len(self.reference_to_candidate) # exclusive
    #print self.reference_to_candidate, candidate_pos
    while end - start > 1:
      mid = start + ( end - start ) / 2
      if self.reference_to_candidate[mid] == candidate_pos: # match
        return mid
      elif self.reference_to_candidate[mid] < candidate_pos: # mid+
        start = mid
      else: # pos > value
        end = mid
    if start == 0 and self.reference_to_candidate[0] > candidate_pos:
      return None
    else:
      return start 

  def candidate_position_to_reference_position(self, candidate_position):
    '''
      returns the corresponding reference position, plus any offset back to the reference position
      e.g. ref = GA; candidate = GTTA; indel should be { pos 2 before GA after GTTA }
      0 -> (0, 0)
      1 -> (0, 1)
      2 -> (0, 2)
      3 -> (1, 0)
    '''
    # look for next indel
    ceil = self.bisect_candidate_position( candidate_position )
    if ceil is None: # before any indels
      return (candidate_position, 0, 0)
    else:
      indel_length = len( self.indel_list[ceil].after ) - len( self.indel_list[ceil].before )
      past = candidate_position - self.reference_to_candidate[ceil]
      if past >= indel_length: # past last indel
        return (self.indel_list[ceil].pos + past - indel_length, 0, 0 )
      else: # inside indel
        return (self.indel_list[ceil].pos - 1, past + 1, indel_length)
      
class SNPVariation(object):
  pass

class IndelVariation(object):
  '''
    represents an indel
  '''

  def __init__( self, pos, before, after ):
    self.pos = pos
    self.before = before
    self.after = after

  def matches( self, other ):
    '''
      first look for an exact match of variation, then an overlapping match
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
