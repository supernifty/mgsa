
import math

import features
import helpers

# describes evidence for a variation in terms of read depth, etc

class BreakpointPrefixEvidence (object):
  
  def __init__( self, fasta ):
    '''
      @fasta: ProbabilisticFasta
    '''
    self.fasta = fasta

  def evidence( pos, prior=0.5 ):
    '''
      @return a probability that the given artifact is present at this position. Possible result is 0..1 for probability
    '''
    return prior

class LinearShapeDifference (object):
  '''
    look for a specified, absolute linear shape
  '''
  def __init__( self, start_pos, start_depth, end_pos, end_depth ):
    self.start_pos = start_pos
    self.start_depth = start_depth
    self.end_pos = end_pos
    self.end_depth = end_depth
    # calculate max error
    self.max_error = 0
    for x_diff in xrange( start_pos, end_pos ):
        x_percent = 1. * ( x_diff - self.start_pos ) / ( self.end_pos - self.start_pos )
        expected_depth = self.start_depth + x_percent * ( self.end_depth - self.start_depth )
        self.max_error += ( 1. * start_depth - expected_depth ) ** 2
    self.max_error = math.sqrt( self.max_error / ( end_pos - start_pos ) )

  def evidence( self, fasta, pos, mean_depth=1., prior=0.5 ):
    feature = features.MapFeature( fasta )
    samples = 0
    error = 0.
    for x_diff in xrange( self.start_pos, self.end_pos + 1 ):
      current_pos = pos + x_diff
      if current_pos >= 0 and current_pos < fasta.length:
        x_percent = 1. * ( x_diff - self.start_pos ) / ( self.end_pos - self.start_pos )
        expected_depth = self.start_depth + x_percent * ( self.end_depth - self.start_depth )
        error += ( 1. * feature.depth( current_pos ) / mean_depth - expected_depth ) ** 2
        samples += 1
    if samples == 0:
      return prior
    average_error = math.sqrt( error / samples )
    if average_error > self.max_error:
      return 0
    else:
      return 1 - ( average_error / self.max_error )
    # convert to probability
    #return 1 - ( helpers.logistic( average_error ) - 0.5 ) * 2 # convert from 0.5->1.0 --> 0.0 --> 1.0 --> 1.0 --> 0.0

# predefined evidence
BREAKPOINT_PREFIX_READ_DEPTH_SHAPE = LinearShapeDifference( -30, 1.0, 0, 0.7  )
BREAKPOINT_SUFFIX_READ_DEPTH_SHAPE = LinearShapeDifference( 0, 0.5, 50, 1.0 )

