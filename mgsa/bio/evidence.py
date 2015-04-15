
import math

import bio
import features
import helpers

# describes evidence for a variation in terms of read depth, etc

class LinearShapeDifference (object):
  '''
    look for a specified, absolute linear shape
  '''
  def __init__( self, start_pos, start_depth, end_pos, end_depth, log=bio.log_quiet ):
    self.start_pos = start_pos
    self.start_depth = start_depth
    self.end_pos = end_pos
    self.end_depth = end_depth
    # calculate max error
    self.max_error = {}
    self.log = log

  def max_error_list( self, mean_depth ):
    if mean_depth not in self.max_error:
      result = []
      no_evidence_depth = 1.
      for x_diff in xrange( self.start_pos, self.end_pos ):
        x_percent = 1. * ( x_diff - self.start_pos ) / ( self.end_pos - self.start_pos )
        expected_depth = self.start_depth + x_percent * ( self.end_depth - self.start_depth )
        expected_probability = helpers.binomial_prob( int( mean_depth * expected_depth ), int( mean_depth * 100 ), 100 )
        result.append( ( no_evidence_depth - expected_depth ) / expected_probability )
      self.max_error[mean_depth] = result
      self.log( 'max error: %s' % result )
    return self.max_error[mean_depth]
      
  def max_errors( self, mean_depth ):
    result = self.max_error_list( mean_depth )
    return math.sqrt( sum( [ x ** 2 for x in result ] ) / len( result ) )
 
  def errors( self, fasta, pos, mean_depth=1. ):
    '''
      return list of all error terms
    '''
    feature = features.MapFeature( fasta )
    samples = 0
    result = []
    for x_diff in xrange( self.start_pos, self.end_pos + 1 ):
      current_pos = pos + x_diff
      if current_pos >= 0 and current_pos < fasta.length:
        x_percent = 1. * ( x_diff - self.start_pos ) / ( self.end_pos - self.start_pos )
        expected_depth = self.start_depth + x_percent * ( self.end_depth - self.start_depth )
        actual_depth = feature.depth( current_pos )
        expected_probability = helpers.binomial_prob( int( mean_depth * expected_depth ), int( mean_depth * 100 ), 100 )
        result.append( ( 1. * actual_depth / mean_depth - expected_depth ) / expected_probability )
    self.log( 'result: %s' % result )
    return result
 
  def evidence( self, fasta, pos, mean_depth=1., prior=0.5 ):
    errors = self.errors( fasta, pos, mean_depth )
    average_error = math.sqrt( sum( [ x ** 2 for x in errors ] ) / len( errors ) )
    return max( 0, 1 - ( average_error / self.max_errors( mean_depth ) ) )
    # convert to probability
    #return 1 - ( helpers.logistic( average_error ) - 0.5 ) * 2 # convert from 0.5->1.0 --> 0.0 --> 1.0 --> 1.0 --> 0.0

  def shift( self, amount ):
    '''
      return a LinearShapeDifference with pos shifted
    '''
    return LinearShapeDifference( self.start_pos + amount, self.start_depth, self.end_pos + amount, self.end_depth )

class MultipleLinearShapeDifference (object):
  def __init__( self, shapes, log=bio.log_quiet ):
    '''
      @shapes: list of linear shapes
    '''
    self.shapes = shapes
    self.log = log
    self.max_error = {}

  def max_errors( self, depth ):
    if depth not in self.max_error:
      result = []
      for shape in self.shapes:
        shape_max_errors = shape.max_error_list( depth )
        result.extend( shape_max_errors  )
      self.max_error[depth] = math.sqrt( sum( [ x ** 2 for x in result ] ) / len( result ) )
    return self.max_error[depth]

  def evidence( self, fasta, pos, mean_depth=1., prior=0.5 ):
    errors = []
    for shape in self.shapes:
      errors.extend( shape.errors( fasta, pos, mean_depth ) )
    #self.log( 'errors: %s' % errors )
    average_error = math.sqrt( sum( [ x ** 2 for x in errors ] ) / len( errors ) )
    #self.log( 'average error: %.2f max error: %.2f' % ( average_error, self.max_error ) )
    return max( 0, 1 - ( average_error / self.max_errors( mean_depth ) ) )

# predefined evidence
# note that these #s are a function of the aligner algorithm (bwa) and read length (100) but not coverage (I believe)
BREAKPOINT_PREFIX_READ_DEPTH_SHAPE = LinearShapeDifference( -30, 1.0, 0, 0.7  )
BREAKPOINT_SUFFIX_READ_DEPTH_SHAPE = LinearShapeDifference( 0, 0.5, 50, 1.0 )

class LongDeletionDetector(object):
  '''
    read depth signature for a deletion that's long enough such that no read spans it (for default config bwa 100coverage 100rl this is ~45)
  '''
  def __init__( self, deletion_size, vcf_writer=None, threshold=0.5, reference_fasta=None, log=bio.log_quiet ):
    self.deletion_size = deletion_size
    self.shape = MultipleLinearShapeDifference( [ BREAKPOINT_PREFIX_READ_DEPTH_SHAPE, LinearShapeDifference( 1, 0, deletion_size, 0 ), BREAKPOINT_SUFFIX_READ_DEPTH_SHAPE.shift( deletion_size + 1 ) ], log=log  )
    self.vcf_writer = vcf_writer
    self.threshold = threshold
    self.last = [ None, None ] # previous 2 results
    self.reference_fasta = reference_fasta # to get deleted bases
    self.log = log

  def evidence( self, fasta, pos, mean_depth=1., prior=0.5 ):
    '''
    '''
    result = self.shape.evidence( fasta, pos, mean_depth, prior )
    # write to vcf if this one is below previous and previous is at least previous[-2]
    if self.threshold is not None and self.vcf_writer is not None:
      if self.last[0] is not None and self.last[1] is not None and self.last[1] >= self.last[0] and result < self.last[1] and self.reference_fasta.read_to( pos + self.deletion_size ):
        self.log( 'potential long deletion at %i: evidence %.2f' % ( pos, self.last[1] ) )
        if self.last[1] >= self.threshold:
          before = '%s%s' % ( fasta.consensus( start = pos - 1, end = pos - 1 + 1 ), self.reference_fasta.fragment( pos, pos + self.deletion_size ) ) # ref has the extra bases
          after = '%s' % ( fasta.consensus( start=pos - 1, end=pos - 1 + 1 ) ) #, fasta.consensus( start=pos + 1 + self.deletion_size, end=pos + 1 + self.deletion_size + 1 ) ) # cand has the missing bases
          self.vcf_writer.indel( pos - 1, before, after, self.last[1] )
      # update previous
      self.last[0] = self.last[1]
      self.last[1] = result
      

    return result
