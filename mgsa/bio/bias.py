
'''
  classes to calculate biasses
'''
import datetime
import os
import sys

import bio
import sam

class Mappability (object):
  def __init__( self, pos_generator, mapper_function, read_length, min_pos=0, max_pos=-1, variation=None, error=None, log=bio.log_stderr, remove_files=True, calculate_bias=False ):
    if max_pos == -1:
      max_pos = len(pos_generator)
    self.accuracy = []
    self.summary = { 'max_accuracy': 0, 'min_accuracy': 1., 'max_incorrect': 0, 'max_unmapped': 0, 'genome_length': max_pos - min_pos, 'read_length': read_length }
    if calculate_bias:
      self.bias = []
      self.summary['max_bias'] = -1
      self.summary['min_bias'] = 1
    target_fq_fn = 'tmp_mappability_%s.fastq' % datetime.datetime.now().strftime("%H%M%S") 
    target_sam_fn = 'tmp_mappability_%s.sam' %  datetime.datetime.now().strftime("%H%M%S")
    is_first = True
    for pos in xrange(min_pos, max_pos):
      with open( target_fq_fn, 'w' ) as target_fq_fh: # make a temp file
        pos_generator.write( target_fq_fh, pos, read_length, variation=variation, error=error ) # write reads to it
      mapper_function( target_fq_fn, target_sam_fn, index=is_first ) # map the reads to the reference
      is_first = False
      # correct vs wrong/unmapped
      # generate stats 
      evaluator = sam.SamAccuracyEvaluator( open( target_sam_fn, 'r' ) )
      denominator = evaluator.stats['correct'] + evaluator.stats['incorrect'] + evaluator.stats['unmapped']
      if denominator == 0:
        alt_accuracy = 0
      else:
        alt_accuracy = 1. * evaluator.stats['correct'] / denominator
      self.accuracy.append( alt_accuracy )
      self.summary['max_accuracy'] = max( self.summary['max_accuracy'], alt_accuracy )
      self.summary['min_accuracy'] = min( self.summary['min_accuracy'], alt_accuracy )
      self.summary['max_incorrect'] = max( self.summary['max_incorrect'], evaluator.stats['incorrect'] )
      self.summary['max_unmapped'] = max( self.summary['max_unmapped'], evaluator.stats['unmapped'] )
      # bias
      if calculate_bias:
        with open( target_fq_fn, 'w' ) as target_fq_fh: # make a temp file
          pos_generator.write( target_fq_fh, pos, read_length, variation=None, error=error ) # write reads to it
        mapper_function( target_fq_fn, target_sam_fn, index=is_first ) # map the reads to the reference
        evaluator = sam.SamAccuracyEvaluator( open( target_sam_fn, 'r' ) )
        denominator = evaluator.stats['correct'] + evaluator.stats['incorrect'] + evaluator.stats['unmapped']
        if denominator == 0:
          ref_accuracy = 0
        else:
          ref_accuracy = 1. * evaluator.stats['correct'] / denominator # without variation
        bias = ref_accuracy - alt_accuracy
        self.summary['max_bias'] = max( self.summary['max_bias'], bias )
        self.summary['min_bias'] = min( self.summary['min_bias'], bias )
        self.bias.append( bias ) 
      # ticker
      if pos < 10 or int( 1. * pos / max_pos * 100 ) > int( 1. * ( pos -1 ) / max_pos * 100 ):
        log( 'Mappability pos %i (%i%%) accuracy %f' % ( pos, int( 1. * pos / max_pos * 100 ), alt_accuracy ) )
        #log( 'Mappability: incorrect: %s' % evaluator.incorrect )
    self.summary['mean_accuracy'] = float( sum(self.accuracy) / len(self.accuracy) )
    if calculate_bias:
      self.summary['mean_bias'] = float( sum(self.bias) / len(self.bias) )
    if remove_files:
      os.remove( target_fq_fn )
      os.remove( target_sam_fn )
