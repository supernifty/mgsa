
import collections
import math
import sys
import zlib

import numpy as np

import bio
# before error_bias_no_variation
#y_map_old = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }

def make_ymap(bias_report, has_no_variation=False):
  y_map = { 'unmapped': -22, 'incorrect': -21, 'read_precision': -20, 'read_recall': -19, 'read_f1': -18, 'vcf_tp': -17, 'vcf_fp': -16, 'vcf_fn': -15, 'vcf_precision': -14, 'vcf_recall': -13, 'vcf_f1': -12, 'vcf_covered': -11, 'vcf_bucket_tp': -10, 'vcf_bucket_fp': -9, 'vcf_bucket_fn': -8, 'reference_bias': -7, 'error_bias': -6, 'unmapped_variations': -5, 'total_variations': -4, 'mean_reference': -3, 'mean_error': -2, 'error_bias_no_variation': -1 }
  if not bias_report:
    for key in y_map:
      y_map[key] += 6
  if not has_no_variation:
    for key in y_map:
      y_map[key] += 1
  return y_map

def series_from_pipeline_batch( fh, x, y, bias_report=True, column_offset=0, has_no_variation=False ):
  '''
    given a pipeline batch result file, return a list of result arrays
    @fh: file .out
    @x: e.g. insert_prob
    @y: e.g. vcf_f1
    @column_offset: use if older version of output with missing columns
  '''
  y_map = make_ymap(bias_report, has_no_variation)
  xs = []
  ys = []
  for line in fh:
    if line.startswith( 'command' ) or line.startswith( '#' ):
      continue
    fields = line.strip().split(',')
    for field in fields:
      field = field.strip()
      if field.startswith( x ):
        fieldname, value = field.split(' ')
        xs.append( float(value) )
        break
    ys.append( float( fields[y_map[y] + column_offset] ) )
  return xs, ys

def series_from_line( line, y, bias_report=False, column_offset=0 ):
  y_map = make_ymap(bias_report)
  xs = []
  ys = []
  # found the line
  fields = line.strip().split(',')
  ys = [ float(x) for x in fields[y_map[y] + column_offset].split('|' ) ] # gets each y field
  for x in xrange(len(ys)):
    xs.append( 1. * x / (len(ys)-1) ) # proportion of y
  return xs, ys # e.g. [ 0.00 0.25 0.5 0.75 ], [ values ]

 
def series_from_pipeline_result( fh, y, bias_report=False, item=1, column_offset=0, has_no_variation=False ):
  '''
    given a pipeline batch result file, look at each row starting from item and return [xs, ys], where
    xs is the proportion; ys is the field specified by y 
  '''
  y_map = make_ymap(bias_report, has_no_variation)
  xs = []
  ys = []
  for line in fh:
    if line.startswith( 'command' ) or line.startswith( '#' ):
      continue
    if item > 1:
      item -= 1
      continue
    # found the line
    fields = line.strip().split(',')
    ys = [ float(x) for x in fields[y_map[y] + column_offset].split('|' ) ] # gets each y field
    for x in xrange(len(ys)):
      xs.append( 1. * x / (len(ys)-1) ) # proportion of y
    return xs, ys # e.g. [ 0.00 0.25 0.5 0.75 ], [ values ]

def item_from_pipeline_result( fh, y, bias_report=False, item=1, column_offset=0 ):
  y_map = make_ymap(bias_report)
  for line in fh:
    if line.startswith( 'command' ) or line.startswith( '#' ):
      continue
    if item > 1:
      item -= 1
      continue
    fields = line.strip().split(',')
    return float(fields[y_map[y] + column_offset])

def add_f1_to_batch_results(h):
  '''
    given batch file from pipeline_batch, add f1-score as a new column
  '''
  for line in h:
    if line.startswith( '#' ):
      sys.stdout.write( line )
    else:
      fields = line.strip().split( ',' )
      precision = float( fields[-2] ) / 100.
      recall = float( fields[-1] ) / 100.
      if precision != 0 or recall != 0:
        f1 = 2 * ( precision * recall ) / ( precision + recall ) * 100
      else:
        f1 = 0
      sys.stdout.write( '%s,%f\n' % ( line.strip(), f1 ) )

def calculate_tpr():
  pass # TODO

def repeats( s, k ):
  '''
    count repeats of length k in string s
>>> repeats( 'abbb', 2 )
1
>>> repeats( 'ababab', 3 )
2
  '''
  elements = set()
  duplicates = 0
  for i in xrange(0, len(s) - k + 1):
    #element = s[i:i+k]
    element = buffer( s, i, k )
    if element not in elements:
      elements.add( element )
    else:
      duplicates += 1
  
  return duplicates

def repeats_update( s, k, elements, duplicates ):
  '''
    count repeats of length k in string s
>>> repeats( 'abbb', 2 )
1
>>> repeats( 'ababab', 3 )
2
  '''
  for i in xrange(0, len(s) - k + 1):
    #element = s[i:i+k]
    #element = buffer( s, i, k )
    element = to_binary( buffer( s, i, k ) )
    if element not in elements:
      elements.add( element )
    else:
      duplicates += 1
  
  return duplicates

def to_binary( sequence ):
   return zlib.compress(sequence) 

def combinations_choose_without_order( n, k ):
  '''
>>> combinations_choose_without_order( 5, 3 )
10
  '''
  return combinations_choose_with_order( n, k ) / math.factorial( k )

def combinations_choose_with_order( n, k ):
  '''
>>> combinations_choose_with_order( 5, 3 ) # 5 * 4 * 3
60
  '''
  result = n
  base = max( 1, n - k )
  for mult in range( n - 1, base, -1 ):
    result *= mult
  #return math.factorial( n ) / math.factorial( n - k )
  return result

def chance_of_no_repeats( n, k ):
  '''
>>> chance_of_no_repeats( 5, 3 )
0.95361328125
  '''
  sequences = 4 ** k # total # of different sequences
  trials = n - k + 1 # need to choose this many uniquely
  unique_combinations = combinations_choose_with_order( sequences, trials ) # ways of selecting unique sequences
  total_combinations = sequences ** trials # all possible selections
  #print "total_combinations", total_combinations, "trials", trials, "this", unique_combinations
  return math.exp( math.log( unique_combinations ) - math.log( total_combinations ) )

def chance_of_repeats( n, k, r ):
  '''
    not quite right
    #>>> chance_of_at_least_repeats( 5, 3, 1 )
    #0.95361328125
  '''
  if r == 0:
    return chance_of_no_repeats( n, k )
  sequences = 4 ** k # total # of different sequences
  trials = n - k + 1 # need to choose this many uniquely
  #print "sequences", sequences, "trials", trials
  #combinations = sequences * trials # put any sequence anywhere
  #print "combinations", combinations
  #for repeats in xrange( 0, r ):
  #  combinations *= ( trials - repeats - 1 ) # same sequence in one of remaining position
  #  print "combinations", combinations
  # how many ways of ending up with these placements? 
  #combinations /= math.factorial( r + 1 )
  # place anything in the remaining
  #combinations *= sequences ** ( trials - r - 1 ) # anything in remaining
  #combinations = sequences * sequences ** ( trials - r - 1 ) * combinations_choose_with_order( trials, r )
  #combinations = sequences * math.factorial( r + 1 ) * math.factorial( trials - r - 1 )
  #combinations = binomial( trials, 1/sequences, 1-1/sequences)
  #combinations = combinations_choose_without_order( trials, r ) * 1 / sequences ** r * ( 1 - 1 / sequences ) ** ( trials - r )
  #combinations = combinations_choose_without_order( trials, r ) * 1 / sequences ** r * ( 1 - 1 / sequences ) ** ( trials - r )
  chance_one = math.exp( math.log( combinations_choose_without_order( trials, r ) ) - r * math.log( sequences ) + ( trials - r ) * math.log( 1 - 1./sequences ) )
  print chance_one
  chance_any = 1 - math.exp( sequences * math.log( 1 - chance_one ) )
  return chance_any
  #print "combinations", combinations
  #total_combinations = sequences ** trials
  #print "total_combinations", total_combinations, "selections", trials, "this", combinations
  #return math.exp( math.log( combinations ) - math.log( total_combinations ) )

def binomial( n, k ):
  return math.factorial( n ) / math.factorial( k ) / math.factorial( n - k )

def miscall_probability( coverage, wrong_prob ):
  p = 0.0
  for k in xrange( int( math.ceil( coverage / 2. ) ), coverage + 1 ):
    prob = binomial( coverage, k ) * ( wrong_prob ** k ) * ( ( 1 - wrong_prob ) ** ( coverage - k ) )
    print k, coverage, binomial( coverage, k ), prob
    p += prob
  return p

def quick_view( fh, bias_report=False ):
  '''
    quick summary of the output from an experiment
  '''
  y_map = make_ymap(bias_report)
  rev_map = dict((v,k) for k,v in y_map.iteritems())
  result = {}
  for line in fh:
    if line.startswith( '#' ): 
      continue
    fields = line.split(',')
    for index, field in enumerate(fields):
      dist = index - len(fields)
      if dist in rev_map:
        key = rev_map[dist]
        if key not in result:
          result[key] = []
        result[key].append( field.strip() )
  keys = sorted( result.keys() )
  
  for key in keys:
    print key
    for item in result[key]:
      print "  %s" % item

if __name__ == '__main__':
  # TODO really should use argparse
  if len(sys.argv) < 2:
    print "%s (quick_view|fasta_stats|fasta_width width|fasta_edit command width|fasta_filter command|fasta_find query|vcf_consensus vcf <fasta)" % sys.argv[0]
  elif sys.argv[1] == 'quick_view': # pipeline batch output
    quick_view( sys.stdin )
  elif sys.argv[1] == 'fasta_stats':
    s = bio.fasta_stats( sys.stdin )
    stats = s.stats
    bases = s.base_counter
    print "%i bases" % stats['count']
    print "Bases:\t%s" % bases
    print "Measure\t| Mean\t| SD"
    print "-------\t|-----\t|---"
    print "GC %%\t| %.1f\t| %.1f" % ( np.mean( stats['gc'] ), np.std( stats['gc'] ) )
    print "Entropy\t| %.1f\t| %.1f" % ( np.mean( stats['entropy'] ), np.std( stats['entropy'] ) )
  elif sys.argv[1] == 'fasta_width': # resize to a provided width
    bio.fasta_set_width( sys.stdin, sys.stdout, int(sys.argv[2]) )
  elif sys.argv[1] == 'fasta_edit':
    width = 50
    if len(sys.argv) > 3:
      width = int(sys.argv[3])
    bio.fasta_edit( sys.stdin, sys.stdout, commands=sys.argv[2], width=width )
  elif sys.argv[1] == 'fasta_filter':
    width = 50
    if len(sys.argv) > 3:
      width = int(sys.argv[3])
    bio.fasta_filter( sys.stdin, sys.stdout, commands=sys.argv[2] )
  elif sys.argv[1] == 'fasta_find':
    sys.stdout.write( '%s\n' % bio.fasta_find( sys.stdin, query=sys.argv[2] ) )
  elif sys.argv[1] == 'vcf_consensus':
    fr = bio.FastaReaderFromVCF( bio.FastaReader( sys.stdin ), bio.VCF( open( sys.argv[2], 'r' ) ) )
    for line in fr.items():
      sys.stdout.write( '%s\n' % line )
  elif sys.argv[1] == 'poisson': # k, m, n: probability of k items in a slot, when there are m items and n slots 
    sys.stdout.write( '%f\n' % bio.binomial_prob( int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) ) )
  else:
    print "unrecognized command"
  #import doctest
  #doctest.testmod()
  #print "done"

  #python_plot_from_pipeline_batch( open( 'pipeline_batch_read_length.out' ) )
  #add_f1_to_batch_results( sys.stdin )
  #print "10, 0.01:", miscall_probability( 10, 0.01 )
  #print "10, 0.1:", miscall_probability( 10, 0.1 )
  #print "5, 0.1:", miscall_probability( 5, 0.1 )
  #print "3, 0.1:", miscall_probability( 3, 0.1 )
  #print "4, 0.02:", miscall_probability( 4, 0.02 )

