
import collections
import math
import sys
import zlib

def series_from_pipeline_batch( fh, x, y, bias_report=True, column_offset=0 ):
  '''
    given a pipeline batch result file, return a list of result arrays
    @fh: file .out
    @x: e.g. insert_prob
    @y: e.g. vcf_f1
    @column_offset: use if older version of output with missing columns
  '''
  y_map = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }
  if not bias_report:
    for key in y_map:
      y_map[key] += 4
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
  y_map = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }
  if not bias_report:
    for key in y_map:
      y_map[key] += 4
  xs = []
  ys = []
  # found the line
  fields = line.strip().split(',')
  ys = [ float(x) for x in fields[y_map[y] + column_offset].split('|' ) ] # gets each y field
  for x in xrange(len(ys)):
    xs.append( 1. * x / (len(ys)-1) ) # proportion of y
  return xs, ys # e.g. [ 0.00 0.25 0.5 0.75 ], [ values ]

 
def series_from_pipeline_result( fh, y, bias_report=False, item=1, column_offset=0 ):
  '''
    given a pipeline batch result file, look at each row starting from item and return [xs, ys], where
    xs is the proportion; ys is the field specified by y 
  '''
  y_map = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }
  if not bias_report:
    for key in y_map:
      y_map[key] += 4
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
  y_map = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }
  if not bias_report:
    for key in y_map:
      y_map[key] += 4
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
  y_map = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }
  if not bias_report:
    for key in y_map:
      y_map[key] += 6
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
  for key in result:
    print key
    for item in result[key]:
      print "  %s" % item
    
if __name__ == '__main__':
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

  quick_view( sys.stdin )
