
import collections
import math
import zlib

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

if __name__ == '__main__':
  import doctest
  doctest.testmod()
  print "done"
