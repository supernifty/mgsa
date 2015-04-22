
def bucket( values, buckets ):
  '''
    puts values into buckets
    speed this up by sorting values
    @values: list of values
    @buckets: list of bucket separators with each being the minimum value that can go in that bucket
  '''
  result = [ 0 ] * len(buckets)
  current_bucket = 0
  for value in sorted( values ):
    while current_bucket < len(buckets) - 1 and value >= buckets[current_bucket + 1]:
      current_bucket += 1
    result[current_bucket] += 1
  return result

def sign_bucket( values ):
  '''
    returns proportions of +ve, -ve, 0 and the proportion of +ve to -ve
  '''
  result = {}
  result['positive'] = 1. * len( [ x for x in values if x > 0 ] ) / len( values )
  result['negative'] = 1. * len( [ x for x in values if x < 0 ] ) / len( values )
  result['zero'] = 1. * len( [ x for x in values if x == 0 ] ) / len( values )
  if result['positive'] + result['negative'] > 0:
    result['positive_proportion'] = result['positive'] / ( result['positive'] + result['negative'] )
  return result
