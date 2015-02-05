
def bucket( values, buckets ):
  '''speed this up by sorting values'''
  result = [ 0 ] * len(buckets)
  for value in values:
    for i in xrange(len(buckets)):
      if i == len(buckets) - 1 or value >= buckets[i] and value < buckets[i+1]:
        result[i] += 1
        break
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
