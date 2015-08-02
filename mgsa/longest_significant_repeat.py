# experimentally measures the longest repeated sequence expected by chance

import argparse
from itertools import imap, izip, starmap, tee
from os.path   import commonprefix
import random

a = ('a', 'c', 'g', 't')
r = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}

def pairwise(iterable): # itertools recipe
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def longest_duplicate_small(data):
    suffixes = sorted(data[i:] for i in xrange(len(data))) # O(n*n) in memory
    return max(imap(commonprefix, pairwise(suffixes)), key=len)

def longest_duplicate_buffer(data):
    n = len(data)
    sa = sorted(xrange(n), key=lambda i: buffer(data, i)) # suffix array
    def lcp_item(i, j):  # find longest common prefix array item
        start = i
        while i < n and data[i] == data[i + j - start]:
            i += 1
        return i - start, start
    size, start = max(starmap(lcp_item, pairwise(sa)), key=lambda x: x[0])
    return data[start:start + size]

def reverse_complement( s ):
  return ''.join( [ r[c] for c in reversed(s) ] )

def find_repeat( l, p, rc=False ):
  '''
    l is the string length
    p is the percentile
  '''
  r = []
  for trial in xrange(100):
    s = ''.join( [ random.choice(a) for _ in xrange(l) ] )
    if rc:
      s = s + '*' + reverse_complement( s )
    r.append( len( longest_duplicate_buffer( s ) ) )
  r.sort()
  #print r
  return r[p]


    
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='random longest significant repeat')
  parser.add_argument('--reverse_complement', dest='reverse_complement', type=bool, default=False, help='include reverse complement')
  args = parser.parse_args() 
  for l in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384):
    print l, find_repeat( l, 95, rc=args.reverse_complement ) # prints length and longest repeat at 95th percentile

