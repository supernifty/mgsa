# experimentally measures the longest repeated sequence expected by chance

from itertools import imap, izip, starmap, tee
from os.path   import commonprefix
import random

a = ('a', 'c', 'g', 't')

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

def find_repeat( l, p ):
  '''
    l is the string length
    p is the percentile
  '''
  r = []
  for trial in xrange(100):
    s = ''.join( [ random.choice(a) for _ in xrange(l) ] )
    r.append( len( longest_duplicate_buffer( s ) ) )
  r.sort()
  #print r
  return r[p]


    
if __name__ == "__main__":
  for l in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384):
    print l, find_repeat( l, 95 ) # prints length and longest repeat at 95th percentile
  
