
import collections
import itertools
import sys


def count_homopolymers( fh ):
  s = []
  print "building..."
  for line in fh:
    if line.startswith( '>' ):
      continue
    s.append( line.strip() )
  print "counting..."
  runs = collections.defaultdict(int)
  best = collections.defaultdict(int)
  last = None
  for c in ''.join(s):
    runs[c] += 1
    best[c] = max(best[c], runs[c])
    if last is not None and c != last:
      runs[last] = 0
    last = c
  return best
  

if __name__ == '__main__':
  print count_homopolymers( sys.stdin )
