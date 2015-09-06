# a simple implementation to find the longest repeat in a string
import collections
import sys

def extend(a, b, k, s, mm=0):
  # extend maximum distance with remaining mismatches
  lhs_mm = mm
  lhs_dist = []
  for offset in xrange(-1, max(-a, -b), -1):
    if s[a+offset] != s[b+offset]:
      lhs_dist.append( -offset )
      lhs_mm -= 1
      if lhs_mm == 0:
        break
  rhs_mm = mm
  rhs_dist = []
  for offset in xrange(1, min( len(s) - a, len(s) - b)):
    if s[a+k+offset] != s[b+k+offset]:
      rhs_dist.append( offset )
      rhs_mm -= 1
      if rhs_mm == 0:
        break
      
  return a[:i]

# build suffixes, sort them
def lrs(s, k, mm):
  kmers = collections.defaultdict(set)
  print "building kmers for %i..." % len(s)
  for i in xrange(len(s)-k):
    kmers[s[i:i+k].append(i)
    if i % 10000 == 0:
      print i, "..."
  best = ''
  print "%i k-mers" % len(kmers)
  for kmer in kmers:
    if len(kmers[kmer]) > 0:
      for a, b in itertools.combinations( kmers[kmer], 2 ):
        cand = extend( a, b, k, s, mm )
        if len(cand) > len(best):
          best = cand
          print "new best:", len(best):
          
  return best

# build kmers, compare
def lrs_mm(s, mm=0):
  pass

s = []
for line in sys.stdin:
  if line.startswith('>'):
    continue
  s.append( line.strip() )

if len(sys.argv) < 2:
  print "Usage: %s k mm" % sys.argv[0]
  print "Note: set k to expected-len / mm"
  sys.exit(1)


k = int(sys.argv[1])
mm = int(sys.argv[2])
best = lrs( ''.join(s), k, mm )
print best
print len(best)
