# a simple implementation to find the longest repeat in a string
import sys

def lcp(a, b):
  for i in xrange(min(len(a), len(b))):
    if a[i] != b[i]:
      break
  return a[:i]

# build suffixes, sort them
def lrs(s):
  suffixes = []
  print "building suffixes for %i..." % len(s)
  for i in xrange(len(s)):
    #suffixes.append( s[i:] )
    suffixes.append( buffer( s, i ) )
    if i % 10000 == 0:
      print i, "..."
  print "sorting suffixes..."
  suffixes.sort()
  print "searching longest..."
  best = ''
  for i in xrange(0, len(s) - 1):
    cand = lcp(suffixes[i], suffixes[i+1]);
    if len(cand) > len(best):
      best = cand
      print "new best: %i" % len(best)
  return best

s = []
for line in sys.stdin:
  if line.startswith('>'):
    continue
  s.append( line.strip() )

best = lrs( ''.join(s) )
print best
print len(best)
