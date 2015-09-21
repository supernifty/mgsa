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
  total = len(s)
  print "building suffixes for %i..." % total
  for i in xrange(total):
    #suffixes.append( s[i:] )
    suffixes.append( buffer( s, i ) )
    if i % 10000 == 0:
      print i, "..."
  print "sorting suffixes..."
  suffixes.sort()
  print "searching longest..."
  best = ''
  best_pos = None
  for i in xrange(0, total - 1):
    cand = lcp(suffixes[i], suffixes[i+1]);
    if len(cand) > len(best):
      best = cand
      best_pos = (total - len(suffixes[i]), total - len(suffixes[i+1]))
      print "new best: %i" % len(best)
  return best, best_pos

s = []
name = ''
for line in sys.stdin:
  if line.startswith('>'):
    if len(s) > 0:
      best = lrs( ''.join(s) )
      print best
      print name, len(best)
      s = []
    name = line[1:]
    continue
  s.append( line.strip() )

best, best_pos = lrs( ''.join(s) )
print best
print len(best), best_pos
