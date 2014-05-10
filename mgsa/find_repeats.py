
import sys

import helpers

if len(sys.argv) != 2:
  print "Usage: %s k" % sys.argv[0]
  sys.exit(1)

k = int(sys.argv[1])

print "reading fasta"
s = ''
for line in sys.stdin:
  if line.startswith( '>' ):
    pass
  else:
    s += line.strip()
  
print "finding duplicates in fasta with k=%i, n=%i" % ( k, len(s) )
duplicates = helpers.repeats( s, k )
print "%i duplicates" % duplicates
