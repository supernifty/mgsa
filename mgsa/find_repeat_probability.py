
import sys

import helpers

if len(sys.argv) != 3:
  print "Usage: %s n k" % sys.argv[0]
  sys.exit(1)

n = int(sys.argv[1])
k = int(sys.argv[2])

print helpers.chance_of_no_repeats( n, k )
