
# compute some coverage stats for an xmfa file

import collections
import sys

count = 0
covered = collections.defaultdict(int)
unique = collections.defaultdict(int)
for line in sys.stdin:
  if count == 0:
    count += 1
    continue
  fields = line.split()
  bases = collections.defaultdict(int)
  for i in xrange(0, len(fields), 2):
    base_count = int(fields[i+1]) - int(fields[i])
    covered[i / 2] += base_count + 1
    bases[i / 2] += base_count
  if bases[0] == 0:
    unique[1] += bases[1] + 1
  if bases[1] == 0:
    unique[0] += bases[0] + 1
  count += 1
  if count % 1000 == 0:
    print "processed %i lines..." % count

for key in covered.keys():
  print "total", key, covered[key]

for key in unique.keys():
  print "unique", key, unique[key]
