
import sys
import helpers

if len(sys.argv) != 2:
  print "Usage: %s k" % sys.argv[0]
  sys.exit(1)

k = int(sys.argv[1])
elements = set()
processed = 0
last = 0
duplicates = 0

print "finding duplicates in fasta with k=%i" % k
print "reading fasta"
unprocessed = ''
count = 0
for line in sys.stdin:
  if line.startswith( '>' ):
    pass
  else:
    s = line.strip()
    to_process = unprocessed + s
    duplicates = helpers.repeats_update( to_process, k, elements, duplicates )
    if len(to_process) >= k:
      unprocessed = to_process[len(to_process)-k+1:] # up to k last characters
    else:
      unprocessed = to_process
    processed += len(s)
    #if count < 10:
    #  print "unprocessed:", len(unprocessed), "to_process:", len(to_process), "s:", len(s)
    if processed - last > 1e7:
      print "processed %i with %i duplicates" % ( processed, duplicates )
      last = processed
    count += 1
  
print "processed %i with %i duplicates and %i uniques" % ( processed, duplicates, len(elements) )
