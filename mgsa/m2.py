import sys

for line in sys.stdin:
  f = line.strip().split('\t')
  for i, v in enumerate(f[1:]):
    print 'x %s,y %i,z %s' % ( f[0], i, v )
