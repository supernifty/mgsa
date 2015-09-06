
# find repeats with mismatches
import collections
import itertools
import sys

def compare( s1, s2, mm ):
  for i in xrange(len(s1)):
    if s1[i] != s2[i]:
      mm -= 1
      if mm < 0:
        return False
  return True
    
def is_match( a, b, s, k, mm ):
  if mm == 0:
    if a > b:
      return (b, a)
    return (a, b)
  # try the block in each position
  for x in xrange(0, mm):
    offset = x * k / (mm+1)
    if a - offset < 0 or b - offset < 0 or a - offset + k > len(s) or b - offset + k > len(s):
      continue
    s1 = s[a - offset:a - offset + k]
    s2 = s[b - offset:b - offset + k]
    if compare(s1, s2, mm):
      if a > b:
        return (b, a)
      return (a, b)
  return None

def find_mm( s, k, mm=0 ):
  split_size = k / ( mm + 1 )
  # make kmers
  print "building %i %i-mers..." % ( len(s) - split_size + 1, split_size )
  kmers = collections.defaultdict(set)
  for i in xrange(0, len(s) - split_size + 1):
    kmers[s[i:i+split_size]].add( i )
    if i % 100000 == 0:
      print "%i added..." % i

  print "searching %i %i-mers for up to %i mismatches..." % ( len(kmers), k, mm )
  found = set()
  for idx, kmer in enumerate(kmers):
    if len(kmers[kmer]) > 1:
      for a, b in itertools.combinations( kmers[kmer], 2 ):
        #print "trying", a, b
        result = is_match( a, b, s, k, mm )
        if result is not None:
          found.add( result )
    if idx % 1000000 == 0:
      print "%i searched..." % idx
  return found

if __name__ == '__main__':
  k = int(sys.argv[1])
  mm = int(sys.argv[2])
  s = []
  for line in sys.stdin:
    if line.startswith( '>' ):
      continue
    s.append( line.strip() )
  f = ''.join( s )
  result = find_mm( f, k=k, mm=mm )
  print "-" * 60
  print "%i matches" % len(result)
  max_show = 10
  print "Showing up to 10 examples..."
  for r in result:
    print '%09i:%s\n%09i:%s\n' % ( r[0], f[r[0]:r[0]+k], r[1], f[r[1]:r[1]+k] )
    max_show -= 1
    if max_show < 0:
      break
