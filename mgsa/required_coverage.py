
# given a poisson dist, what coverage is required to get a uniform coverage of n

import math
import numpy as np
import random

def poisson( genome_len, read_len, depth, trials=100 ):
  results = []
  for trial in xrange(trials):
    genome = [0] * genome_len
    count = 0
    bases = 0
    next_trial = False
    passed = set()
    while not next_trial:
      p = random.randint(-read_len + 1, genome_len - 1)
      count += 1
      for x in xrange(max(0, p), min(p+read_len, genome_len)):
        genome[x] += 1 
        bases += 1
        if genome[x] >= depth and x not in passed:
          passed.add( x )
          if len(passed) == genome_len:
            results.append( bases )
            next_trial = True
            break
  return results

def no_gaps( genome_len, read_len, reads ):
  return math.exp( -reads * math.exp( -reads * 1. * read_len / genome_len ) )

gl = 4600000
rl = 100
depth = 2
trials = 10
bases = poisson( gl, rl, depth, trials )
coverage = [ 1. * b / gl for b in bases ] #[ 1. * y * rl / gl for y in x ]
print 'max %.1f min %.1f avg %.1f sd %.1f' % ( max(coverage), min(coverage), np.mean(coverage), np.std(coverage) )
average_reads = np.mean(bases) / rl
print 'no gaps for %.1f: %.1f' % ( average_reads, no_gaps( gl, rl, average_reads ) )
