
import argparse
import collections
import numpy
import sys

import bio
import config

parser = argparse.ArgumentParser(description='Analyze FASTA')
parser.add_argument('fasta', metavar='fasta', help='fasta file to analyze')
parser.add_argument('--buckets', metavar='buckets', type=int, default=10, help='number of buckets')

args = parser.parse_args()
fs = bio.FastaStats( open( args.fasta, 'r' ) )

# gc
buckets = numpy.linspace(0, 1, args.buckets + 1)
gc_hist = bio.bucket( filter( None, fs.stats['gc'] ), buckets )
entropy_hist = bio.bucket( filter( None, fs.stats['entropy'] ), buckets )
total_gc = sum(gc_hist)
total_entropy = sum(entropy_hist)

print '========== GC content =========='
print 'GC %%:     %s' % '\t'.join( [ '%.2f' % bucket for bucket in buckets ] )
print 'dist:     %s' % '\t'.join( [ '%.1f' % ( 100. * x / total_gc ) for x in gc_hist ] )

# entropy
print '\n========== Entropy =========='
print 'Entropy:  %s' % '\t'.join( [ '%.2f' % bucket for bucket in buckets ] )
print 'dist:     %s' % '\t'.join( [ '%.1f' % ( 100. * x / total_entropy ) for x in entropy_hist ] )


