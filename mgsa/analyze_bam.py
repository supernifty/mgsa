
import argparse
import collections
import numpy
import sys

import bio
import config

parser = argparse.ArgumentParser(description='Analyze BAM')
parser.add_argument('bam', metavar='bam', help='bam file to analyze')
parser.add_argument('--buckets', metavar='buckets', type=int, default=10, help='number of buckets')

args = parser.parse_args()
bam = bio.BamReaderExternal( config.BAM_TO_SAM, args.bam )

# gc
stats = bio.SamStats( bam )
buckets = numpy.linspace(0, 1, args.buckets + 1)
mapped_buckets = bio.bucket( stats.mapped['gc'], buckets )
unmapped_buckets = bio.bucket( stats.unmapped['gc'], buckets )
total_mapped = sum( mapped_buckets )
total_unmapped = sum( unmapped_buckets )

print '========== GC content =========='
print 'GC %%:     %s' % '\t'.join( [ '%.2f' % bucket for bucket in buckets ] )
print 'mapped:   %s' % '\t'.join( [ '%.1f' % ( 100. * x / total_mapped ) for x in mapped_buckets ] )
print 'unmapped: %s' % '\t'.join( [ '%.1f' % ( 100. * x / total_unmapped ) for x in unmapped_buckets ] )

# entropy
mapped_buckets = bio.bucket( stats.mapped['entropy'], buckets )
unmapped_buckets = bio.bucket( stats.unmapped['entropy'], buckets )
total_mapped = sum( mapped_buckets )
total_unmapped = sum( unmapped_buckets )

print '\n========== Entropy =========='
print 'Mapped:   min: %.2f max: %.2f' % ( min( stats.mapped['entropy'] ), max( stats.mapped['entropy'] ) )
print 'Unmapped: min: %.2f max: %.2f' % ( min( stats.unmapped['entropy'] ), max( stats.unmapped['entropy'] ) )
print 'Entropy:  %s' % '\t'.join( [ '%.2f' % bucket for bucket in buckets ] )
print 'mapped:   %s' % '\t'.join( [ '%.1f' % ( 100. * x / total_mapped ) for x in mapped_buckets ] )
print 'unmapped: %s' % '\t'.join( [ '%.1f' % ( 100. * x / total_unmapped ) for x in unmapped_buckets ] )

