
import argparse
import collections
import sys

import bio
import config

parser = argparse.ArgumentParser(description='Compare BAMs')
parser.add_argument('bams', metavar='bams', nargs='+', help='bam files to analyze')
parser.add_argument('--mapq', dest='mapq_min', type=int, default=-1, help='only consider alignments with at least this mapq value')
parser.add_argument('--compare_position', dest='compare_position', type=bool, default=False, help='compare position of alignments')

args = parser.parse_args()

diff = bio.SamDiff( [ bio.BamReaderExternal( config.BAM_TO_SAM, sam_file ) for sam_file in args.bams ], mapq_min=args.mapq_min, compare_position=args.compare_position )

print "mapq stats\n=========="
for idx, stats in enumerate( diff.mapq_stats ):
  print '%i: count: %i\tmax: %.2f\tmin: %.2f\tmean: %.2f\tsd: %.2f\t%s' % ( idx, stats['mapped'], stats['max'], stats['min'], stats['mean'], stats['sd'], args.bams[idx]  )

print "\nmapped vs unmapped commonality\n==================="
for key in sorted( diff.totals.keys() ):
  print ( "{0:0%ib}: {1}" % ( len(args.bams) ) ).format( key, diff.totals[key] )

if args.compare_position:
  print "\nmapped vs unmapped commonality including position differences\n==================="
  for key in sorted( diff.position_totals.keys() ):
    print ( "{0:0%ib}: {1}" % ( len(args.bams) ) ).format( key, diff.position_totals[key] )
