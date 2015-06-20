
import argparse
import collections
import sys

import bio
import config

parser = argparse.ArgumentParser(description='Compare BAMs')
parser.add_argument('bams', metavar='bams', nargs='+', help='bam files to analyze')
parser.add_argument('--mapq', dest='mapq_min', type=int, default=-1, help='only consider alignments with at least this mapq value')
parser.add_argument('--compare_position', dest='compare_position', type=bool, default=False, help='compare position of alignments')
parser.add_argument('--subset_detail', dest='subset_detail', type=bool, default=False, help='more detail about subsets')
parser.add_argument('--mismatch_detail', dest='mismatch_detail', type=int, default=-1, help='more detail about mismatches for this bam index')

args = parser.parse_args()

diff = bio.SamDiff( [ bio.BamReaderExternal( config.BAM_TO_SAM, sam_file ) for sam_file in args.bams ], mapq_min=args.mapq_min, compare_position=args.compare_position, subset_detail=args.subset_detail, mismatch_detail=None if args.mismatch_detail == -1 else args.mismatch_detail )

print "mapq stats\n=========="
print "i: n\tmax\tmin\tmean\tsd\tfilename"
for idx, stats in enumerate( diff.mapq_stats ):
  print '%i: %i\t%.2f\t%.2f\t%.2f\t%.2f\t%s' % ( idx, stats['mapped'], stats['max'], stats['min'], stats['mean'], stats['sd'], args.bams[idx]  )

print "\nmapped vs unmapped commonality\n==================="
for key in sorted( diff.totals.keys() ):
  print ( "{0:0%ib}: {1}" % ( len(args.bams) ) ).format( key, diff.totals[key] )

if args.compare_position:
  print "\nmapped vs unmapped commonality including position differences\n==================="
  for key in sorted( diff.position_totals.keys() ):
    print ( "{0:0%ib}: {1}" % ( len(args.bams) ) ).format( key, diff.position_totals[key] )

if args.subset_detail:
  print "\nmapq vs position differences\n==================="
  print "i: max\tmin\tmean\tsd\thist"
  for key, value in diff.mapq_subset_stats.items():
    bin_key = ( "{0:0%ib}" % ( len(args.bams) ) ).format( key )
    print '%s: %.2f\t%.2f\t%.2f\t%.2f\t%s' % ( bin_key, value['max'], value['min'], value['mean'], value['sd'], value['hist'] )

if args.mismatch_detail > -1:
  print "\nmismatch details\n==================="
  print "pos,alt"
  mismatch_count = 0
  unpaired_count = 0
  for read, value in diff.mismatch_stats.items():
    if 'pos' in value and 'alt' in value:
      print '%i,%i' % ( value['pos'], value['alt'] )
      mismatch_count += 1
    else:
      unpaired_count += 1
      #print "missing values", read, value
  print mismatch_count, "mismatches with incorrect alternatives;", unpaired_count, "unpaired reads"
