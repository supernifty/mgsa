
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
parser.add_argument('--xmfa', dest='xmfa', help='include mapping nearness in mismatch detail')
parser.add_argument('--origin', type=int, default=1, help='index of fasta that bam was mapped to')
parser.add_argument('--target', type=int, default=2, help='index of fasta that bam should be mapped to')

args = parser.parse_args()

diff = bio.SamDiff( [ bio.BamReaderExternal( config.BAM_TO_SAM, sam_file ) for sam_file in args.bams ], mapq_min=args.mapq_min, compare_position=args.compare_position, subset_detail=args.subset_detail, mismatch_detail=None if args.mismatch_detail == -1 else args.mismatch_detail )

sys.stdout.write( "mapq stats\n==========\n" )
sys.stdout.write( "i:\tn\tmax\tmin\tmean\tsd\tfilename\n" )
for idx, stats in enumerate( diff.mapq_stats ):
  sys.stdout.write( '%i:\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n' % ( idx, stats['mapped'], stats['max'], stats['min'], stats['mean'], stats['sd'], args.bams[idx] ) )

sys.stdout.write( "\nmapped vs unmapped commonality\n===================\n" )
for key in sorted( diff.totals.keys() ):
  sys.stdout.write( ( "{0:0%ib}: {1}\n" % ( len(args.bams) ) ).format( key, diff.totals[key] ) )

if args.compare_position:
  sys.stdout.write( "\nmapped vs unmapped commonality including position differences\n===================\n" )
  for key in sorted( diff.position_totals.keys() ):
    sys.stdout.write( ( "{0:0%ib}: {1}\n" % ( len(args.bams) ) ).format( key, diff.position_totals[key] ) )

if args.subset_detail:
  sys.stdout.write( "\nmapq vs position differences\n===================\n" )
  sys.stdout.write( "i:\tmax\tmin\tmean\tsd\thist\n" )
  for key, value in diff.mapq_subset_stats.items():
    bin_key = ( "{0:0%ib}" % ( len(args.bams) ) ).format( key )
    sys.stdout.write( '%s:\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n' % ( bin_key, value['max'], value['min'], value['mean'], value['sd'], value['hist'] ) )

if args.mismatch_detail > -1:
  sys.stdout.write( "\nmismatch details\n===================\n" )
  if args.xmfa is None:
    sys.stdout.write( "pos,alt,read_id,wrongness\n" )
  else:
    mauve_map = bio.MauveMap( open(args.xmfa, 'r'), src_strand=args.origin, target_strand=args.target )
    sys.stdout.write( "pos,alt,read_id,wrongness,nearestmap\n" )
  mismatch_count = 0
  unpaired_count = 0
  for read, value in diff.mismatch_stats.items():
    if 'p' in value and 'a' in value:
      if args.xmfa is None:
        sys.stdout.write( '%i,%i,%s,%i\n' % ( value['p'], value['a'], read, value['p'] - value['a'] ) )
      else:
        nearest = mauve_map.find_nearest_target( int(value['a']) )
        sys.stdout.write( '%i,%i,%s,%i,%i\n' % ( value['p'], value['a'], read, value['p'] - value['a'], nearest ) )
      mismatch_count += 1
    else:
      unpaired_count += 1
      #print "missing values", read, value
  bio.log_stderr( "%i mismatches with incorrect alternatives; %i unpaired reads" % ( mismatch_count, unpaired_count ) )

