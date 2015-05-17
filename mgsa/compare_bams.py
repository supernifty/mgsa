
import collections
import sys

import bio
import config

if len(sys.argv) < 2:
  print "Usage: %s bam1 bam2 ..." % sys.argv[0]
  sys.exit(1)

diff = bio.SamDiff( [ bio.BamReaderExternal( config.BAM_TO_SAM, sam_file ) for sam_file in sys.argv[1:] ] )

print "mapq stats\n=========="
for idx, stats in enumerate( diff.mapq_stats ):
  print '%i: %s\tcount: %i\tmax: %.2f\tmin: %.2f\tmean: %.2f\tsd: %.2f' % ( idx, sys.argv[idx + 1], stats['count'], stats['max'], stats['min'], stats['mean'], stats['sd'] )

print "\nmapping commonality\n==================="
for key in sorted( diff.totals.keys() ):
  print ( "{0:0%ib}: {1}" % ( len(sys.argv) - 1 ) ).format( key, diff.totals[key] )
