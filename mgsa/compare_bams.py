
import collections
import sys

import bio
import config

if len(sys.argv) < 2:
  print "Usage: %s bam1 bam2 ..." % sys.argv[0]
  sys.exit(1)

diff = bio.SamDiff( [ bio.BamReaderExternal( config.BAM_TO_SAM, sam_file ) for sam_file in sys.argv[1:] ] )

for key in sorted( diff.totals.keys() ):
  print ( "{0:0%ib}: {1}" % ( len(sys.argv) - 1 ) ).format( key, diff.totals[key] )
