
# converts bams to features for classification
#
# example usage:
#   extract_features.py bam1 bam2 ... bamAll

import argparse
import sys

import config
import bio

parser = argparse.ArgumentParser(description="Extract features from bam files")
parser.add_argument( 'bams', metavar='bams', nargs='+', help='bam files to extract reads' )
opts = parser.parse_args()

classes = []
for bam in opts.bams[:-1]:
  sys.stderr.write( "processing %s...\n" % bam )
  tags = bio.SamTags( bio.BamReaderExternal( config.BAM_TO_SAM, bam ) )
  classes.append( tags.tags )
  sys.stderr.write( "processing %s: found %i tags\n" % ( bam, len(tags.tags) ) )

sys.stderr.write( "building feature file from %s...\n" % opts.bams[-1] )

bio.SamFeatures( bio.BamReaderExternal( config.BAM_TO_SAM, opts.bams[-1] ), sys.stdout, classes, exclude_unmapped=True, log=bio.log_stderr )
