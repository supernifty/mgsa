
# filter bams with a list of tagnames
# writes out a filtered sam

import argparse
import sys

import bio
import config

parser = argparse.ArgumentParser( description='Filter BAM on tag' )
parser.add_argument('bam', help='input bam file')
parser.add_argument('--tags', dest='tags', required=True, help='file containing tags')
parser.add_argument('--exclude', dest='exclude', type=bool, required=False, default=False, help='exclude tags instead of include')
args = parser.parse_args()

bio.log_stderr( 'Building list of tags...' )
allowed_tags = set()
for line in open( args.tags, 'r' ):
  allowed_tags.add( line.strip() )
bio.log_stderr( '%i allowed tags' % len(allowed_tags) )

bio.SamFilter( sam_fh=bio.BamReaderExternal( config.BAM_TO_SAM, args.bam ), target_fh=sys.stdout, allowed_tags=allowed_tags, exclude=args.exclude, log=bio.log_stderr )

