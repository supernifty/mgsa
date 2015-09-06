
# using output of compare_bams mismatch, find nearest target from mauve

# first read the xmfa

import argparse
import sys

import bio

parser = argparse.ArgumentParser(description='Amend mismatch information with nearest')
parser.add_argument('--xmfa', required=True, help='xmfa file')
parser.add_argument('--origin', type=int, default=1, help='index of fasta that bam was mapped to')
parser.add_argument('--target', type=int, default=2, help='index of fasta that bam should be mapped to')

args = parser.parse_args()

# measure coverage
bio.log_stderr( "building coverage map..." )
#print args.origin, dest
mauve_map = bio.MauveMap( open(args.xmfa, 'r'), src_strand=args.origin, target_strand=args.target )
bio.log_stderr( "%i positions mapped" % len(mauve_map.coverage) )

bio.log_stderr( "processing CSV..." )
for idx, line in enumerate(sys.stdin):
  if idx == 0:
    continue # skip header
  fields = line.strip().split(',')
  if len(fields) > 2:
    alt = int(fields[1])
    nearest = mauve_map.find_nearest_target( alt )
    sys.stdout.write( '%s,%i\n' % ( line.strip(), nearest ) )
    if idx % 100 == 0:
      bio.log_stderr( "processed %i lines..." % idx )
  else:
    bio.log_stderr( "skipped line %i" % idx )
bio.log_stderr( "processing CSV: done" )
