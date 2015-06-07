
# given an xmfa, remap a bam file coordinates
# usage remap_bam --xmfa xmfa --origin 0 --target 1 srcbamfile > targetsamfile
# writes coverage stats to stderr

import argparse
import config
import sys

import bio

parser = argparse.ArgumentParser(description='Remap a bam file using XMFA')
parser.add_argument('bam', help='bam file to analyze')
parser.add_argument('--xmfa', required=True, help='xmfa file')
parser.add_argument('--origin', type=int, default=1, help='index of fasta that bam was mapped to')
parser.add_argument('--target', type=int, default=2, help='index of fasta that bam should be mapped to')
parser.add_argument('--output', required=True, help='sam output file')

args = parser.parse_args()

# measure coverage
print "building coverage map..."
#print args.origin, dest
mauve_map = bio.MauveMap( open(args.xmfa, 'r') )
print len(mauve_map.coverage), "positions mapped"

print "analyzing bam..."
if args.bam.endswith( '.bam' ):
  sam_fh = bio.BamReaderExternal( config.BAM_TO_SAM, args.bam )
else:
  sam_fh = open( args.bam, 'r' )

mauve_map.remap( sam_fh, open( args.output, 'w' ) )
 
print "====== Coverage Stats ====="
for key in mauve_map.stats.keys():
  print "%20s: %i" % ( key, mauve_map.stats[key] )