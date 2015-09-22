
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
parser.add_argument('--new_reference', required=False, help='new reference name for sam output')
parser.add_argument('--remap_cigar', required=False, default=False, help='attempt to remap the cigar string')
parser.add_argument('--output_not_covered', required=False, help='sam output file for aligned reads with no map')
parser.add_argument('--output_target_coverage', required=False, help='bed file of covered bases')

args = parser.parse_args()

# measure coverage
print "building coverage map..."
#print args.origin, dest
mauve_map = bio.MauveMap( open(args.xmfa, 'r'), src_strand=args.origin, target_strand=args.target, new_reference=args.new_reference )
print len(mauve_map.coverage), "positions mapped"

if args.output_target_coverage:
  print "generating target bed..."
  target = set()
  for k, v in mauve_map.coverage.iteritems():
    target.add( v )
  print "writing target bed..."
  with open( args.output_target_coverage, 'w' ) as ofh:
    range_start = None
    for x in xrange(mauve_map.genome_stats['ymin'], mauve_map.genome_stats['ymax']+2): # +2 to make sure range ends
      if x in target:
        if range_start is None:
          range_start = x
      else:
        if range_start is not None:
          ofh.write( '%s\t%i\t%i\n' % ( args.new_reference, range_start, x ) )
          range_start = None

print "analyzing bam..."
if args.bam.endswith( '.bam' ):
  sam_fh = bio.BamReaderExternal( config.BAM_TO_SAM, args.bam )
else:
  sam_fh = open( args.bam, 'r' )

if args.output_not_covered:
    mauve_map.remap( sam_fh, open( args.output, 'w' ), remap_cigar=args.remap_cigar, not_covered_output=open( args.output_not_covered, 'w' ) )
else:
    mauve_map.remap( sam_fh, open( args.output, 'w' ), remap_cigar=args.remap_cigar )
 
sys.stdout.write( "====== Mapping Stats =====\n" )
for key in mauve_map.genome_stats.keys():
  sys.stdout.write( "%20s: %i\n" % ( key, mauve_map.genome_stats[key] ) )

sys.stdout.write( "====== Coverage Stats =====\n" )
for key in mauve_map.stats.keys():
  sys.stdout.write( "%20s: %i\n" % ( key, mauve_map.stats[key] ) )
