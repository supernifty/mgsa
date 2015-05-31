
# given an xmfa, remap a bam file coordinates
# usage remap_bam --xmfa xmfa --origin 0 bam 
# currently just calculates coverage stats

import argparse
import config

import bio

parser = argparse.ArgumentParser(description='Measure coverage')
parser.add_argument('bam', help='bam file to analyze')
parser.add_argument('--xmfa', required=True, help='xmfa file')
parser.add_argument('--origin', type=int, default=0, help='index of fasta that bam was mapped to')

args = parser.parse_args()

# measure coverage
print "building coverage map..."
dest = ( (args.origin + 1) % 2 ) * 2
#print args.origin, dest
covered = set()
for pos, line in enumerate(open(args.xmfa, 'r')):
  if pos == 0:
    continue
  fields = line.strip().split()
  start = int( fields[args.origin * 2])
  stop = int( fields[args.origin * 2 + 1] )
  alt_start = int( fields[dest] )
  alt_stop = int( fields[dest + 1] )
  if stop > start and alt_start != 0:
    for pos in xrange(start, stop+1):
      covered.add(pos)
print len(covered), "positions covered"

print "analyzing bam..."
if args.bam.endswith( '.bam' ):
  sam_fh = bio.BamReaderExternal( config.BAM_TO_SAM, args.bam )
else:
  sam_fh = open( args.bam, 'r' )
stats = { 'total': 0, 'unmapped': 0, 'mapped': 0, 'mapped_covered': 0, 'mapped_uncovered': 0, 'mapped_partial': 0 }
for pos, line in enumerate(sam_fh):
  fields = line.split()
  if len(fields) > 9:
    #print fields
    flag = int( fields[1] )
    stats['total'] += 1
    if flag & 0x04 != 0:
      stats['unmapped'] += 1
    else:
      stats['mapped'] += 1
      # assess coverage
      start_pos = int(fields[3])
      read_len = len(fields[9]) # TODO not always correct
      if all( [ x in covered for x in xrange( start_pos, start_pos + read_len ) ] ):
        stats['mapped_covered'] += 1
      elif any( [ x in covered for x in xrange( start_pos, start_pos + read_len ) ] ):
        stats['mapped_partial'] += 1
      else:
        stats['mapped_uncovered'] += 1
  if pos < 10 or pos % 100000 == 0:
    print "processed %i lines: %s..." % ( pos, stats )

print "====== Coverage Stats ====="
for key in stats.keys():
  print "%20s: %i" % ( key, stats[key] )
