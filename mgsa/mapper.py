
# a simple baseline aligner with simple configurable behaviour
# given reference fasta and fastq reads, generates a sam alignment file

import bio
import sys
import re

def find_all_exact_matches( sub, fasta ):
  return [m.start() for m in re.finditer('(?=%s)' % sub, fasta)]

if len(sys.argv) != 4:
  print "Usage: %s reference.fasta reads.fastq output.sam" % sys.argv[0]
  sys.exit(0)

# read the fasta into memory
bio.log_stderr( 'reading fasta' )
fasta = ''

for line in open( sys.argv[1], 'r' ):
  if line.startswith( '>' ):
    continue
  fasta += line.strip()

# align reads
fastq = open( sys.argv[2], 'r' )
sam = bio.SamWriter( open( sys.argv[3], 'w' ) )

# sam header
sam.write_header( fasta )

while True:
  name = fastq.readline()
  if not name:
    break
  name = name.strip()
  dna = fastq.readline().strip()
  extra = fastq.readline().strip()
  confidence = fastq.readline().strip()
  # map it
  matches = find_all_exact_matches( dna, fasta )
  #bio.log_stderr( 'matches for %s: %s' % ( dna, matches ) )
  if len(matches) == 0: # no exact match
    sam.write_no_match( name, dna, extra, confidence )
  elif len(matches) == 1: # one exact match
    sam.write_match( name, dna, extra, confidence, matches[0] )
  else:
    sam.write_multi_match( name, dna, extra, confidence, matches )

