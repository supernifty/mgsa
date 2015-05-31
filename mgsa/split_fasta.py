
import argparse
import sys

current = 0
fh = None

parser = argparse.ArgumentParser()
parser.add_argument( '--prefix', type=str, default='split', required=False, help='File prefix')
args = parser.parse_args()

for line in sys.stdin:
  if line.startswith( '>' ):
    if fh is not None:
      fh.close()
    current += 1
    fn = '%s-%02i.fasta' % (args.prefix, current)
    fh = open( fn, 'w' )
    sys.stderr.write( 'Writing to %s...\n' % fn )
  fh.write( line )
