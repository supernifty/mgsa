
import sys
import bio

if len(sys.argv) != 3:
  print "Usage: %s fasta multiple" % sys.argv[0]
  sys.exit(0)

fasta = bio.FastaReader( open( sys.argv[1] ) )
multiplier = int( sys.argv[2] )

sys.stdout.write( '>generated fasta %ix\n' % multiplier )
sequence = ''
for item in fasta.items():
  if len(item) > 0:
    sys.stdout.write( '%s\n' % item )
    sequence += '%s\n' % item

for i in xrange(0, multiplier - 1):
  sys.stdout.write( sequence )

