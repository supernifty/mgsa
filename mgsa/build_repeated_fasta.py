# uses mult and mult_snp_prob from config

import sys
import bio

if len(sys.argv) != 3:
  print "Usage: %s config fasta" % sys.argv[0]
  sys.exit(0)

cfg = bio.Config().read_config_file( open( sys.argv[1], 'r' ) )

multiplier = int(cfg['mult'])

sys.stdout.write( '>generated fasta %ix\n' % multiplier )

# write out unmodified sequence
#sequence = ''
reader = bio.FastaReader( open( sys.argv[2], 'r' ) )
for item in reader.items():
  if len(item) > 0:
    sys.stdout.write( '%s\n' % item ) # write unmodified
    #sequence += '%s\n' % item # remember

for i in xrange(0, multiplier - 1):
  reader = bio.FastaReader( open( sys.argv[2], 'r' ) )
  mutator = bio.FastaMutate( reader, snp_prob=float(cfg['mult_snp_prob']), insert_prob=0, delete_prob=0 )
  for item in mutator.items():
    if len(item) > 0:
      sys.stdout.write( '%s\n' % item ) # write unmodified
   #sys.stdout.write( sequence )

