# uses mult and mult_snp_prob from config

import sys
import bio

if len(sys.argv) != 3:
  print "Usage: %s config fasta" % sys.argv[0]
  sys.exit(0)

cfg = bio.Config().read_config_file( open( sys.argv[1], 'r' ) )
multiplier = int(cfg['mult'])

multi_reader = bio.MultiFastaReader( open( sys.argv[2], 'r' ) )
bio.RepeatedMultiFastaGenerator( multi_reader, sys.stdout, multiplier, cfg )
