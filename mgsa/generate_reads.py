
import sys
import bio

if len(sys.argv) < 2:
  print "Usage: %s config_file [vcf variations]" % sys.argv[0]
  sys.exit(1)

config_helper = bio.Config()
cfg = config_helper.read_config_file( open( sys.argv[1], 'r' ) )

if len(sys.argv) >= 3:
  vcf = bio.VCF( reader=open(sys.argv[2], 'r') )
else:
  vcf = None

if len(sys.argv) == 4:
  variation_map = open(sys.argv[3], 'w')
else:
  variation_map = None

bio.FastqGenerator( sys.stdin, cfg, vcf, variation_map )

