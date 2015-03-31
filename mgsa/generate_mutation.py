#
# given a fasta file, mutate it into another fasta file
#
import bio
import sys

if len(sys.argv) < 3:
  print "Usage: %s reference.fasta output.vcf mutation_profile.cfg" % sys.argv[0]
  sys.exit(0)

# read configuration
if len(sys.argv) > 3:
  cfg = bio.Config().read_config_file( open( sys.argv[3] ) )
else:
  bio.Config.default_config

mutator = bio.FastaMutate( bio.FastaReader( open( sys.argv[1], 'r' ) ), vcf_file=open( sys.argv[2], 'w' ), snp_prob=cfg['snp_prob'], insert_prob=cfg['insert_prob'], delete_prob=cfg['delete_prob'], min_insert_len=cfg['min_insertion_len'], max_insert_len=cfg['max_insertion_len'], min_delete_len=cfg['min_deletion_len'], max_delete_len=cfg['max_deletion_len'], min_variation_start=cfg['min_variation_start'], min_variation_dist=cfg['min_variation_distance'], insert_source=cfg['insertion_source'] )
sys.stdout.write( '>mutated %s\n' % sys.argv[1] )
for fragment in mutator.items():
  sys.stdout.write( '%s\n' % fragment )

