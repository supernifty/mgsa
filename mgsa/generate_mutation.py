#
# given a fasta file, mutate it into another fasta file
#
import bio
import sys

if len(sys.argv) < 3:
  print "Usage: %s reference.fasta output.vcf [snp_prob] [insert_prob] [delete_prob]" % sys.argv[0]
  sys.exit(0)

snp_prob = 0.01
insert_prob = 0
delete_prob = 0

if len(sys.argv) > 3:
  snp_prob = float(sys.argv[3])
if len(sys.argv) > 4:
  insert_prob = float(sys.argv[4])
if len(sys.argv) > 5:
  delete_prob = float(sys.argv[5])

mutator = bio.FastaMutate( bio.FastaReader( open( sys.argv[1], 'r' ) ), vcf_file=open( sys.argv[2], 'w' ), snp_prob=snp_prob, insert_prob=insert_prob, delete_prob=delete_prob )
sys.stdout.write( '>mutated %s\n' % sys.argv[1] )
for fragment in mutator.items():
  sys.stdout.write( '%s\n' % fragment )

