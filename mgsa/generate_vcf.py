
import sys

import bio

if len(sys.argv) != 2:
  print "Usage: %s reference_fasta < sam_file > vcf_file" % sys.argv[0]
  sys.exit(0)

# python generate_fasta .reference < .sam > .vcf
converter = bio.SamToVCF( sam=sys.stdin, reference=open( sys.argv[1], 'r' ), target_vcf=bio.VCF( writer=bio.VCFWriter( sys.stdout ) ), log=bio.log_stderr )
