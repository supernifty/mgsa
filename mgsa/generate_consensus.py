
import sys

import bio

# python generate_fasta < .sam > .fasta
converter = bio.SamToFasta( sys.stdin, bio.log_stderr )
converter.write( sys.stdout )
