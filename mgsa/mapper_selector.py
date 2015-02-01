#
# build a samfile with bowtie or other aligner
#
import os
import sys

import bio
import config

if __name__ == "__main__":
  if len(sys.argv) != 5:
    print "Usage: %s [baseline|bowtie2|bwa|soap] reference.fasta reads.fastq output.sam" % sys.argv[0]
    sys.exit(0)

  bio.get_mapper( config, sys.argv[1], sys.argv[2] )( sys.argv[3], sys.argv[4] )
