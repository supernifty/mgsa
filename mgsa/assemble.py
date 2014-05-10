#
# build a samfile with bowtie
#
import os
import sys

def bowtie_assemble( fasta_file, fastq_file, sam_file ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  os.system( '../../tools/bowtie2-2.1.0/bowtie2-build %s %s-bt2' % ( fasta_file, fasta_file ) )
  os.system( '../../tools/bowtie2-2.1.0/bowtie2 --local -p 16 -x %s-bt2 -U %s -t -S %s' % ( fasta_file, fastq_file, sam_file ) )

if len(sys.argv) != 4:
  print "Usage: %s reference.fasta reads.fastq output.sam" % sys.argv[0]
  sys.exit(0)

bowtie_assemble( sys.argv[1], sys.argv[2], sys.argv[3] )
