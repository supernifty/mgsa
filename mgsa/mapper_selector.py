#
# build a samfile with bowtie or other aligner
#
import os
import sys

def baseline_map( fasta_file, fastq_file, sam_file, params=None ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  os.system( 'python mapper.py %s %s %s' % ( fasta_file, fastq_file, sam_file ) )

def bowtie_map( fasta_file, fastq_file, sam_file, params ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  os.system( '../../tools/bowtie2-2.1.0/bowtie2-build %s %s-bt2' % ( fasta_file, fasta_file ) )
  os.system( '../../tools/bowtie2-2.1.0/bowtie2 %s -p 16 -x %s-bt2 -U %s -t -S %s' % ( params, fasta_file, fastq_file, sam_file ) )

def bwa_map( fasta_file, fastq_file, sam_file, params ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  os.system( '../../tools/bwa-0.7.5a/bwa index -p %s-bwa %s' % ( fasta_file, fasta_file ) )
  os.system( '../../tools/bwa-0.7.5a/bwa %s %s-bwa %s > %s' % ( params, fasta_file, fastq_file, sam_file ) )

if len(sys.argv) != 5:
  print "Usage: %s [bowtie2|bwa] reference.fasta reads.fastq output.sam" % sys.argv[0]
  sys.exit(0)

if sys.argv[1] == 'baseline':
  baseline_map( sys.argv[2], sys.argv[3], sys.argv[4] )
elif sys.argv[1] == 'bowtie2':
  bowtie_map( sys.argv[2], sys.argv[3], sys.argv[4], '--local' )
elif sys.argv[1] == 'bowtie2_e2e':
  bowtie_map( sys.argv[2], sys.argv[3], sys.argv[4], '--end-to-end' )
elif sys.argv[1] == 'bwa':
  bwa_map( sys.argv[2], sys.argv[3], sys.argv[4], 'mem' )
elif sys.argv[1] == 'bwasw':
  bwa_map( sys.argv[2], sys.argv[3], sys.argv[4], 'bwasw' )
else:
  print "Unrecognized mapper: %s" % sys.argv[1]
