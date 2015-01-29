#
# build a samfile with bowtie or other aligner
#
import os
import sys

import config

def baseline_map( fasta_file, fastq_file, sam_file, params=None ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  def result():
    os.system( 'python mapper.py %s %s %s' % ( fasta_file, fastq_file, sam_file ) )
  return result

def bowtie_map( fasta_file, params ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  def result( fastq_file, sam_file, index=True ):
    if index:
      os.system( '%s-build %s %s-bt2' % ( config.BOWTIE_PATH, fasta_file, fasta_file ) )
    os.system( '%s %s -p 16 -x %s-bt2 -U %s --quiet -S %s' % ( config.BOWTIE_PATH, params, fasta_file, fastq_file, sam_file ) ) # -t adds time
  return result

def bwa_map( fasta_file, params ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  def result( fastq_file, sam_file, index=True ):
    if index:
      os.system( '%s index -p %s-bwa %s' % ( config.BWA_PATH, fasta_file, fasta_file ) )
    os.system( '%s %s %s-bwa -v 2 %s > %s' % ( config.BWA_PATH, params, fasta_file, fastq_file, sam_file ) )
  return result

def soap_map( fasta_file ):
  def result( fastq_file, sam_file, index=True ):
    if index:
      os.system( '%s/2bwt-builder %s' % ( config.SOAP_PATH, fasta_file ) )
    os.system( '%s/soap -a %s -D %s.index -o %s' % ( config.SOAP_PATH, fastq_file, fasta_file, sam_file ) )
  return result

if __name__ == "__main__":
  if len(sys.argv) != 5:
    print "Usage: %s [baseline|bowtie2|bwa|soap] reference.fasta reads.fastq output.sam" % sys.argv[0]
    sys.exit(0)
  
  if sys.argv[1] == 'baseline':
    baseline_map( sys.argv[2], sys.argv[3], sys.argv[4] )()
  elif sys.argv[1] == 'bowtie2':
    bowtie_map( sys.argv[2], '--local' )(sys.argv[3], sys.argv[4])
  elif sys.argv[1] == 'bowtie2_e2e':
    bowtie_map( sys.argv[2], '--end-to-end' )(sys.argv[3], sys.argv[4])
  elif sys.argv[1] == 'bwa':
    bwa_map( sys.argv[2], 'mem' )(sys.argv[3], sys.argv[4])
  elif sys.argv[1] == 'bwasw':
    bwa_map( sys.argv[2], 'bwasw' )(sys.argv[3], sys.argv[4] )
  elif sys.argv[1] == 'soap':
    soap_map( sys.argv[2] )(sys.argv[3], sys.argv[4] )
  else:
    print "Unrecognized mapper: %s" % sys.argv[1]
