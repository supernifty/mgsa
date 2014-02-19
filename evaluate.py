# evaluate alignment software
import datetime
import os
import sys

import bio

REFERENCE_PREFIX = 'test1k'
CANDIDATE_PREFIX = 'test1ksnp'

# assume snp is there TODO later
def generate_mutant():
  m = bio.FastaMutate( bio.FastaReader( open( '../data/test/%s.fasta' % REFERENCE_PREFIX, 'r' ) ), log, 0.1 )
  o = open( '../data/test/%s.fasta' % CANDIDATE_PREFIX, 'w' ) 
  o.write( '>%s\n' % CANDIDATE_PREFIX )
  for line in m.items():
    o.write( '%s\n' % line )
  o.close()
  log( "%i mutations" % m.mutations )

def generate_reads():
  # generate reads (.sam)
  os.system( 'java -jar -Xmx2048m ../tools/SimSeq/SimSeqNBProject/store/SimSeq.jar -1 100 -2 100 --error ../tools/SimSeq/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt --error2 ../tools/SimSeq/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt -l 300 -s 30 -n 5000 -r ../data/test/%s.fasta -o ../data/test/%s.sam -u 0.01' % ( CANDIDATE_PREFIX, CANDIDATE_PREFIX ) )

  # generate reads (.fastq)
  os.system( 'java -jar ../tools/picard-tools-1.96/SamToFastq.jar I=../data/test/%s.sam FASTQ=../data/test/%s_R1.fastq F2=../data/test/%s_R2.fastq VALIDATION_STRINGENCY=SILENT' % ( CANDIDATE_PREFIX, CANDIDATE_PREFIX, CANDIDATE_PREFIX ) )

def bowtie_assemble():
  # evaluate bowtie - build bowtie sam
  os.system( '../tools/bowtie2-2.1.0/bowtie2-build ../data/test/%s.fasta ../data/test/%s-bt2' % ( REFERENCE_PREFIX, REFERENCE_PREFIX ) )
  os.system( '../tools/bowtie2-2.1.0/bowtie2 --local -p 16 -x ../data/test/%s-bt2 -1 ../data/test/%s_R1.fastq -2 ../data/test/%s_R2.fastq -S ../data/test/%s_bowtie.sam -t' % ( REFERENCE_PREFIX, CANDIDATE_PREFIX, CANDIDATE_PREFIX, CANDIDATE_PREFIX ) )

def log(msg):
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  sys.stderr.write( '%s: %s\n' % ( when, msg ) )

def bowtie_evaluate():
  # sam to fasta
  candidate_fasta = bio.SamToFasta( open( '../data/test/%s_bowtie.sam' % CANDIDATE_PREFIX, 'r' ), log )
  log( 'size is %i' % candidate_fasta.fasta.size )

  #log( s.fasta.prior )
  #for i in xrange(1, 10):
  #  print i, s.fasta.count(i), s.fasta.confidence(i)

  #pos = 1
  #while pos < s.fasta.size:
  #  print s.fasta.consensus( pos, pos + 70 )
  #  pos += 70

  candidate_diff = bio.FastaDiff( bio.FastaReader( open( '../data/test/%s.fasta' % CANDIDATE_PREFIX, 'r' ) ), candidate_fasta.fasta, log )
  log( "%i errors to candidate" % candidate_diff.error_total )
  log( candidate_diff.errors )

if __name__ == '__main__':
  #generate_mutant()  
  #generate_reads()  
  #bowtie_assemble()
  bowtie_evaluate()

