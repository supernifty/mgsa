# evaluate alignment software
import ConfigParser as configparser
import datetime
import os
import sys

import bio

import matplotlib.pyplot as plt

cfg = configparser.ConfigParser()

#config('reference_prefix') = 'test1k'
#config('reference_prefix') = 'test10k'
#config('reference_prefix') = 'test100k'

#config('candidate_prefix') = 'test1ksnp'
#config('candidate_prefix') = 'test10ksnp'
#config('candidate_prefix') = 'test100ksnp'
#config('candidate_prefix') = 'test100kindel'
#config('candidate_prefix') = 'test100kins'
#config('candidate_prefix') = 'test100kdel'
#config('candidate_prefix') = 'test1kins'
#config('candidate_prefix') = 'test100del'

#config('snp_prob') = 0.00
#config('insert_prob') = 0.01
#config('delete_prob') = 0.01

#config_int('read_count') = 5000
#config_int('read_length') = 100
#config_int('read_library_length') = 300
# coverage = 5000 * 100 * 2 / length
# - 5000 * 100 * 2 / 100k = 10

ALLOW_SAM_INDELS = True

def config( name ):
  return cfg.get( 'evaluate', name )

def config_float( name ):
  return cfg.getfloat( 'evaluate', name )

def config_int( name ):
  return cfg.getint( 'evaluate', name )

def config_bool( name ):
  return cfg.getboolean( 'evaluate', name )

# assume snp is there TODO later
def generate_mutant():
  '''generates a mutated fasta file (candidate)'''
  m = bio.FastaMutate( bio.FastaReader( open( '../data/test/%s.fasta' % config('reference_prefix'), 'r' ) ), log, open( '../data/test/%s.vcf' % config('candidate_prefix'), 'w' ), config_float('snp_prob'), config_float('insert_prob'), config_float('delete_prob') )
  o = open( '../data/test/%s.fasta' % config('candidate_prefix'), 'w' ) 
  o.write( '>%s\n' % config('candidate_prefix') )
  for line in m.items():
    o.write( '%s\n' % line )
  o.close()
  log( "%i mutations" % m.mutations )

def generate_reads():
  '''candidate.fasta to .sam to candidate.fastq'''
  # generate reads (.sam)
  # -1 -2 read lengths
  # -n reads
  # -l library insert size mean
  # -s library insert size stdev
  # -u duplicate probability
  os.system( 'java -jar -Xmx2048m ../tools/SimSeq/SimSeqNBProject/store/SimSeq.jar -1 %i -2 %i --error ../tools/SimSeq/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt --error2 ../tools/SimSeq/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt -l %i -s 30 -n %i -r ../data/test/%s.fasta -o ../data/test/%s.sam -u 0.01' % ( config_int('read_length'), config_int('read_length'), config_int('read_library_length'), config_int('read_count'), config('candidate_prefix'), config('candidate_prefix') ) )

  # generate reads (.fastq)
  os.system( 'java -jar ../tools/picard-tools-1.96/SamToFastq.jar I=../data/test/%s.sam FASTQ=../data/test/%s_R1.fastq F2=../data/test/%s_R2.fastq VALIDATION_STRINGENCY=SILENT' % ( config('candidate_prefix'), config('candidate_prefix'), config('candidate_prefix') ) )

def bowtie_assemble():
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  os.system( '../tools/bowtie2-2.1.0/bowtie2-build ../data/test/%s.fasta ../data/test/%s-bt2' % ( config('reference_prefix'), config('reference_prefix') ) )
  os.system( '../tools/bowtie2-2.1.0/bowtie2 --local -p 16 -x ../data/test/%s-bt2 -1 ../data/test/%s_R1.fastq -2 ../data/test/%s_R2.fastq -S ../data/test/%s_bowtie.sam -t' % ( config('reference_prefix'), config('candidate_prefix'), config('candidate_prefix'), config('candidate_prefix') ) )

def log(msg):
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  sys.stderr.write( '%s: %s\n' % ( when, msg ) )
  logfile.write( '%s: %s\n' % ( when, msg ) )

def bowtie_evaluate():
  # sam to fasta
  candidate_fasta = bio.SamToFasta( open( '../data/test/%s_bowtie.sam' % config('candidate_prefix'), 'r' ), log, allow_indels=ALLOW_SAM_INDELS )
  log( 'size is %i' % candidate_fasta.fasta.size )

  # write fasta
  o = open( '../data/test/%s_bowtie.fasta' % config('candidate_prefix'), 'w' )
  o.write( '>reassembled-bowtie-consensus\n')
  pos = 1
  while pos < candidate_fasta.fasta.size:
    o.write( '%s\n' % candidate_fasta.fasta.consensus( pos, pos + 70 ) )
    pos += 70
  o.close()

  log( 'overall counts %s' % candidate_fasta.fasta.prior )
  total = sum( [ candidate_fasta.fasta.prior[k] for k in candidate_fasta.fasta.prior ] )
  log( 'overall %% %s' % [ '%s: %f' % ( k, 100.0 * candidate_fasta.fasta.prior[k] / total ) for k in candidate_fasta.fasta.prior ] )
  #for i in xrange(1, 10):
  #  print i, s.fasta.count(i), s.fasta.confidence(i)

  #pos = 1
  #while pos < s.fasta.size:
  #  print s.fasta.consensus( pos, pos + 70 )
  #  pos += 70

  # compare candidate fasta to what we just built
  candidate_diff = bio.FastaDiff( bio.FastaReader( open( '../data/test/%s.fasta' % config('candidate_prefix'), 'r' ) ), candidate_fasta.fasta, log )
  log( "%i insertions" % len(candidate_fasta.fasta.inserted) )
  log( "%i deletions" % len(candidate_fasta.fasta.deleted) )
  log( "%i differences to candidate" % candidate_diff.error_total )
  log( candidate_diff.errors )
  #log( candidate_diff.error_positions )

  # compare reference fasta to what we just built
  #reference_diff = bio.FastaDiff( bio.FastaReader( open( '../data/test/%s.fasta' % config('reference_prefix'), 'r' ) ), candidate_fasta.fasta, log )
  #log( "%i insertions" % len(candidate_fasta.fasta.inserted) )
  #log( "%i deletions" % len(candidate_fasta.fasta.deleted) )
  #log( "%i differences to reference" % reference_diff.error_total )
  #log( reference_diff.errors )

  # see what vcfs were found
  vcf_diff = bio.VCFDiff( bio.VCF(reader=open('../data/test/%s.vcf' % config('candidate_prefix'), 'r')), candidate_fasta.fasta, log )
  #log( candidate_fasta.fasta.deleted )
  log( 'snps: %s' % vcf_diff.snp_stats )
  log( 'inserts: %s' % vcf_diff.ins_stats )
  log( 'deletes: %s' % vcf_diff.del_stats )
  
  # where are the errors?
  n, bins, patches = plt.hist(list(candidate_diff.error_positions), bins=100)
  plt.show()

if __name__ == '__main__':
  if len(sys.argv) != 2:
    print "Usage: %s config_file" % sys.argv[0]
    exit(0)

  cfg.readfp(open(sys.argv[1]))
  logfile = open( '../data/test/%s-%s.log' % ( config('candidate_prefix'), datetime.datetime.now().strftime("%Y%m%d-%H%M%S") ), 'w' )
  if config_bool( 'generate_mutant' ):
    generate_mutant()  
  if config_bool( 'generate_reads' ):
    generate_reads()  
  if config_bool( 'bowtie_assemble' ):
    bowtie_assemble()
  if config_bool( 'bowtie_evaluate' ):
    bowtie_evaluate()

