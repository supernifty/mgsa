
#####################################
# Pipeline for evaluating an aligner
# 
# requires a config file on the command line
# config file is of the form (all params optional):
# mult 1
# snp_prob 0.01
# insert_prob 0
# delete_prob 0
# max_insertion_len 1
# max_deletion_len 1
#####################################

import bio
import datetime
import os
import sys

def run( cmd ):
  bio.log_stderr( cmd )
  os.system( cmd )

if len(sys.argv) < 2:
  #print "Usage: %s multiplier (snp|ins|del|none) (bowtie2|bowtie2_e2e|bwa|bwasw)" % sys.argv[0]
  print "Usage: %s config_file" % sys.argv[0]
  sys.exit(0)

# read configuration
cfg = { 'insert_prob': 0, 'delete_prob': 0, 'snp_prob': 0, 'max_insertion_len': 1, 'max_deletion_len': 1, 'mult': 1, 'mapper': 'bowtie2' }
floats = ( 'insert_prob', 'delete_prob', 'snp_prob' )
ints = ( 'max_insertion_len', 'max_deletion_len', )

mutation_type = ''
for line in open( sys.argv[1], 'r' ):
  if line.startswith( '#' ):
    continue
  key, value = line.strip().split()
  if key in floats:
    cfg[key] = float(value)
  elif key in ints:
    cfg[key] = int(value)
  else:
    cfg[key] = value
  # generate brief description
  if key in ( 'insert_prob', 'delete_prob', 'snp_prob' ):
    mutation_type = mutation_type + '_' + key[:3]
  if key in ( 'max_insertion_len', 'max_deletion_len' ):
    mutation_type = mutation_type + '_' + value

mutation_type = mutation_type[1:]

# write configuration for generate_mutation
when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
cfg_file = '../../data/pipeline_mutation_%s.cfg' % when
with open( cfg_file, 'w' ) as cfg_out:
  for key in cfg:
    cfg_out.write( '%s %s\n' % ( key, str(cfg[key]) ) )

# generates a reference sequence that consists of repeated copies of the provided sequence
#run( "python build_repeated_fasta.py ../../data/circoviridae.fasta %s > ../../data/circoviridae_x%s.fasta" % ( cfg['mult'], cfg['mult'] ) )

# generates a donor sequence that contains mutations from the reference sequence
run( "python generate_mutation.py ../../data/circoviridae_x%s.fasta ../../data/circoviridae_%s_x%s.vcf %s > ../../data/circoviridae_%s_x%s.fasta" % ( cfg['mult'], mutation_type, cfg['mult'], cfg_file, mutation_type, cfg['mult'] ) )

# generates reads from the donor sequence
run( "python generate_reads.py %s ../../data/circoviridae_%s_x%s.vcf < ../../data/circoviridae_%s_x%s.fasta > ../../data/circoviridae_%s_x%s.fastq" % ( cfg_file, mutation_type, cfg['mult'], mutation_type, cfg['mult'], mutation_type, cfg['mult'] ) )

# runs an aligner
run( "python mapper_selector.py %s ../../data/circoviridae_x%s.fasta ../../data/circoviridae_%s_x%s.fastq ../../data/circoviridae_%s_x%s_%s.sam" % ( cfg['mapper'], cfg['mult'], mutation_type, cfg['mult'], mutation_type, cfg['mult'], cfg['mapper'] ) )

# evaluates the alignment accuracy
run( "python evaluate_reads.py %s < ../../data/circoviridae_%s_x%s_%s.sam > ../../data/circoviridae_%s_x%s_%s_%s_evaluation.txt" % ( cfg_file, mutation_type, cfg['mult'], cfg['mapper'], mutation_type, cfg['mult'], cfg['mapper'], when ) )

