
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

if len(sys.argv) < 3:
  print "Usage: %s config_file output_file" % sys.argv[0]
  sys.exit(0)

target = open( sys.argv[2], 'w' )
target.write( '# cfg, unmapped, incorrect, read_precision, read_recall, read_f1, vcf_tp, vcf_fp, vcf_fn, vcf_precision, vcf_recall, vcf_f1, vcf_bucket_tp, vcf_bucket_fp, vcf_bucket_fn' )
target.write( ', reference_bias, error_bias' )
target.write( '\n' )
first = True
config_helper = bio.Config()
for line in open( sys.argv[1], 'r' ):
  if line.startswith( '#' ) or line.strip() == '':
    target.write( line )
    continue
  cfg = config_helper.read_config_line( line )
  cfg['mutation_type'] = cfg['mutation_type'][1:]

  # write configuration for generate_mutation
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  # for generate_mutation
  cfg_file = '../../data/pipeline_mutation_%s.cfg' % when
  with open( cfg_file, 'w' ) as cfg_out:
    for key in cfg:
      cfg_out.write( '%s %s\n' % ( key, str(cfg[key]) ) )


  # generates a reference sequence that consists of repeated copies of the provided sequence
  if cfg['command'] == 'fasta':
    run( "python build_repeated_fasta.py ../../data/%s.fasta %s > ../../data/%s_x%s.fasta" % ( cfg['fasta'], cfg['mult'], cfg['fasta'], cfg['mult'] ) )
    continue # next line

  fasta_file = "../../data/%s_x%s.fasta" % (cfg['fasta'], cfg['mult'] )
  vcf_file = "../../data/%s_%s_x%s.vcf" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'] )

  # generates a donor sequence that contains mutations from the reference sequence
  if cfg['command'] == 'mutations':  
    run( "python generate_mutation.py %s %s %s > ../../data/%s_%s_x%s.fasta" % ( fasta_file, vcf_file, cfg_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'] ) )
    continue

  # generates reads from the donor sequence
  variation_map_file = '../../data/%s_%s_x%s.vmap' % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'] )
  if cfg['command'] == 'reads':  
    run( "python generate_reads.py %s ../../data/%s_%s_x%s.vcf %s < ../../data/%s_%s_x%s.fasta > ../../data/%s_%s_x%s.fastq" % ( cfg_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'], variation_map_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['fasta'], cfg['mutation_type'], cfg['mult'] ) )
    continue

  # runs an aligner
  sam_file =  "../../data/%s_%s_x%s_%s.sam" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'] )
  run( "python mapper_selector.py %s %s ../../data/%s_%s_x%s.fastq %s" % ( cfg['mapper'], fasta_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'], sam_file ) )

  # evaluates the alignment accuracy
  eval_file = "../../data/%s_%s_x%s_%s_%s_evaluation.txt" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
  run( "python evaluate_reads.py %s %s < %s > %s" % ( cfg_file, variation_map_file, sam_file, eval_file ) )

  # parse results
  unmapped = 0
  incorrect = 0
  recall = 0
  for eval_line in open( eval_file, 'r' ):
    if eval_line.startswith( '%mapped_unmapped' ):
      unmapped = float( eval_line.strip().split()[1] )
    if eval_line.startswith( '%mapped_incorrectly' ):
      incorrect = float( eval_line.strip().split()[1] )
    if eval_line.startswith( '%mapped_recall' ):
      recall = float( eval_line.strip().split()[1] ) / 100.

  precision = ( 100 - incorrect ) / 100.

  f1 = 0
  if recall != 0 or precision != 0:
    f1 = 2 * ( precision * recall ) / ( precision + recall )
      
  # generate consensus dna
  #consensus_file = "../../data/%s_%s_x%s_%s_%s_recovered.fasta" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
  #run( "python generate_consensus.py < %s > %s" % ( sam_file, consensus_file ) )
  # generate vcf
  recovered_vcf_file = "../../data/%s_%s_x%s_%s_%s_recovered.vcf" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
  #run( "python generate_vcf.py %s < %s > %s" % ( fasta_file, sam_file, recovered_vcf_file ) )

  # build candidate vcf
  candidate_vcf = bio.VCF( writer=bio.VCFWriter( open( recovered_vcf_file, 'w' ) ) )
  bio.SamToVCF( sam=open( sam_file, 'r' ), reference=open( fasta_file, 'r' ), target_vcf=candidate_vcf, log=bio.log_stderr )
  # compare correct vcf to candidate vcf
  vcf_diff = bio.VCFDiff( vcf_correct=bio.VCF(reader=open( vcf_file, 'r' ), log=bio.log_stderr), vcf_candidate=candidate_vcf, log=bio.log_stderr )
  #print vcf_diff.stats
  vcf_precision = 0
  vcf_recall = 0
  vcf_f1 = 0
  if vcf_diff.stats['tp'] > 0 or vcf_diff.stats['fp'] > 0:
    vcf_precision = 1. * vcf_diff.stats['tp'] / ( vcf_diff.stats['tp'] + vcf_diff.stats['fp'] )
  if vcf_diff.stats['tp'] > 0 or vcf_diff.stats['fn'] > 0:
    vcf_recall = 1. * vcf_diff.stats['tp'] / ( vcf_diff.stats['tp'] + vcf_diff.stats['fn'] )
  if vcf_precision != 0 or vcf_recall != 0:
    vcf_f1 = 2. * ( vcf_precision * vcf_recall ) / ( vcf_precision + vcf_recall ) 
  target.write( '%s,%f,%f,%f,%f,%f,%i,%i,%i,%f,%f,%f' % ( line.strip(), unmapped, incorrect, precision * 100, recall * 100, f1 * 100, vcf_diff.stats['tp'], vcf_diff.stats['fp'], vcf_diff.stats['fn'], vcf_precision * 100, vcf_recall * 100, vcf_f1 * 100 ) )
  target.write( ',%s,%s,%s' % ( '|'.join( [ str( x['tp'] ) for x in vcf_diff.buckets ] ), '|'.join( [ str( x['fp'] ) for x in vcf_diff.buckets ] ),'|'.join( [ str( x['fn'] ) for x in vcf_diff.buckets ] ) ) )

  # bias report
  if cfg['reports'] is not None and cfg['reports'].find('bias') != -1:
    report = bio.BiasReport( 
      candidate_fasta=bio.SamToFasta( sam=open( sam_file, 'r' ), log=bio.log_stderr ).fasta, # what the aligner has built
      reference_fasta=open( fasta_file, 'r' ), # the reference
      donor_vcf=bio.VCF(reader=open( vcf_file, 'r' ), log=bio.log_stderr), # variations to donor
      log=bio.log_stderr, 
      buckets=cfg['bias_report_buckets'] )
    #bias_report_file = "out/%s_%s_x%s_%s_%s_bias.txt" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
    #print "stats", report.stats, "reference", report.reference_histogram, "error", report.error_histogram
    target.write( ',%s,%s' % ( '|'.join( [str(x) for x in report.reference_histogram ] ), '|'.join( [ str(x) for x in report.error_histogram ] ) ) )

  target.write( '\n' )
  target.flush()
