
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
# see bio.Config for full list...
#####################################

import bio
import config
import datetime
import os
import sys

def run( cmd ):
  '''
    run a system command
  '''
  bio.log_stderr( cmd )
  os.system( cmd )

if len(sys.argv) < 3:
  print "Usage: %s config_file output_file" % sys.argv[0]
  sys.exit(0)

# open config file
# if a field is added, update helpers.py
target = open( sys.argv[2], 'w' )
target.write( '# cfg, unmapped, incorrect, read_precision, read_recall, read_f1, vcf_tp, vcf_fp, vcf_fn, vcf_precision, vcf_recall, vcf_f1, vcf_bucket_tp, vcf_bucket_fp, vcf_bucket_fn' )
target.write( ', reference_bias, error_bias, unmapped_variations, total_variations, mean_reference, mean_error, error_bias_no_variation' ) # part of BiasReport
target.write( '\n' )
first = True
config_helper = bio.Config()
config_file = sys.argv[1]

# evaluate each line and create config object
for line in open( config_file, 'r' ):
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
  if cfg['fasta_reference'] == 'generated':
    fasta_file = "../../data/%s_x%s.fasta" % (cfg['fasta'], cfg['mult'] )

    if cfg['command'] == 'fasta':
      run( "python build_repeated_fasta.py %s ../../data/%s.fasta > %s" % ( cfg_file, cfg['fasta'], fasta_file ) )
      continue # next line
  else:
    fasta_file = "../../data/%s" % cfg['fasta_reference']

  bio.log_stderr( "fasta reference: %s" % fasta_file )

  if cfg['vcf_source'] == 'generated':
    vcf_file = "../../data/%s_%s_x%s.vcf" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'] )
  else:
    vcf_file = "../../data/%s" % cfg['vcf_source']

  # generates a donor sequence that contains mutations from the reference sequence and writes the mutations to vcf_file
  if cfg['command'] == 'mutations':  
    run( "python generate_mutation.py %s %s %s > ../../data/%s_%s_x%s.fasta" % ( fasta_file, vcf_file, cfg_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'] ) )
    continue

  # generates reads from the donor sequence
  variation_map_file = '../../data/%s_%s_x%s.vmap' % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'] )
  if cfg['command'] == 'reads':  
    run( "python generate_reads.py %s ../../data/%s_%s_x%s.vcf %s < ../../data/%s_%s_x%s.fasta > ../../data/%s_%s_x%s.fastq" % ( cfg_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'], variation_map_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['fasta'], cfg['mutation_type'], cfg['mult'] ) )
    continue

  # if sam not provided, run an aligner and evaluate read alignment
  # default read alignment results
  unmapped = 0
  incorrect = 0
  recall = 0

  if cfg['sam_source'] == 'generated':
    sam_file =  "../../data/%s_%s_x%s_%s.sam" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'] )
    run( "python mapper_selector.py %s %s ../../data/%s_%s_x%s.fastq %s" % ( cfg['mapper'], fasta_file, cfg['fasta'], cfg['mutation_type'], cfg['mult'], sam_file ) )

    # evaluates the alignment accuracy
    eval_file = "../../data/%s_%s_x%s_%s_%s_evaluation.txt" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
    run( "python evaluate_reads.py %s %s < %s > %s" % ( cfg_file, variation_map_file, sam_file, eval_file ) )
    # parse results
    for eval_line in open( eval_file, 'r' ):
      if eval_line.startswith( '%mapped_unmapped' ):
        unmapped = float( eval_line.strip().split()[1] )
      if eval_line.startswith( '%mapped_incorrectly' ):
        incorrect = float( eval_line.strip().split()[1] )
      if eval_line.startswith( '%mapped_recall' ):
        recall = float( eval_line.strip().split()[1] ) / 100.

  else:
    sam_file =  "../../data/%s" % cfg['sam_source'] # existing sam file

  # calculate extra stats
  precision = ( 100 - incorrect ) / 100.
  f1 = 0
  if recall != 0 or precision != 0:
    f1 = 2 * ( precision * recall ) / ( precision + recall )
      
  # generate consensus dna
  #consensus_file = "../../data/%s_%s_x%s_%s_%s_recovered.fasta" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
  #run( "python generate_consensus.py < %s > %s" % ( sam_file, consensus_file ) )
  # generate vcf
  if cfg['vcf_out'] == 'generated':
    recovered_vcf_file = "../../data/%s_%s_x%s_%s_%s_recovered.vcf" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
  else:
    recovered_vcf_file = "../../data/%s" % ( cfg['vcf_out'] )
  #run( "python generate_vcf.py %s < %s > %s" % ( fasta_file, sam_file, recovered_vcf_file ) )

  # build candidate vcf
  candidate_vcf = bio.VCF( writer=bio.VCFWriter( open( recovered_vcf_file, 'w' ) ) ) # empty vcf

  if cfg['command'] is None or cfg['command'].find( '-novcf' ) == -1:
    if sam_file.endswith('.bam'):
      sam_fh = bio.BamReaderExternal( config.BAM_TO_SAM, sam_file )
    else:
      sam_fh = open( sam_file, 'r' )
    if cfg['chromosomes'] == 'true':
      bio.SamToMultiChromosomeVCF( sam=sam_fh, multi_fasta_reference=open( fasta_file, 'r' ), target_vcf=candidate_vcf, log=bio.log_stderr, call_strategy=cfg['call_strategy'] ) # write to vcf
    else:
      bio.SamToVCF.instance( sam=sam_fh, reference=open( fasta_file, 'r' ), target_vcf=candidate_vcf, log=bio.log_stderr, call_strategy=cfg['call_strategy'] ) # write to vcf

  if cfg['command'] == 'vcf':
    # don't compare correct vcf to candidate vcf, we just wanted to generate a vcf
    target.write( '%s' % ( line.strip() ) )

  if cfg['command'] is not None and cfg['command'].startswith( 'features' ):
    fasta = bio.SamToFasta( sam=open( sam_file, 'r' ), log=bio.log_stderr ).fasta
    feature = bio.MapFeature( fasta ) # what the aligner has built
    if cfg['command'].find( '-depth' ) != -1:
      depths = [ str( feature.depth( x ) ) for x in xrange(0, fasta.length) ]
      target.write( '%s\n' % ','.join( depths ) )
    if cfg['command'].find( '-breakpoints' ) != -1:
      breakpoints = [ '%.2f' % bio.BREAKPOINT_PREFIX_READ_DEPTH_SHAPE.evidence( fasta, x, mean_depth=cfg['coverage'] ) for x in xrange(0, fasta.length) ]
      target.write( '%s\n' % ','.join( breakpoints ) )
      breakpoints = [ '%.2f' % bio.BREAKPOINT_SUFFIX_READ_DEPTH_SHAPE.evidence( fasta, x, mean_depth=cfg['coverage'] ) for x in xrange(0, fasta.length) ]
      target.write( '%s\n' % ','.join( breakpoints ) )
    if cfg['command'].find( '-deletions' ) != -1:
      detector = bio.LongDeletionDetector( 50, vcf_writer=candidate_vcf, threshold=0.6, reference_fasta=bio.Fasta( bio.FastaReader( open( fasta_file, 'r' ) ) ) ) # TODO size hardcoded
      deletions = [ '%.2f' % detector.evidence( fasta, x, mean_depth=cfg['coverage'] ) for x in xrange(0, fasta.length) ]
      target.write( '%s\n' % ','.join( deletions ) )
    # also write vcf sites
    if cfg['command'].find( '-novcf' ) == -1:
      vcf = bio.VCF( reader=open( vcf_file, 'r' ) )
      target.write( ','.join( [ str(x) for x in vcf.breakpoints( all_affected=True ) ] ) )

  if cfg['command'] is None or cfg['command'].startswith( 'vcfdiff' ): # else: # if no command, compare correct vcf to candidate vcf
    if cfg['command'].find( '-deletions' ) != -1:
      fasta = bio.SamToFasta( sam=open( sam_file, 'r' ), log=bio.log_stderr ).fasta
      detector = bio.LongDeletionDetector( 50, vcf_writer=candidate_vcf, threshold=0.5, reference_fasta=bio.Fasta( bio.FastaReader( open( fasta_file, 'r' ) ) ), log=bio.log_stderr ) # TODO size hardcoded
      [ detector.evidence( fasta, x, mean_depth=cfg['coverage'] ) for x in xrange(0, fasta.length) ] # writes to vcf_writer
    vcf_diff = bio.VCFDiff( vcf_correct=bio.VCF(reader=open( vcf_file, 'r' ), log=bio.log_stderr), vcf_candidate=candidate_vcf, log=bio.log_stderr )
    #print vcf_diff.stats
    vcf_precision = 0
    vcf_recall = 0
    vcf_f1 = 0
    if vcf_diff.stats['tp'] > 0 or vcf_diff.stats['fp'] > 0:
      vcf_precision = 1. * vcf_diff.stats['tp'] / ( vcf_diff.stats['tp'] + vcf_diff.stats['fp'] )
    elif vcf_diff.stats['fp'] == 0: 
      vcf_precision = 1 # nothing wrong
    if vcf_diff.stats['tp'] > 0 or vcf_diff.stats['fn'] > 0:
      vcf_recall = 1. * vcf_diff.stats['tp'] / ( vcf_diff.stats['tp'] + vcf_diff.stats['fn'] )
    if vcf_precision != 0 or vcf_recall != 0:
      vcf_f1 = 2. * ( vcf_precision * vcf_recall ) / ( vcf_precision + vcf_recall ) 
    target.write( '%s,%f,%f,%f,%f,%f,%i,%i,%i,%f,%f,%f' % ( line.strip(), unmapped, incorrect, precision * 100, recall * 100, f1 * 100, vcf_diff.stats['tp'], vcf_diff.stats['fp'], vcf_diff.stats['fn'], vcf_precision * 100, vcf_recall * 100, vcf_f1 * 100 ) )
    target.write( ',%s,%s,%s' % ( '|'.join( [ str( x['tp'] ) for x in vcf_diff.buckets ] ), '|'.join( [ str( x['fp'] ) for x in vcf_diff.buckets ] ),'|'.join( [ str( x['fn'] ) for x in vcf_diff.buckets ] ) ) )

  # bias report - may as well always have this on (?)
  if cfg['reports'] is not None and cfg['reports'].find('bias') != -1:
    report = bio.BiasReport( 
      candidate_fasta=bio.SamToFasta( sam=open( sam_file, 'r' ), log=bio.log_stderr ).fasta, # what the aligner has built
      reference_fasta=open( fasta_file, 'r' ), # the reference
      donor_vcf=bio.VCF(reader=open( vcf_file, 'r' ), log=bio.log_stderr), # variations to donor
      log=bio.log_stderr, 
      buckets_count=cfg['bias_report_buckets'] )
    #bias_report_file = "out/%s_%s_x%s_%s_%s_bias.txt" % ( cfg['fasta'], cfg['mutation_type'], cfg['mult'], cfg['mapper'], when )
    #print "stats", report.stats, "reference", report.reference_histogram, "error", report.error_histogram
    #target.write( ',%s,%s,%.2f,%.2f,%f,%f,%s' % ( '|'.join( [str(x) for x in report.reference_histogram ] ), '|'.join( [ str(x) for x in report.error_histogram ] ), report.stats['unmapped'], report.stats['total'], report.stats['mean_reference'], report.stats['mean_error'], '|'.join( [str(x) for x in report.error_histogram_no_variation ] ) ) )
    target.write( ',%s,%s,%.2f,%.2f,%f,%f' % ( '|'.join( [str(x) for x in report.reference_histogram ] ), '|'.join( [ str(x) for x in report.error_histogram ] ), report.stats['unmapped'], report.stats['total'], report.stats['mean_reference'], report.stats['mean_error'] ) )

  target.write( '\n' )
  target.flush()
