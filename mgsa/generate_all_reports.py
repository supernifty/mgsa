#
# build the charts used in the final report
# generally these assume the analysis has been done; the command to do this is included in the function
#

import datetime
import math
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys

import bio
import helpers

REPORT_DIRECTORY = '../../report/plots'

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
pylab.rcParams['figure.figsize'] = 6, 4

##########################################################################
# helpers
def run( cmd ):
  bio.log_stderr( cmd )
  os.system( cmd )

def is_numeric( field ):
  return re.match( '^[0-9\\-.]+$', field ) is not None

columns = [ 'unmapped', 'incorrect', 'read_precision', 'read_recall', 'read_f1', 'vcf_tp', 'vcf_fp', 'vcf_fn', 'vcf_precision', 'vcf_recall', 'vcf_f1', 'vcf_bucket_tp', 'vcf_bucket_fp', 'vcf_bucket_fn', 'reference_bias', 'error_bias', 'unmapped_variations', 'total_variations', 'mean_reference', 'mean_error' ]
def find_column( line, name ):
  fields = line.split(',')
  # find first numeric
  for idx in xrange(0, len(fields)):
    if is_numeric( fields[idx] ):
      break
  # find referenced column
  for col in xrange(0, len(columns)):
    if columns[col] == name:
      return fields[idx + col]
  return None

def find_parameter( line, name ):
  '''
    find a parameter e.g. read_length 50
  '''
  fields = line.split( ',' )
  for field in fields:
    keyvalue = field.strip().split(' ')
    if keyvalue[0] == name and len(keyvalue) == 2:
      return keyvalue[1]
  return None
  
def get_value( column, fh, default=None ):
  while True:
    line = fh.readline()
    if line == "":
      return None
    if line.startswith( '#' ):
      continue
    value = find_column( line, column )
    if value is not None:
      if default is not None and default[0] == float(value):
        return default[1]
      else:
        return value

def get_values( count, column, fh ):
  result = []
  while len(result) < count:
    value = get_value( column, fh )
    result.append( float(value ) )
  return result

def get_multicolumn_value( columns, fh ):
  while True: # find first valid line
    line = fh.readline()
    if line == "":
      return None
    if line.startswith( '#' ):
      continue
    result = []
    for column in columns: # find each column
      columns.append( find_column( line, column ) )
    return result

##########################################################################
# run job
def mutation_hiv():
  out_file = "out/pipeline_batch_mutation_hiv_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_mutation_hiv.cfg %s" % out_file )
  
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_base = []
  y_bow = []
  y_bwa = []
  for x in xrange(0, 11):
    y_base.append( get_value( 'unmapped', fh ) )
    y_bow.append( get_value( 'unmapped', fh ) )
    y_bwa.append( get_value( 'unmapped', fh ) )
  
  # draw graph
  x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y_base, label='Baseline', color='b')
  ax.plot(x, y_bow, label='Bowtie2', color='g')
  ax.plot(x, y_bwa, label='BWA', color='r')
  ax.set_ylabel('Unmapped %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/mutation-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
  
def mutation_hiv_snp_f1():
  out_file = "out/pipeline_batch_mutation_hiv_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_mutation_hiv.cfg %s" % out_file )
  
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_base = []
  y_bow = []
  y_bwa = []
  for x in xrange(0, 11):
    y_base.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    y_bow.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    y_bwa.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    
  # draw graph
  x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y_base, label='Baseline', color='b')
  ax.plot(x, y_bow, label='Bowtie2', color='g')
  ax.plot(x, y_bwa, label='BWA', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/mutation-f1-snp-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def mutation_ecoli_snp_f1():
  out_file = "out/pipeline_batch_mutation_ecoli_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_mutation_ecoli_140812.cfg %s" % out_file )
  
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_bow = []
  y_bwa = []
  for x in xrange(0, 11):
    y_bow.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    y_bwa.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    
  # draw graph
  x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y_bow, label='Bowtie2', color='g')
  ax.plot(x, y_bwa, label='BWA', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutations-snps-2.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
 
def mutation_ecoli_snp_histogram():
  out_file = "out/pipeline_batch_mutation_ecoli_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  bio.log_stderr( 'extracting values from %s...' % out_file )

  fh = open( out_file, 'r' )
  y_bow_tp = []
  y_bwa_tp = []
  for x in xrange(0, 11):
    y_bow_tp.append( get_value( 'vcf_bucket_tp', fh) )
    y_bwa_tp.append( get_value( 'vcf_bucket_tp', fh) )
    
  fh = open( out_file, 'r' )
  y_bow_fn = []
  y_bwa_fn = []
  for x in xrange(0, 11):
    y_bow_fn.append( get_value( 'vcf_bucket_fn', fh) )
    y_bwa_fn.append( get_value( 'vcf_bucket_fn', fh) )

  fh = open( out_file, 'r' )
  y_bow_fp = []
  y_bwa_fp = []
  for x in xrange(0, 11):
    y_bow_fp.append( get_value( 'vcf_bucket_fp', fh) )
    y_bwa_fp.append( get_value( 'vcf_bucket_fp', fh) )
    
  # bowtie
  # just 0.4 i.e. [5]
  #tp, fp, fn = [ x.split('|') for x in bowtie_0_4.split(',') ]
  GENOME = 4641652. 
  tp = [ int(x) / GENOME * 20 * 100 for x in y_bow_tp[5].split('|') ]
  fp = [ int(x) / GENOME * 20 * 100 for x in y_bow_fp[5].split('|') ]
  fn = [ int(x) / GENOME * 20 * 100 for x in y_bow_fn[5].split('|') ]
  x = np.linspace( 0, 95, 20 )
  #x
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.bar(x, fp, width=4, color='r', label='False Positives')
  ax.bar(x, fn, width=4, color='b', label='False Negatives', bottom=fp)
  ax.set_ylabel('Errors %')
  ax.set_xlabel('Genome Position %')
  ax.set_ylim(ymax=0.7)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutations-bowtie-histogram.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

  # bwa
  tp = [ int(x) / GENOME * 20 * 100 for x in y_bwa_tp[5].split('|') ]
  fp = [ int(x) / GENOME * 20 * 100 for x in y_bwa_fp[5].split('|') ]
  fn = [ int(x) / GENOME * 20 * 100 for x in y_bwa_fn[5].split('|') ]
  x = np.linspace( 0, 95, 20 )
  #x
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.bar(x, fp, width=4, color='r', label='False Positives')
  ax.bar(x, fn, width=4, color='b', label='False Negatives', bottom=fp)
  ax.set_ylabel('Errors %')
  ax.set_xlabel('Genome Position %')
  ax.set_ylim(ymax=0.7)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutations-bwa-histogram.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

  bio.log_stderr( 'extracting values from %s: done' % out_file )

def coverage_ecoli():
  out_file = "out/pipeline_batch_coverage_ecoli_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_coverage_ecoli.cfg %s" % out_file )
  
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  b5 = []
  b10 = []
  b20 = []
  for x in xrange(0, 11):
    b5.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    b10.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    b20.append( get_value( 'vcf_f1', fh, default=(0,100) ) )

  # draw graph
  x = np.linspace( 1, 10, 10 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, b5[1:], label='5x coverage', color='b')
  ax.plot(x, b10[1:], label='10x coverage', color='g')
  ax.plot(x, b20[1:], label='20x coverage', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-coverage-snp.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
 
def read_length_ecoli():
  out_file = "out/pipeline_batch_read_length_ecoli_%s_inchworm.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_read_length_ecoli.cfg %s" % out_file )

  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  b100 = []
  b200 = []
  b400 = []
  for x in xrange(0, 10):
    b100.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    b200.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    b400.append( get_value( 'vcf_f1', fh, default=(0,100) ) )

  #print b100
  # draw graph
  x = np.linspace( 1, 10, 10 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, b100, label='100bp', color='b')
  ax.plot(x, b200, label='200bp', color='g')
  ax.plot(x, b400, label='400bp', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  ax.set_ylim(ymin=75)
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-read-length-bowtie-snp.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
  
def read_length_ecoli_map():
  out_file = "out/pipeline_batch_read_length_ecoli_%s_inchworm.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_read_length_ecoli.cfg %s" % out_file )

  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  b100 = []
  b200 = []
  b400 = []
  for x in xrange(0, 10):
    b100.append( get_value( 'read_f1', fh, default=(0,100) ) )
    b200.append( get_value( 'read_f1', fh, default=(0,100) ) )
    b400.append( get_value( 'read_f1', fh, default=(0,100) ) )

  #print b100
  # draw graph
  x = np.linspace( 1, 10, 10 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, b100, label='100bp', color='b')
  ax.plot(x, b200, label='200bp', color='g')
  ax.plot(x, b400, label='400bp', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-read-length-bowtie-mapping.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
 
def read_length_hiv():
  out_file = "out/pipeline_batch_mutation_hiv_read_length_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_mutation_hiv_read_length.cfg %s" % out_file )

  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  b100 = []
  t100 = []
  t200 = []
  t400 = []
  for x in xrange(0, 11):
    b100.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    t100.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    t200.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    t400.append( get_value( 'vcf_f1', fh, default=(0,100) ) )

  print b100
  print t100
  print t200
  print t400
  # draw graph
  x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, b100, label='Baseline 100bp', color='#ff9933')
  ax.plot(x, t100, label='Bowtie2 100bp', color='b')
  ax.plot(x, t200, label='Bowtie2 200bp', color='g')
  ax.plot(x, t400, label='Bowtie2 400bp', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/hiv-read-length-vcf-snps.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def snp_vs_map_ecoli():
  out_file = "out/pipeline_batch_read_length_ecoli_%s_inchworm.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_read_length_ecoli.cfg %s" % out_file )

  bio.log_stderr( 'extracting values from %s...' % out_file )

  fh = open( out_file, 'r' )
  s100 = []
  s200 = []
  s400 = []
  for x in xrange(0, 10):
    s100.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    s200.append( get_value( 'vcf_f1', fh, default=(0,100) ) )
    s400.append( get_value( 'vcf_f1', fh, default=(0,100) ) )

  fh = open( out_file, 'r' )
  m100 = []
  m200 = []
  m400 = []
  for x in xrange(0, 10):
    m100.append( get_value( 'read_f1', fh, default=(0,100) ) )
    m200.append( get_value( 'read_f1', fh, default=(0,100) ) )
    m400.append( get_value( 'read_f1', fh, default=(0,100) ) )

  #print 400bp
  # draw graph
  x = np.linspace( 1, 10, 10 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, m400, label='Align', color='b')
  ax.plot(x, s400, label='SNP', color='g')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-map_vs_snp-400.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

  #print 100bp
  # draw graph
  x = np.linspace( 1, 10, 10 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, m100, label='Align', color='b')
  ax.plot(x, s100, label='SNP', color='g')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  #ax.set_ylim(ymin=75)
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-map_vs_snp-100.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
 
# experimental
def read_length_ecoli_map_low_mutation():
  #out_file = "out/pipeline_batch_read_length_ecoli_low_mutation_%s_inchworm.out" % datetime.datetime.now().strftime("%Y%m%d")
  out_file = "out/pipeline_batch_read_length_ecoli_low_mutation_20141030_inchworm.out"
  #run( "python pipeline_batch.py batch/pipeline_batch_read_length_ecoli_low_mutation.cfg %s" % out_file )

  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  b50 = []
  b100 = []
  b200 = []
  b400 = []
  for x in xrange(0, 9):
    b50.append( get_value( 'read_f1', fh, default=(0,100) ) )
    b100.append( get_value( 'read_f1', fh, default=(0,100) ) )
    b200.append( get_value( 'read_f1', fh, default=(0,100) ) )
    b400.append( get_value( 'read_f1', fh, default=(0,100) ) )

  #print b100
  # draw graph
  x = np.linspace( 0.1, 0.9, 9 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, b50, label='50bp', color='m')
  ax.plot(x, b100, label='100bp', color='b')
  ax.plot(x, b200, label='200bp', color='g')
  ax.plot(x, b400, label='400bp', color='r')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-read-length-bowtie-mapping-low-mutation.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
 
def plot_reference_bias_ecoli_example(include_zero=False):
  out_file = "out/pipeline_batch_bias_ecoli_6.out"
  #run( "python pipeline_batch.py batch/pipeline_batch_bias_ecoli_2.cfg %s" % out_file )

  xs, ys = helpers.series_from_pipeline_result( fh=open(out_file), y='reference_bias', bias_report=True )
  xs = [ x * 100 for x in xs ]

  if include_zero:
    plt.xlim( xmin=-5, xmax=110 )
  else:
    plt.xlim( xmin=5, xmax=110 )

  #print xs, ys

  # weird crap
  ys = [ y * 100 + 1e-7 for y in ys ]
  width = 100. / len(ys) * 0.75

  #plt.bar( left=xs, height=ys, label='Bias', width=width )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlim(xmax=105)
  if include_zero:
    ax.bar( left=[x for x in xs], height=ys, label='SNPs affected (%)', width=width, color='b', log=False, align='center')
  else:
    ax.bar( left=[x for x in xs][1:], height=ys[1:], label='SNPs affected (%)', width=width, color='b', log=False, align='center')
    ax.set_xlim(xmin=5)
  ax.set_ylabel('SNPs (%)')
  ax.set_xlabel('Reference Bias (%)')
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/reference-bias-profile.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def plot_reference_bias_ecoli_example_with_errors(include_zero=False):
  out_file = "out/bias_ecoli_141127.out"
  #run( "python pipeline_batch.py batch/pipeline_batch_bias_ecoli_2.cfg %s" % out_file )

  xs, ys = helpers.series_from_pipeline_result( fh=open(out_file), y='reference_bias', bias_report=True)
  xs, ys_err = helpers.series_from_pipeline_result( fh=open(out_file), y='error_bias', bias_report=True)
  xs = [ x * 100 for x in xs ]

  if include_zero:
    plt.xlim( xmin=-5, xmax=110 )
  else:
    plt.xlim( xmin=5, xmax=110 )

  #print xs, ys

  # weird crap
  ys = [ y * 100 + 1e-7 for y in ys ]
  ys_err = [ y * 100 + 1e-7 for y in ys_err ]
  width = 100. / len(ys) * 0.75 / 2

  #plt.bar( left=xs, height=ys, label='Bias', width=width )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlim(xmax=105)
  if include_zero:
    ax.bar( left=[x for x in xs], height=ys, label='Reference bias (%)', width=width, color='b', log=False, align='center')
    ax.bar( left=[x + width/2 for x in xs], height=ys_err, label='Error bias (%)', width=width, color='b', log=False, align='center')
  else:
    ax.bar( left=[x for x in xs][1:], height=ys[1:], label='Reference bias (%)', width=width, color='b', log=False, align='center')
    ax.bar( left=[x + width for x in xs][1:], height=ys_err[1:], label='Error bias (%)', width=width, color='r', log=False, align='center')
    ax.set_xlim(xmin=5)
  ax.set_ylabel('SNPs (%)')
  ax.set_xlabel('Bias (%)')
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/reference-bias-profile-with-errors.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def plot_reference_bias( out_file, name ):
  xs, ys = helpers.series_from_pipeline_batch( fh=open(out_file), x='mult_snp_prob', y='mean_reference', bias_report=True, has_no_variation=False)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot( xs, ys )
  fig.savefig('%s/reference-bias-repeat-%s.pdf' % (REPORT_DIRECTORY, name), format='pdf', dpi=1000)
  

def plot_read_length_vs_alignment_ecoli():
  out_file = "out/pipeline_batch_many_read_lengths_ecoli_20141119.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'read_length' ) )
    y.append( find_column( line, 'read_f1' ) )
    y_snp.append( find_column( line, 'vcf_f1' ) )

  #print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Alignment', color='r')
  ax.plot(x, y_snp, label='SNP Call', color='b')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Read Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-read-length-vs-alignment-ecoli.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_read_length_vs_alignment_ecoli_detailed():
  out_file = "out/pipeline_batch_many_read_lengths_ecoli_150107.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'read_length' ) )
    y.append( find_column( line, 'read_f1' ) )
    y_snp.append( find_column( line, 'vcf_f1' ) )

  #print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Alignment', color='r')
  ax.plot(x, y_snp, label='SNP Call', color='b')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Read Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-read-length-vs-alignment-ecoli-detailed.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_coverage_vs_alignment_circoviridae():
  out_file = "out/pipeline_batch_many_read_coverages_circoviridae_141127.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'coverage' ) )
    y.append( find_column( line, 'read_f1' ) )
    y_snp.append( find_column( line, 'vcf_f1' ) )

  #print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Alignment', color='r')
  ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Coverage')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/circoviridae-coverage-vs-alignment-no-errors.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_coverage_vs_alignment_circoviridae_errors():
  out_file = "out/pipeline_batch_many_read_coverages_errors_circoviridae_141127.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'coverage' ) )
    y.append( find_column( line, 'read_f1' ) )
    y_snp.append( find_column( line, 'vcf_f1' ) )

  #print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Alignment', color='r')
  ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Coverage')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/circoviridae-coverage-vs-alignment-with-errors.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_coverage_vs_alignment_circoviridae_high_errors():
  out_file = "out/pipeline_batch_many_read_coverages_high_errors_circoviridae_141127.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'coverage' ) )
    y.append( find_column( line, 'read_f1' ) )
    y_snp.append( find_column( line, 'vcf_f1' ) )

  #print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Alignment', color='r')
  ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Coverage')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/circoviridae-coverage-vs-alignment-high-errors.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_coverage_vs_alignment_ecoli_high_errors():
  out_file = "out/pipeline_batch_many_read_coverages_high_errors_ecoli_141127.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'coverage' ) )
    y.append( find_column( line, 'read_f1' ) )
    y_snp.append( find_column( line, 'vcf_f1' ) )

  #print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Alignment', color='r')
  ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Coverage')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-coverage-vs-alignment-high-errors.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )


def plot_mutation_vs_bias_ecoli():
  out_file = "out/pipeline_batch_bias_vs_mutation_ecoli_141119.out"
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_snp = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    #total_snps = find_column( line, 'total_variations' )
    x.append( find_parameter( line, 'snp_prob' ) )
    reference_bias_hist = find_column( line, 'reference_bias' )
    invisible = reference_bias_hist.split('|')[-1]
    percent_invisible = float(invisible)
    y.append( percent_invisible )

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Invisible SNPs (%)', color='r')
  ax.set_ylabel('%')
  ax.set_xlabel('Mutation rate (%)')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutation-vs-invisible-ecoli.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_alignment_circoviridae():
  out_file = "out/insert_experiment_141206b.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'max_insertion_len' ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )
    #y_snp.append( find_column( line, 'vcf_f1' ) )

  print "x", x
  print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r')
  ax.plot(x, y_precision, label='Precision', color='g')
  ax.plot(x, y_recall, label='Recall', color='b')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('%')
  ax.set_xlabel('Insertion Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/circoviridae-insertion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_alignment_ecoli():
  out_file = "out/insert_experiment_141209.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'max_insertion_len' ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )

  print "x", x
  print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r')
  ax.plot(x, y_precision, label='Precision', color='g')
  ax.plot(x, y_recall, label='Recall', color='b')
  ax.set_ylabel('%')
  ax.set_xlabel('Insertion Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-insertion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_deletion_vs_alignment_circoviridae():
  out_file = "out/delete_experiment_circoviridae_141209.out"
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'max_deletion_len' ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )
    #y_snp.append( find_column( line, 'vcf_f1' ) )

  print "x", x
  print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r')
  ax.plot(x, y_precision, label='Precision', color='g')
  ax.plot(x, y_recall, label='Recall', color='b')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('%')
  ax.set_xlabel('Deletion Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/circoviridae-deletion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_readlength_circoviridae():
  out_file = "out/insert_readlength_150108.out"
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'read_length' ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )
    #y_snp.append( find_column( line, 'vcf_f1' ) )

  print "x", x
  print "y", y
  print "y_p", y_precision
  print "y_r", y_recall
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r')
  ax.plot(x, y_precision, label='Precision', color='g')
  ax.plot(x, y_recall, label='Recall', color='b')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('%')
  ax.set_xlabel('Read Length')
  ax.set_title( 'Alignment accuracy with insertion of length 20' )
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/circoviridae-insertion-vs-readlength.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_readlength_ecoli():
  out_file = "out/insert_readlength_ecoli_150108a.out"
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'read_length' ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )
    #y_snp.append( find_column( line, 'vcf_f1' ) )

  print "x", x
  print "y", y
  print "y_p", y_precision
  print "y_r", y_recall
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r')
  ax.plot(x, y_precision, label='Precision', color='g')
  ax.plot(x, y_recall, label='Recall', color='b')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('%')
  ax.set_xlabel('Read Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-insertion-vs-readlength.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )


def plot_deletion_vs_alignment_ecoli():
  out_file = "out/delete_experiment_ecoli_141211.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'max_deletion_len' ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )
    #y_snp.append( find_column( line, 'vcf_f1' ) )

  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r')
  ax.plot(x, y_precision, label='Precision', color='g')
  ax.plot(x, y_recall, label='Recall', color='b')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('%')
  ax.set_xlabel('Deletion Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-deletion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_vcf_parent_vs_child_chromosomes( parent_fn, child_fn, true_fn, out_fn, ymax=None ):
  parent = bio.MultiChromosomeVCF( reader=open( '../../data/%s' % parent_fn ) )
  child = bio.MultiChromosomeVCF( reader=open( '../../data/%s' % child_fn ) )
  truth = bio.MultiChromosomeVCF( reader=open( '../../data/%s' % true_fn ) )

  result = bio.multi_chromosome_comparison( parent, child, truth )

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.scatter(result['x_fp'], result['y_fp'], color='r', label='Present only in child (%i)' % len(result['x_fp']), edgecolor = "none")#, log=True)
  ax.scatter(result['x_tp'], result['y_tp'], color='g', label='Present in both (%i)' % len(result['x_tp']), edgecolor="none", alpha=0.2)#, log=True)
  ax.scatter(result['x_fn'], result['y_fn'], color='b', label='Present only in parent (%i)' % len(result['x_fn']), edgecolor = "none", alpha=0.4)#, log=True)
  ax.scatter(result['x_t'], result['y_t'], color='#000000', marker='s', s=80, facecolors='none')#, log=True)

  ax.set_ylabel('Quality')
  ax.set_xlabel('Position')
  ax.set_xlim( xmin=0, xmax=result['offset'] )
  if ymax is None:
    ax.set_ylim( ymin=0 )#, ymax=60 )
  else:
    ax.set_ylim( ymin=0, ymax=ymax )
  leg = ax.legend(loc='upper left', bbox_to_anchor=(0., 0.9), prop={'size':8})
  leg.get_frame().set_alpha(0.8)

  ymin, ymax = ax.get_ylim()
  ax.vlines(x=result['chromosome_offsets'][:-1], ymin=ymin, ymax=ymax-0.1, color='#909090')
  prop = ymax / 30
  for item in result['chromosome_names']:
    ax.text(item[0], ymax - prop - prop * item[1], item[2], color='k', fontsize=6 )

  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

  
def plot_vcf_parent_vs_child( parent_fn, child_fn ):
  parent = bio.VCF( reader=open( '../../data/%s' % parent_fn ) )
  child = bio.VCF( reader=open( '../../data/%s' % child_fn ) )
  diff = bio.VCFDiff( vcf_correct=parent, vcf_candidate=child, generate_positions = True )
  
  x_tp, y_tp = diff.confidence_map( 'tp' )
  x_fp, y_fp = diff.confidence_map( 'fp' )
  x_fn, y_fn = diff.confidence_map( 'fn' )

  fig = plt.figure()
  ax = fig.add_subplot(111)
  
  ax.scatter(x_fp, y_fp, color='r', label='Present only in child (%i)' % len(x_fp), edgecolor = "none")#, log=True)
  ax.scatter(x_tp, y_tp, color='g', label='Present in both (%i)' % len(x_tp), edgecolor = "none")#, log=True)
  ax.scatter(x_fn, y_fn, color='b', label='Present only in parent (%i)' % len(x_fn), edgecolor = "none")#, log=True)

  ax.set_ylabel('Quality')
  ax.set_xlabel('Position')
  ax.set_xlim( xmin=0 )
  leg = ax.legend(loc='upper left', prop={'size':10})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

from matplotlib.ticker import FuncFormatter
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    #s = str(10 * y)
    s = '%.0f' % (1000 * y)
    return s

    # The percent symbol needs escaping in latex
    #if rcParams['text.usetex'] == True:
    #    return s + r'$\%$'
    #else:
    #    return s + '%'

def autolabel(ax, rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        #ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height*1000), ha='center', va='bottom')
        ax.text(rect.get_x()+rect.get_width()/2., height + 0., '%.1f'%(height*1000), ha='center', va='bottom')

def plot_mappability( src, name, short_name ):
  lines = open( src ).readlines()
  accuracy_list = lines[1].split(':')[1]
  accuracy = [float(x) * 100 for x in accuracy_list.strip( ' []\n' ).split(',')]

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(range(0, len(accuracy)), accuracy, s=2, label=name, color='b' )#, s=2)
  ax.set_ylabel('Mapping Accuracy')
  ax.set_xlabel('Position')
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  plt.xlim(xmin=0, xmax=8204)
  plt.ylim(ymax=100.01, ymin=0.)
  fig.savefig('%s/mappability-circoviridae-%s-50.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

def plot_bias( src, name, short_name ):
  lines = open( src ).readlines()
  bias_list = lines[2].split(':')[1]
  bias = [float(x) * 100 for x in bias_list.strip( ' []\n' ).split(',')]

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(range(0, len(bias)), bias, s=2, label=name, color='b' )#, s=2)
  ax.set_ylabel('Bias')
  ax.set_xlabel('Position')
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  plt.xlim(xmin=0, xmax=8204)
  plt.ylim(ymax=100.01) #, ymin=0.)
  fig.savefig('%s/bias-circoviridae-%s-50.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

def plot_mappability_comparison():
  lines = open( 'out/mappability_circoviridae_bwa.out' ).readlines()
  accuracy_list = lines[1].split(':')[1]
  accuracy_bwa = [float(x) * 100 for x in accuracy_list.strip( ' []\n' ).split(',')]

  lines = open( 'out/mappability_circoviridae_bwa_sw.out' ).readlines()
  accuracy_list = lines[1].split(':')[1]
  accuracy_bwa_sw = [float(x) * 100 for x in accuracy_list.strip( ' []\n' ).split(',')]

  lines = open( 'out/mappability_circoviridae_bowtie.out' ).readlines()
  accuracy_list = lines[1].split(':')[1]
  accuracy_bowtie = [float(x) * 100 for x in accuracy_list.strip( ' []\n' ).split(',')]

  fig = plt.figure()
  ax = fig.add_subplot(111)

  #ax.plot(accuracy)
  ax.scatter(range(0, len(accuracy_bwa)), accuracy_bwa, s=2, label="BWA", color='g' )#, s=2)
  ax.scatter(range(0, len(accuracy_bowtie)), accuracy_bowtie, s=2, label="Bowtie", color='b' )#, s=2)
  ax.scatter(range(0, len(accuracy_bwa_sw)), accuracy_bowtie, s=2, label="BWA-SW", color='r' )#, s=2)
  ax.set_ylabel('Mapping Accuracy')
  ax.set_xlabel('Position')
  ax.set_title( 'Mappability by genome position' )
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  plt.xlim(xmin=0, xmax=8204)
  plt.ylim(ymax=100.01, ymin=0.)
  fig.savefig('%s/mappability-circoviridae-comparison-50.pdf' % (REPORT_DIRECTORY), format='pdf', dpi=1000)

def plot_mappability_hist( src, mapper ):
  #lines = open( 'results.out' ).readlines()
  lines = open( src ).readlines()
  accuracy_list = lines[1].split(':')[1]
  accuracy = [float(x) * 100 for x in accuracy_list.strip( ' []\n' ).split(',')]
  #print len(n), len(bins)
  
  fig = plt.figure()
  ax = fig.add_subplot(111)
  #ax.plot( bins, n )
  bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
  #n, bins, patches = ax.hist(accuracy, 10, normed=1)#, normed=1, histtype='stepfilled')
  n, bins, patches = ax.hist(accuracy, bins=bins, normed=1)#, normed=1, histtype='stepfilled')
  formatter = FuncFormatter(to_percent)
  #hist = np.histogram( accuracy, bins=50, density=True ) # sums to read length
  plt.gca().yaxis.set_major_formatter(formatter)
  ax.set_ylabel('%')
  ax.set_xlabel('Mapping Accuracy')
  ax.set_title( mapper )
  ax.set_xlim(xmin=0, xmax=100)
  ax.set_ylim(ymax=70./1000)
  autolabel( ax, patches ) 
  fig.savefig('%s/mappability-hist-circoviridae-repeated-%s-50.pdf' % ( REPORT_DIRECTORY, mapper ), format='pdf', dpi=1000)

def plot_bias_hist( src, mapper ):
  #lines = open( 'results.out' ).readlines()
  lines = open( src ).readlines()
  bias_list = lines[2].split(':')[1]
  bias = [float(x) * 100 for x in bias_list.strip( ' []\n' ).split(',')]
  #print len(bias)
  #bias = [ x for x in bias if x != 0.0 ]
  #print len(bias)
  
  fig = plt.figure()
  ax = fig.add_subplot(111)
  #ax.plot( bins, n )
  #n, bins, patches = ax.hist(accuracy, 10, normed=1)#, normed=1, histtype='stepfilled')
  bins = [-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, -0.01, +0.01, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
  n, bins, patches = ax.hist(bias, bins=bins, normed=1)#, normed=1, histtype='stepfilled')
  #print n, bins, patches
  print "*** Bias for %s with %s: %s ***" % ( src, mapper, bio.sign_bucket( bias ) )
  values = bio.bucket( bias, bins )
  for i in xrange(len(bins)):
    print "%.2f\t%i" % ( bins[i], values[i] )
  formatter = FuncFormatter(to_percent)
  #hist = np.histogram( accuracy, bins=50, density=True ) # sums to read length
  plt.gca().yaxis.set_major_formatter(formatter)
  ax.set_ylabel('%')
  ax.set_xlabel('Bias')
  ax.set_title( mapper )
  ax.set_xlim(xmin=-100, xmax=100)
  ax.set_ylim(ymax=100./1000)
  autolabel( ax, patches ) 
  fig.savefig('%s/bias-hist-circoviridae-repeated-%s-50.pdf' % ( REPORT_DIRECTORY, mapper ), format='pdf', dpi=1000)

# no longer used
# TODO ecoli-mutations-snps-unmapped-2.pdf (not used)
# TODO mutation-f1-hiv.pdf (not used)
# TODO ecoli-mutations-snps-f1-2.pdf (not used)
#plot_read_length_vs_alignment_ecoli()

##### in report
#mutation_hiv() # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_hiv_snp_f1() # mutation-f1-snp-hiv.pdf
#mutation_ecoli_snp_f1() # ecoli-mutations-snps-2.pdf
#mutation_ecoli_snp_histogram() # ecoli-mutations-snps-2.pdf
#coverage_ecoli() # ecoli-coverage-snp.pdf
#read_length_ecoli() # ecoli-read-length-bowtie-snp.pdf
#read_length_ecoli_map() # ecoli-read-length-bowtie-mapping.pdf
#read_length_hiv() # hiv-read-length-vcf-snps.pdf (lots of noise)
#snp_vs_map_ecoli() # ecoli-map_vs_snp-400.pdf plot of f1-scores of snps found and mapping accuracy, vs mutation rate, for e-coli
#plot_coverage_vs_alignment_ecoli_high_errors()
#plot_reference_bias_ecoli_example() # 
#plot_read_length_vs_alignment_ecoli_detailed()

##### experimental
#read_length_ecoli_map_low_mutation() # ecoli-read-length-bowtie-mapping.pdf
#plot_mutation_vs_bias_ecoli()
#plot_reference_bias_ecoli_example_with_errors() # 
#plot_insertion_vs_alignment_circoviridae()
#plot_insertion_vs_alignment_ecoli()
#plot_deletion_vs_alignment_ecoli()
#plot_deletion_vs_alignment_circoviridae()
#plot_insertion_vs_readlength_circoviridae()
#plot_insertion_vs_readlength_ecoli()
#old plot_vcf_parent_vs_child()
#plot_vcf_parent_vs_child( 'Plasmodium_falciparum_3d_p_CHR02_recovered.vcf', 'Plasmodium_falciparum_3d7_1q_CHR02_recovered.vcf', 'malaria-3d7-p-vs-3d7-1q-vcf' )
#plot_vcf_parent_vs_child_chromosomes( 'sample_multi.vcf', 'sample_multi_recovered.vcf', 'sample_multi' )
#plot_vcf_parent_vs_child_chromosomes( 'Plasmodium_falciparum_3d_p_resequenced_recovered.vcf', 'Plasmodium_falciparum_3d_1q_resequenced_recovered.vcf', 'Plasmodium_falciparum_3d7_1q.vcf', 'malaria_vcf_1q_snp' )
#plot_vcf_parent_vs_child_chromosomes( 'freebayes_3d_p_filtered_normed.vcf', 'freebayes_1q_filtered_normed.vcf', 'Plasmodium_falciparum_3d7_1q.vcf', 'freebayes_malaria_vcf_1q_snp', ymax=120 )
#plot_vcf_parent_vs_child_chromosomes( 'Plasmodium_falciparum_3d_p_resequenced_recovered.vcf', 'Plasmodium_falciparum_1a_resequenced_recovered.vcf', 'Plasmodium_falciparum_3d7_1a.vcf', 'malaria_vcf_1a_snp' )
#plot_vcf_parent_vs_child_chromosomes( 'freebayes_3d_p_filtered_normed.vcf', 'freebayes_1a_normed.vcf', 'Plasmodium_falciparum_3d7_1a.vcf', 'freebayes_malaria_vcf_1a_snp', ymax=50 )
#plot_mappability_comparison()
#plot_mappability( 'out/mappability_circoviridae_bowtie.out', 'Bowtie', 'bowtie' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw.out', 'BWASW', 'bwasw' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_snp.out', 'BWASW with SNP', 'bwasw-snp' )
#-- plot_mappability( 'out/mappability_circoviridae_bwa_sw_insert_1.out', 'BWASW with insertion', 'bwasw-insertion' )
#-- plot_mappability( 'out/mappability_circoviridae_bwa_sw_delete_1.out', 'BWASW with deletion', 'bwasw-deletion' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_insert_10.out', 'BWASW with 10 insertion', 'bwasw-insertion-10' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_insert_20.out', 'BWASW with 20 insertion', 'bwasw-insertion-20' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_insert_30.out', 'BWASW with 30 insertion', 'bwasw-insertion-30' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_insert_40.out', 'BWASW with 40 insertion', 'bwasw-insertion-40' )
#
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_delete_10.out', 'BWASW with 10 deletion', 'bwasw-deletion-10' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_delete_20.out', 'BWASW with 20 deletion', 'bwasw-deletion-20' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_delete_30.out', 'BWASW with 30 deletion', 'bwasw-deletion-30' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_delete_40.out', 'BWASW with 40 deletion', 'bwasw-deletion-40' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_snp_error_snp.out', 'BWASW with SNP and 1 subst error', 'bwasw-snp-subst-error' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_snp_error_snp_2.out', 'BWASW with SNP and 2 subst error', 'bwasw-snp-subst-2-error' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_snp_error_snp_5.out', 'BWASW with SNP and 5 subst error', 'bwasw-snp-subst-5-error' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_snp_error_ins_1.out', 'BWASW with SNP and 1 insert error', 'bwasw-snp-ins-1-error' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_error_ins_1.out', 'BWASW with 1 insert error', 'bwasw-ins-1-error' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_error_ins_1_test.out', 'BWASW with 1 insert error', 'bwasw-ins-1-error-test' )
#plot_mappability( 'out/mappability_circoviridae_bwa_sw_error_del_1.out', 'BWASW with 1 delete error', 'bwasw-del-1-error' )
#plot_mappability( 'out/mappability_circoviridae_bowtie_error_ins_1.out', 'Bowtie with 1 insert error', 'bowtie-ins-1-error' )
#
#-- plot_mappability_hist( 'out/mappability_circoviridae_bwa.out', 'bwa' )
#plot_mappability_hist( 'out/mappability_circoviridae_bowtie.out', 'bowtie' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw.out', 'bwasw' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_snp.out', 'bwasw-snp' )
#-- plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_insert_1.out', 'bwasw-insertion' )
#-- plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_delete_1.out', 'bwasw-deletion' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_insert_10.out', 'bwasw-insertion-10' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_insert_20.out', 'bwasw-insertion-20' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_insert_30.out', 'bwasw-insertion-30' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_insert_40.out', 'bwasw-insertion-40' )
#
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_delete_10.out', 'bwasw-deletion-10' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_delete_20.out', 'bwasw-deletion-20' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_delete_30.out', 'bwasw-deletion-30' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_delete_40.out', 'bwasw-deletion-40' )
#
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_snp_error_snp.out', 'bwasw-snp-subst-error' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_snp_error_snp_2.out', 'bwasw-snp-subst-2-error' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_snp_error_snp_5.out', 'bwasw-snp-subst-5-error' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_snp_error_ins_1.out', 'bwasw-snp-ins-1-error' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_error_ins_1.out', 'bwasw-ins-1-error' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_error_ins_1_test.out', 'bwasw-ins-1-error-test' )
#plot_mappability_hist( 'out/mappability_circoviridae_bowtie_error_ins_1.out', 'bowtie-ins-1-error' )
#plot_mappability_hist( 'out/mappability_circoviridae_bwa_sw_error_del_1.out', 'bwasw-del-1-error' )
#
#plot_bias( 'out/bias_circoviridae_bwa_sw_snp_error_snp_1.out', 'BWASW with SNP and 1 subst error', 'bwasw-snp-subst-error' )
#plot_bias_hist( 'out/bias_circoviridae_bwa_sw_snp_error_snp_1.out', 'bwasw-snp-subst-error' )

#plot_reference_bias( 'out/hiv_reference_bias_vs_repeats_150212.out', 'hiv' )
plot_reference_bias( 'out/ecoli_reference_bias_vs_repeats_150212.out', 'ecoli' )
