#
# build the charts used in the final report
# generally these assume the analysis has been done; the command to do this is included in the function
#

import collections
import datetime
import math
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab
import sys

import bio
import helpers

REPORT_DIRECTORY = '../../report/plots'

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
pylab.rcParams['figure.figsize'] = 6, 4

MALARIA_CHROMOSOME_MAP = {
  'Pf3D7_01_v3': 'psu|MAL1',
  'Pf3D7_02_v3': 'psu|MAL2',
  'Pf3D7_03_v3': 'psu|MAL3',
  'Pf3D7_04_v3': 'psu|MAL4',
  'Pf3D7_05_v3': 'psu|MAL5',
  'Pf3D7_07_v3': 'psu|MAL7',
  'Pf3D7_06_v3': 'psu|MAL6',
  'Pf3D7_08_v3': 'psu|MAL8',
  'Pf3D7_09_v3': 'psu|MAL9',
  'Pf3D7_10_v3': 'psu|MAL10',
  'Pf3D7_11_v3': 'psu|MAL11',
  'Pf3D7_12_v3': 'psu|MAL12',
  'Pf3D7_13_v3': 'psu|MAL13',
  'Pf3D7_14_v3': 'psu|MAL14',
}


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
def mutation_hiv( out_file ):
  #out_file = "out/pipeline_batch_mutation_hiv_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
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
  ax.plot(x, y_base, label='Baseline', color='b', marker='s')
  ax.plot(x, y_bow, label='Bowtie2', color='g', marker='^')
  ax.plot(x, y_bwa, label='BWA', color='r', marker='o')
  ax.set_ylabel('Unmapped %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/mutation-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )
  
def mutation_hiv_snp_f1(out_file):
  #out_file = "out/pipeline_batch_mutation_hiv_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
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
  ax.plot(x, y_base, label='Baseline', color='b', marker='s')
  ax.plot(x, y_bow, label='Bowtie2', color='g', marker='^')
  ax.plot(x, y_bwa, label='BWA', color='r', marker='o')
  ax.set_ylabel('F1-Score')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/mutation-f1-snp-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def mutation_ecoli_snp_f1(out_file):
  #out_file = "out/pipeline_batch_mutation_ecoli_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
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
  ax.plot(x, y_bow, label='Bowtie2', color='g', marker='^')
  ax.plot(x, y_bwa, label='BWA', color='r', marker='o')
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

def coverage_ecoli(out_file):
  #out_file = "out/pipeline_batch_coverage_ecoli_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
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
  ax.plot(x, b5[1:], label='5x coverage', color='b', marker='s')
  ax.plot(x, b10[1:], label='10x coverage', color='g', marker='^')
  ax.plot(x, b20[1:], label='20x coverage', color='r', marker='o')
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

def snp_vs_map_ecoli(out_file):
  #out_file = "out/pipeline_batch_read_length_ecoli_%s_inchworm.out" % datetime.datetime.now().strftime("%Y%m%d")
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
  ax.plot(x, m400, label='Align', color='b', marker='s')
  ax.plot(x, s400, label='SNV', color='g', marker='^')
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
  ax.plot(x, m100, label='Align', color='b', marker='s')
  ax.plot(x, s100, label='SNV', color='g', marker='^')
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
    ax.bar( left=[x for x in xs], height=ys, label='SNVs affected (%)', width=width, color='b', log=False, align='center')
  else:
    ax.bar( left=[x for x in xs][1:], height=ys[1:], label='SNVs affected (%)', width=width, color='b', log=False, align='center')
    ax.set_xlim(xmin=5)
  ax.set_ylabel('SNVs (%)')
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
  ax.set_ylabel('SNVs (%)')
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
  ax.plot(x, y, label='Alignment', color='r', marker='o')
  ax.plot(x, y_snp, label='SNV Call', color='b', marker='^')
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
  ax.plot(x, y, label='Alignment', color='r', marker='o')
  ax.plot(x, y_snp, label='SNV Call', color='b', marker='s')
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
  ax.plot(x, y, label='Alignment', color='r', marker='o')
  ax.plot(x, y_snp, label='SNV Call', color='b', marker='s')
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
  ax.plot(x, y, label='Invisible SNVs (%)', color='r')
  ax.set_ylabel('%')
  ax.set_xlabel('Mutation rate (%)')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutation-vs-invisible-ecoli.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

#def plot_insertion_vs_alignment_circoviridae():
#  out_file = "out/insert_experiment_141206b.out"
#
#  bio.log_stderr( 'extracting values from %s...' % out_file )
#  x = []
#  y = []
#  y_precision = []
#  y_recall = []
#  fh = open( out_file, 'r' )
#  for line in fh:
#    if line.startswith( '#' ):
#      continue
#    x.append( find_parameter( line, 'max_insertion_len' ) )
#    y.append( find_column( line, 'vcf_f1' ) )
#    y_precision.append( find_column( line, 'vcf_precision' ) )
#    y_recall.append( find_column( line, 'vcf_recall' ) )
#    #y_snp.append( find_column( line, 'vcf_f1' ) )
#
#  print "x", x
#  print "y", y
#  # draw graph
#  fig = plt.figure()
#  ax = fig.add_subplot(111)
#  ax.plot(x, y, label='F1', color='r')
#  ax.plot(x, y_precision, label='Precision', color='g')
#  ax.plot(x, y_recall, label='Recall', color='b')
#  #ax.plot(x, y_snp, label='SNV Call', color='b')
#  ax.set_ylabel('%')
#  ax.set_xlabel('Insertion Length')
#  #ax.set_ylim(ymin=90)
#  leg = ax.legend(loc='upper right', prop={'size':12})
#  leg.get_frame().set_alpha(0.8)
#  fig.savefig('%s/circoviridae-insertion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
#  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_alignment( out_file, x_axis, target_file, x_field ):
  #out_file = "out/insert_experiment_141209.out"

  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, x_field ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )

  print "x", x
  print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r', marker='o')
  ax.plot(x, y_precision, label='Precision', color='g', marker='^')
  ax.plot(x, y_recall, label='Recall', color='b', marker='s')
  ax.set_ylabel('Accuracy (%)')
  #ax.set_xlabel('Insertion Length')
  ax.set_xlabel(x_axis)
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  #fig.savefig('%s/ecoli-insertion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file ), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

#def plot_deletion_vs_alignment_circoviridae():
#  out_file = "out/delete_experiment_circoviridae_141209.out"
#  bio.log_stderr( 'extracting values from %s...' % out_file )
#  x = []
#  y = []
#  y_precision = []
#  y_recall = []
#  fh = open( out_file, 'r' )
#  for line in fh:
#    if line.startswith( '#' ):
#      continue
#    x.append( find_parameter( line, 'max_deletion_len' ) )
#    y.append( find_column( line, 'vcf_f1' ) )
#    y_precision.append( find_column( line, 'vcf_precision' ) )
#    y_recall.append( find_column( line, 'vcf_recall' ) )
#    #y_snp.append( find_column( line, 'vcf_f1' ) )
#
#  print "x", x
#  print "y", y
#  # draw graph
#  fig = plt.figure()
#  ax = fig.add_subplot(111)
#  ax.plot(x, y, label='F1', color='r')
#  ax.plot(x, y_precision, label='Precision', color='g')
#  ax.plot(x, y_recall, label='Recall', color='b')
#  #ax.plot(x, y_snp, label='SNV Call', color='b')
#  ax.set_ylabel('%')
#  ax.set_xlabel('Deletion Length')
#  #ax.set_ylim(ymin=90)
#  leg = ax.legend(loc='upper right', prop={'size':12})
#  leg.get_frame().set_alpha(0.8)
#  fig.savefig('%s/circoviridae-deletion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
#  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_variation_distance( out_file, x_col, x_name, title, report_fn, as_percent=False, legend='upper left' ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y = []
  y_precision = []
  y_recall = []
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( float( find_parameter( line, x_col ) ) )
    y.append( find_column( line, 'vcf_f1' ) )
    y_precision.append( find_column( line, 'vcf_precision' ) )
    y_recall.append( find_column( line, 'vcf_recall' ) )
    #y_snp.append( find_column( line, 'vcf_f1' ) )

  if as_percent:
    x = [ el * 100 for el in x ]

  print "x", x
  print "y", y
  print "y_p", y_precision
  print "y_r", y_recall
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='F1', color='r', marker='o')
  ax.plot(x, y_precision, label='Precision', color='g', marker='^')
  ax.plot(x, y_recall, label='Recall', color='b', marker='s')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('Accuracy (%)')
  ax.set_xlabel(x_name)
  #ax.set_title( title)
  #ax.set_ylim(ymin=90)
  ax.set_ylim(ymin=0)
  leg = ax.legend(loc=legend, prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, report_fn ), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_insertion_vs_readlength( out_file, target_file ):
  #out_file = "out/insert_readlength_150108.out"
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
  ax.plot(x, y, label='F1', color='r', marker='o')
  ax.plot(x, y_precision, label='Precision', color='g', marker='^')
  ax.plot(x, y_recall, label='Recall', color='b', marker='s')
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel('Accuracy (%)')
  ax.set_xlabel('Read Length')
  #ax.set_title( 'Alignment accuracy with insertion of length 20' )
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='upper left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

#def plot_insertion_vs_readlength_ecoli():
#  out_file = "out/insert_readlength_ecoli_150108a.out"
#  bio.log_stderr( 'extracting values from %s...' % out_file )
#  x = []
#  y = []
#  y_precision = []
#  y_recall = []
#  fh = open( out_file, 'r' )
#  for line in fh:
#    if line.startswith( '#' ):
#      continue
#    x.append( find_parameter( line, 'read_length' ) )
#    y.append( find_column( line, 'vcf_f1' ) )
#    y_precision.append( find_column( line, 'vcf_precision' ) )
#    y_recall.append( find_column( line, 'vcf_recall' ) )
#    #y_snp.append( find_column( line, 'vcf_f1' ) )
#
#  print "x", x
#  print "y", y
#  print "y_p", y_precision
#  print "y_r", y_recall
#  # draw graph
#  fig = plt.figure()
#  ax = fig.add_subplot(111)
#  ax.plot(x, y, label='F1', color='r')
#  ax.plot(x, y_precision, label='Precision', color='g')
#  ax.plot(x, y_recall, label='Recall', color='b')
#  #ax.plot(x, y_snp, label='SNV Call', color='b')
#  ax.set_ylabel('%')
#  ax.set_xlabel('Read Length')
#  #ax.set_ylim(ymin=90)
#  leg = ax.legend(loc='upper left', prop={'size':12})
#  leg.get_frame().set_alpha(0.8)
#  fig.savefig('%s/ecoli-insertion-vs-readlength.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
#  bio.log_stderr( 'extracting values from %s: done' % out_file )

#def plot_deletion_vs_alignment_ecoli():
#  out_file = "out/delete_experiment_ecoli_141211.out"
#
#  bio.log_stderr( 'extracting values from %s...' % out_file )
#  x = []
#  y = []
#  y_precision = []
#  y_recall = []
#  fh = open( out_file, 'r' )
#  for line in fh:
#    if line.startswith( '#' ):
#      continue
#    x.append( find_parameter( line, 'max_deletion_len' ) )
#    y.append( find_column( line, 'vcf_f1' ) )
#    y_precision.append( find_column( line, 'vcf_precision' ) )
#    y_recall.append( find_column( line, 'vcf_recall' ) )
#    #y_snp.append( find_column( line, 'vcf_f1' ) )
#
#  # draw graph
#  fig = plt.figure()
#  ax = fig.add_subplot(111)
#  ax.plot(x, y, label='F1', color='r')
#  ax.plot(x, y_precision, label='Precision', color='g')
#  ax.plot(x, y_recall, label='Recall', color='b')
#  #ax.plot(x, y_snp, label='SNV Call', color='b')
#  ax.set_ylabel('%')
#  ax.set_xlabel('Deletion Length')
#  #ax.set_ylim(ymin=90)
#  leg = ax.legend(loc='upper right', prop={'size':12})
#  leg.get_frame().set_alpha(0.8)
#  fig.savefig('%s/ecoli-deletion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
#  bio.log_stderr( 'extracting values from %s: done' % out_file )

import re
def _simplify_chromosome_name( n ):
  x = re.sub( '.*_([^_]*)_.*', '\\1', n )
  return x

def plot_vcf_parent_vs_child_chromosomes( parent_fn, child_fn, true_fn, out_fn, ymax=None, legend_loc='upper left', legend_bbox=(0., 0.9), chromosome_map=None):
  parent = bio.MultiChromosomeVCF( reader=open( '../../data/%s' % parent_fn ) )
  child = bio.MultiChromosomeVCF( reader=open( '../../data/%s' % child_fn ) )
  truth = bio.MultiChromosomeVCF( reader=open( '../../data/%s' % true_fn ) )

  result = bio.multi_chromosome_comparison( parent, child, truth, parent_chromosome_map=chromosome_map )

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
  leg = ax.legend(loc=legend_loc, bbox_to_anchor=legend_bbox, prop={'size':8})
  leg.get_frame().set_alpha(0.8)

  ymin, ymax = ax.get_ylim()
  ax.vlines(x=result['chromosome_offsets'][:-1], ymin=ymin, ymax=ymax-0.1, color='#909090')
  prop = ymax / 30
  for item in result['chromosome_names']:
    ax.text(item[0], ymax - prop - prop * item[1], _simplify_chromosome_name( item[2] ), color='k', fontsize=6 )

  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

  
def plot_vcf( fn, out_fn ):
  vcf = bio.VCF( reader=open( '../../data/%s' % fn ) )
  x = []
  y = []
  for snp in vcf.snp_list:
    x.append( snp['pos'] )
    confidence = snp['conf']
    if confidence >= 1.:
      qual = MAX_QUALITY
    else:
      qual = -10 * math.log( 1. - confidence, 10 )
    y.append( qual )

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(x, y, color='g', edgecolor = "none")#, log=True)
  ax.set_ylabel('Quality')
  ax.set_xlabel('Position')
  ax.set_xlim( xmin=0 )
  #leg = ax.legend(loc='upper left', prop={'size':10})
  #leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

def plot_vcf_parent_vs_child( parent_fn, child_fn, out_fn, parent_label='parent', child_label='child', highlight_differences=False ):
  parent = bio.VCF( reader=open( '../../data/%s' % parent_fn ) )
  child = bio.VCF( reader=open( '../../data/%s' % child_fn ) )
  diff = bio.VCFDiff( vcf_correct=parent, vcf_candidate=child, generate_positions = True )
  
  x_tp, y_tp = diff.confidence_map( 'tp' )
  x_fp, y_fp = diff.confidence_map( 'fp' )
  x_fn, y_fn = diff.confidence_map( 'fn' )

  fig = plt.figure()
  ax = fig.add_subplot(111)
  
  ax.set_xlim( xmin=0, xmax=max(x_fp + x_tp + x_fn ) * 1.05 )
  ax.set_ylim( ymin=0, ymax=max(y_fp + y_tp + y_fn) * 1.05 )

  if highlight_differences:
    ax.scatter(x_fp, y_fp, color='r', label='Present only in %s (%i)' % (child_label, len(x_fp) ), edgecolor = "none")#, log=True)
    ax.scatter(x_tp, y_tp, color='g', label='Present in both (%i)' % len(x_tp), edgecolor = "none", alpha=0.2)#, log=True)
    ax.scatter(x_fn, y_fn, color='b', label='Present only in %s (%i)' % ( parent_label, len(x_fn) ), edgecolor = "none")#, log=True)
  else:
    ax.scatter(x_fp, y_fp, color='g', label='Present only in %s (%i)' % (child_label, len(x_fp) ), edgecolor = "none", alpha=0.4)#, log=True)
    ax.scatter(x_tp, y_tp, color='r', label='Present in both (%i)' % len(x_tp), edgecolor = "none")#, log=True)
    ax.scatter(x_fn, y_fn, color='b', label='Present only in %s (%i)' % ( parent_label, len(x_fn) ), edgecolor = "none", alpha=0.4)#, log=True)

  ax.set_ylabel('Quality')
  ax.set_xlabel('Position')
  ax.set_xlim( xmin=0 )
  leg = ax.legend(loc='upper right', prop={'size':8})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

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

def plot_all_mappability( infile, title, outfile ):
  plot_mappability( infile, title, outfile )
  plot_mappability_hist( infile, title, outfile )
  plot_bias( infile, title, outfile )
  plot_bias_hist( infile, title, outfile )

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
  plt.xlim(xmin=0, xmax=len(accuracy))#, xmax=8204)
  plt.ylim(ymax=100.01, ymin=0.)
  fig.savefig('%s/mappability-%s.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

def plot_entropy( src, short_name ):
  lines = open( src ).readlines()
  # find accuracy
  # find entropy
  for line in lines:
    if line.startswith('accuracy'):
      accuracy_list = [ float(x) for x in line.split(':')[1].strip( ' []\n' ).split(',')]
    if line.startswith('entropy'):
      entropy_list = [ float(x) for x in line.split(':')[1].strip( ' []\n' ).split(',') ]

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(entropy_list, accuracy_list, s=2, color='b' )#, s=2)
  fig.savefig('%s/entropy-%s.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

def plot_depth( src, short_name, range_start=0, range_stop=-1, variation_label='Variation', features=(), show_features=False ):
  '''
    plots either depth or evidence
  '''
  lines = open( src ).readlines()
  lists = []
  for line in lines:
    if line.startswith( '#' ):
      continue
    lists.append( [ float(x) for x in line.split(',') ] )

  # lists[0] = depth, lists[1] = breakpoints
  depths = lists[0]
  breakpoints = lists[1 + len(features)]
  fig = plt.figure()
  ax = fig.add_subplot(111)
  if range_stop < 0:
    range_stop = len(depths)
  first = True
  for breakpoint in breakpoints:
    if range_start <= breakpoint < range_stop:
      if first:
        first = False
        ax.axvline(x=breakpoint, color='r', alpha=0.1, label=variation_label)
      else:
        ax.axvline(x=breakpoint, color='r', alpha=0.1 )

  if not show_features: # do the depth
    plt.ylim(ymax=max(depths)+0.5, ymin=-0.5)
    ax.plot(range(range_start, range_stop), depths[range_start:range_stop], label='Coverage' ) #, s=2)

  # additional features
  else:
    for idx in xrange(0, len(features)):
      ax.plot(range(range_start, range_stop), lists[1 + idx][range_start:range_stop], label=features[idx] ) #, s=2)

  leg = ax.legend( loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/depth-%s.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

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
  plt.xlim(xmin=0, xmax=len(bias))#, xmax=8204)
  #plt.ylim(ymin=max(-100.01, min(bias)), ymax=min(100.01, max(bias))) #, ymin=0.)
  plt.ylim(ymin=-100.01, ymax=100.01) #, ymin=0.)
  fig.savefig('%s/bias-%s.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

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
  plt.xlim(xmin=0, xmax=max(len(accuracy_bwa)))#8204)
  plt.ylim(ymax=100.01, ymin=0.)
  fig.savefig('%s/mappability-circoviridae-comparison-50.pdf' % (REPORT_DIRECTORY), format='pdf', dpi=1000)

def plot_mappability_hist( src, title, fn ):
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
  ax.set_title( title )
  ax.set_xlim(xmin=0, xmax=100)
  ax.set_ylim(ymax=70./1000)
  autolabel( ax, patches ) 
  fig.savefig('%s/mappability-hist-%s.pdf' % ( REPORT_DIRECTORY, fn ), format='pdf', dpi=1000)

def plot_bias_hist( src, title, fn ):
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
  print "*** Bias for %s with %s: %s ***" % ( src, title, bio.sign_bucket( bias ) )
  values = bio.bucket( bias, bins )
  for i in xrange(len(bins)):
    print "%.2f\t%i" % ( bins[i], values[i] )
  formatter = FuncFormatter(to_percent)
  #hist = np.histogram( accuracy, bins=50, density=True ) # sums to read length
  plt.gca().yaxis.set_major_formatter(formatter)
  ax.set_ylabel('%')
  ax.set_xlabel('Bias')
  ax.set_title( title )
  ax.set_xlim(xmin=-100, xmax=100)
  ax.set_ylim(ymax=100./1000)
  autolabel( ax, patches ) 
  fig.savefig('%s/bias-hist-%s.pdf' % ( REPORT_DIRECTORY, fn ), format='pdf', dpi=1000)

##### migrated from ipython
def ecoli_mutations_snps_unmapped():
  # e-coli vs unmapped
  # from pipeline_batch_mutation_ecoli-140805.out
  x = np.linspace( 1, 10, 10 )
  ybase = [39.362935, 63.848565, 78.354124, 86.817800, 92.319459, 95.545825, 97.332601, 98.467877, 99.101603, 99.497479]
  ybow = [0.982743, 4.099772, 9.130041, 15.653223, 23.007691, 31.118041, 39.108282, 46.678516, 53.948746, 60.810065]
  ybwa = [0.430347, 2.065667, 5.465788, 10.584389, 17.202473, 25.366899, 34.136181, 43.066400, 51.827064, 60.169230]
  bowfit = np.poly1d(np.polyfit(x, ybow, 3))(x)
  bwafit = np.poly1d(np.polyfit(x, ybwa, 1))(x)
  yave = [ (z[0]+z[1])/2. for z in zip(ybow,ybwa) ]
  avefit = np.poly1d(np.polyfit(x, yave, 1))(x)
  
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.plot(x, ybase, label='Baseline', color='b', marker='s')
  ax.plot(x, ybow, label='Bowtie2', color='g', marker='^')
  ax.plot(x, ybwa, label='BWA', color='r', marker='o')
  ax.plot(x, avefit, '--', label='Linear Fit', color='#a0a0a0')
  #plt.plot(x, linear, label='Linear', color='#cccccc')
  ax.set_ylim(ymin=0.0)
  
  ax.set_ylabel('Unmapped %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutations-snps-unmapped-2.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def mutation_f1_hiv():
  # F1 Score with HIV mutations
  # from pipeline_batch_mutation_hiv.out
  x = np.linspace( 0, 10, 11 )
  ybase = [ 92.580019, 68.335878, 48.283262, 32.604374, 24.567328, 14.625393, 7.895958, 3.233477, 2.297751, 1.252409, 0.290557 ]
  ybow = [ 91.852570, 91.664642, 90.395481, 87.150838, 83.280085, 81.704669, 75.523673, 73.681143, 64.577398, 57.611241, 47.594936 ]
  ybwa = [ 92.143550, 91.960165, 91.565091, 90.469136, 88.239776, 82.909869, 78.275290, 72.542766, 69.081977, 59.929192, 51.646811 ]
  
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.plot(x, ybase, label='Baseline', marker='s')
  ax.plot(x, ybow, label='Bowtie2', marker='^')
  ax.plot(x, ybwa, label='BWA', marker='o')
  
  ax.set_ylabel('F1')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/mutation-f1-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def ecoli_mutations_snps_f1():
  # e-coli vs f1
  # from pipeline_batch_mutation_ecoli-140805.out
  x = np.linspace( 1, 10, 10 )
  ybase = [74.475096, 52.398053, 35.110156, 22.981784, 14.074148, 8.382036, 5.125477, 2.983206, 1.759230, 0.983296]
  ybow = [97.818423, 96.234812, 93.583440, 89.932263, 85.456786, 80.089299, 74.334322, 68.263784, 61.847011, 55.233382]
  ybwa = [98.421014, 97.596470, 95.854366, 93.101787, 89.298540, 84.248381, 78.309618, 71.517099, 64.036484, 56.116548]

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.plot(x, ybase, label='Baseline', marker='s')
  ax.plot(x, ybow, label='Bowtie2', marker='^')
  ax.plot(x, ybwa, label='BWA', marker='o')
  
  yave = [ (z[0]+z[1])/2. for z in zip(ybow,ybwa) ]
  avefit = np.poly1d(np.polyfit(x, yave, 1))(x)
  ax.set_ylim(ymax=100.)
  
  ax.plot(x, avefit, '--', label='Linear Fit', color='#a0a0a0')
  
  ax.set_ylabel('F1')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12}, loc='lower right')
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/ecoli-mutations-snps-f1-2.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def error_f1_hiv():
  x = np.linspace( 0, 10, 11 )
  # f1 scores with increasing read errors on HIV
  # from pipeline_batch_error_hiv_f1.out
  ybase = [ 91.076625, 69.418499, 49.078014, 31.958349, 21.278421, 13.363432, 7.452258, 4.360190, 2.583731, 1.919386, 0.290557 ]
  ybow = [ 91.852570, 91.081871, 90.360850, 87.757699, 84.312697, 80.472398, 76.288073, 68.799512, 62.196680, 59.044481, 52.192067 ]
  ybwa = [ 92.143550, 91.855093, 91.162790, 89.398853, 87.238387, 83.021890, 78.567462, 73.025161, 67.183142, 59.154930, 52.309830 ]
  
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.plot(x, ybase, label='Baseline', marker='s')
  ax.plot(x, ybow, label='Bowtie2', marker='^')
  ax.plot(x, ybwa, label='BWA', marker='o')
  
  ax.set_ylabel('F1')
  ax.set_xlabel('Error %')
  leg = ax.legend(prop={'size':12}, loc='lower right')
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/error-f1-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def error_unmapped_hiv():
  x = np.linspace( 0, 10, 11 )
  ybase = [ 0.000000, 39.039767, 63.239573, 78.903977, 86.954413, 92.580019, 95.877789, 97.672163, 98.642095, 98.933075, 99.854510 ]
  ybow = [ 0.000000, 0.969932, 3.782735, 9.456838, 15.130941, 23.423860, 29.631426, 40.834142, 48.108632, 52.812803, 60.620757 ]
  ybwa = [ 0.000000, 0.533463, 1.891368, 5.577110, 9.990301, 18.331717, 25.315228, 34.238603, 43.501455, 51.939864, 60.378274 ]
  
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.plot(x, ybase, label='Baseline', marker='s')
  ax.plot(x, ybow, label='Bowtie2', marker='^')
  ax.plot(x, ybwa, label='BWA', marker='o')
  
  ax.set_ylabel('Unmapped %')
  ax.set_xlabel('Error %')
  leg = ax.legend(prop={'size':12}, loc='lower right')
  leg.get_frame().set_alpha(0.8)
  
  fig.savefig('%s/error-unmapped-hiv.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)

def plot_heatmap( out_file, x_field_name, y_field_name, z_field_name, target_file, y_label ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x_fields = set()
  y_fields = set()
  z_values = {}
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x_field = find_parameter( line, x_field_name )
    y_field = find_parameter( line, y_field_name )
    z_field = find_column( line, z_field_name )
    x_fields.add( int(x_field) )
    y_fields.add( int(y_field) )
    z_values[ '%s,%s' % ( x_field, y_field ) ] = z_field

  print z_values
  # make matrix
  m = np.zeros((len(y_fields), len(x_fields)))
  x_pos = 0
  y_pos = 0
  x_field_list = list(x_fields)
  x_field_list.sort()
  y_field_list = list(y_fields)
  y_field_list.sort()
  for x_field in x_field_list:
    y_pos = 0
    for y_field in y_field_list:
      #print "assigning to %i %i" % ( y_pos, x_pos )
      m[y_pos, x_pos] = z_values[ '%s,%s' % ( x_field, y_field ) ]
      y_pos += 1
    x_pos += 1

  # draw graph
  fig = plt.figure()
  #fig.subplots_adjust(right=0.5)
  #ax = fig.add_subplot(111)
  gs = mpl.gridspec.GridSpec(1, 2, height_ratios=(1,), width_ratios=(15, 1))

  ax = plt.subplot(gs[0,0])
  heatmap = ax.pcolor(m, cmap=plt.cm.RdYlGn )
  #fig.subplots_adjust(right=0.9, left=0.4)
  #cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
  cbar_ax = plt.subplot(gs[0,1])
  #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7], frameon=False)
  cb = plt.colorbar(heatmap, cax=cbar_ax)
  #cbar_ax = fig.add_subplot(122)
  #fig.colorbar(im, cax=cbar_ax)
  #cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=plt.cm.RdYlGn, spacing='proportional')#, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
  ax.set_ylabel(y_label)
  ax.set_xlabel('Read length' )

  ax.set_yticks(np.arange(m.shape[0])+0.5, minor=False)
  ax.set_xticks(np.arange(m.shape[1])+0.5, minor=False)
  ax.set_xticklabels(x_field_list)
  ax.set_yticklabels(y_field_list)
  #leg = ax.legend(loc='upper right', prop={'size':12})
  #leg.get_frame().set_alpha(0.8)
  #fig.savefig('%s/ecoli-insertion-vs-alignment.pdf' % REPORT_DIRECTORY, format='pdf', dpi=1000)
  #fig.tight_layout()
  #ax.set_position( ( 0.4, 0.4, 0.3, 0.5 ) )
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file ), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def plot_comparison( out_file, positions, names, x_field, x_name, y_field, y_name, y2_field, y2_name, y3_field, y3_name, target_file, with_precision_recall=True ):
  colors = ('r', 'g', 'b' )
  markers = ('o', '^', 's' )
  y2_marker = '+'
  y3_marker = 'x'
  #out_file = "out/insert_readlength_150108.out"
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y2 = []
  y3 = []
  y = []
  current = 0
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( float( find_parameter( line, x_field ) ) * 100 )
    y.append( float( find_column( line, y_field ) ) )
    y2.append( float( find_column( line, y2_field ) ) )
    y3.append( float( find_column( line, y3_field ) ) )

  print "x", x
  print "y", y
  print "y2", y2
  print "y3", y3
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  for idx in xrange(len(positions)):
    start = positions[idx]
    if idx < len(positions)-1:
      finish = positions[idx+1]
      ax.plot(x[start:finish], y[start:finish], label=names[idx], color=colors[idx], marker=markers[idx] )
      if with_precision_recall:
        ax.plot(x[start:finish], y2[start:finish], color=colors[idx], alpha=0.3, marker=y2_marker )
        ax.plot(x[start:finish], y3[start:finish], color=colors[idx], alpha=0.3, marker=y3_marker )
    else:
      ax.plot(x[start:], y[start:], label=names[idx], color=colors[idx], alpha=0.8, marker=markers[idx] )
      if with_precision_recall:
        ax.plot(x[start:], y2[start:], color=colors[idx], alpha=0.3, marker=y2_marker, label=y2_name )
        ax.plot(x[start:], y3[start:], color=colors[idx], alpha=0.3, marker=y3_marker, label=y3_name )
  #ax.plot(x, y_snp, label='SNV Call', color='b')
  ax.set_ylabel(y_name)
  ax.set_xlabel(x_name)
  #ax.set_title( 'Alignment accuracy with insertion of length 20' )
  #ax.set_xlim(xmax=max(x) * 1.05)
  ax.set_ylim(ymax=max(y) * 1.05)
  #ax.set_xlim(xmax=0.1)
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )


#### experimental
#plot_entropy( 'out/entropy_test_150319.out', 'circoviridae-x2' )
#plot_entropy( 'out/entropy_circ_150319.out', 'circoviridae' )

#plot_depth( 'out/read-depth-deletion-150330-bwa.out', 'circoviridae-100-bwa')
#plot_depth( 'out/read-depth-deletion-150330-bowtie.out', 'circoviridae-100-bowtie')
#plot_depth( 'out/read-depth-deletion-150330-bwa.out', 'circoviridae-100-zoom-bwa', 900, 1200 )
#plot_depth( 'out/read-depth-deletion-150330-bowtie.out', 'circoviridae-100-zoom-bowtie', 900, 1200 )

#plot_depth( 'out/read-depth-deletion-150331-bwa.out', 'circoviridae-10-bwa' )#, 900, 1200 )
#plot_depth( 'out/read-depth-deletion-150331-bwa.out', 'circoviridae-10-zoom-bwa', 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150331-bowtie.out', 'circoviridae-10-bowtie' )#, 900, 1200 )
#plot_depth( 'out/read-depth-deletion-150331-bowtie.out', 'circoviridae-10-zoom-bowtie', 980, 1030 )

# duplication depth
#plot_depth( 'out/read-depth-duplication-150401.out', 'circoviridae-duplication-100-bwa' )#, 980, 1030 )
#plot_depth( 'out/read-depth-duplication-150401-bowtie.out', 'circoviridae-duplication-100-bowtie' )#, 980, 1030 )
#plot_depth( 'out/read-depth-duplication-150401-bowtie.out', 'circoviridae-duplication-100-bowtie-zoom', 1900, 2100 )
#plot_depth( 'out/read-depth-duplication-150401-coverage.out', 'circoviridae-duplication-100-coverage' )#, 980, 1030 )
#plot_depth( 'out/read-depth-duplication-150401-coverage.out', 'circoviridae-duplication-100-coverage-zoom', 900, 1100 )

#plot_depth( 'out/read-depth-deletion-150407-coverage.out', 'circoviridae-deletion-100-coverage-zoom', 900, 1200 )
#plot_depth( 'out/read-depth-50-deletion-150407-coverage.out', 'circoviridae-deletion-50-coverage-zoom', 900, 1200 )
#plot_depth( 'out/read-depth-30-deletion-150407-coverage.out', 'circoviridae-deletion-30-coverage-zoom', 900, 1200 )
#plot_depth( 'out/read-depth-deletion-150407-100x-poisson.out', 'circoviridae-deletion-100x-poisson-zoom', 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150407-100x-poisson.out', 'circoviridae-deletion-100x-poisson', variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150407-15x-poisson.out', 'circoviridae-deletion-15x-poisson', variation_label='Deletion' )

#plot_depth( 'out/read-depth-30-deletion-150407-coverage.out', 'circoviridae-deletion-30-coverage-megazoom', 970, 1010 )
#plot_depth( 'out/read-depth-deletion-150407-coverage.out', 'circoviridae-deletion-100-coverage-megazoom', 960, 1010 )
#plot_depth( 'out/read-depth-deletion-150407-coverage.out', 'circoviridae-deletion-100-coverage-megazoom-suffix', 1091, 1160 )

# breakpoints
#plot_depth( 'out/read-depth-deletion-150408-100x-breakpoints.out', 'circoviridae-deletion-100-coverage-evidence', features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150408-100x-breakpoints.out', 'circoviridae-deletion-100-coverage-evidence-zoom', 900, 1200, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-breakpoints-rl50.out', 'circoviridae-deletion-100-coverage-evidence-zoom-rl50', 950, 1150, features=('Start', 'End'), show_features=False ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector.out', 'circoviridae-deletion-100-coverage-evidence-deletions', variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector.out', 'circoviridae-deletion-100-coverage-evidence-zoom-deletions', 950, 1150, variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-poisson.out', 'circoviridae-deletion-100-coverage-evidence-deletions-poisson', variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-poisson.out', 'circoviridae-deletion-100-coverage-evidence-zoom-deletions-poisson', 950, 1150, variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-10x-deletion-detector-poisson.out', 'circoviridae-deletion-10-coverage-evidence-deletions-poisson', variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-10x-deletion-detector-poisson.out', 'circoviridae-deletion-10-coverage-evidence-zoom-deletions-poisson', 950, 1150, variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-40del.out', 'circoviridae-deletion-100-coverage-evidence-deletions-40del', variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-40del.out', 'circoviridae-deletion-100-coverage-evidence-zoom-deletions-40del', 950, 1150, variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-40del.out', 'circoviridae-deletion-100-coverage-depth-zoom-deletions-40del', 950, 1150, variation_label='Deletion', features=('Evidence',), show_features=False ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-20del.out', 'circoviridae-deletion-100-coverage-depth-zoom-deletions-20del', 950, 1100, variation_label='Deletion', features=('Evidence',), show_features=False ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-30del.out', 'circoviridae-deletion-100-coverage-depth-zoom-deletions-30del', 950, 1100, variation_label='Deletion', features=('Evidence',), show_features=False ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-10del.out', 'circoviridae-deletion-100-coverage-depth-zoom-deletions-10del', 950, 1100, variation_label='Deletion', features=('Evidence',), show_features=False ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-100x-deletion-detector-45del.out', 'circoviridae-deletion-100-coverage-depth-zoom-deletions-45del', 950, 1150, variation_label='Deletion', features=('Evidence',), show_features=False ) #, features=('Start', 'End') )#, 980, 1030 )

#plot_depth( 'out/read-depth-deletion-150409-10x-deletion-detector-poisson-50del.out', 'circoviridae-deletion-100-coverage-evidence-zoom-deletions-poisson-50del', 900, 1150, variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-10x-deletion-detector-poisson-50del.out', 'circoviridae-deletion-100-coverage-evidence-deletions-poisson-50del', variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-10x-deletion-detector-50del.out', 'circoviridae-deletion-100-coverage-evidence-zoom-deletions-50del', 900, 1150, variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )
#plot_depth( 'out/read-depth-deletion-150409-10x-deletion-detector-50del.out', 'circoviridae-deletion-100-coverage-evidence-deletions-50del', variation_label='Deletion', features=('Evidence',), show_features=True ) #, features=('Start', 'End') )#, 980, 1030 )

# novel
#plot_depth( 'out/read-depth-novel-150401-coverage.out', 'circoviridae-novel-100-coverage' )
#plot_depth( 'out/read-depth-novel-150401-coverage.out', 'circoviridae-novel-100-coverage-zoom', 900, 1100 )

#### in report
#plot_depth( 'out/read-depth-deletion-150407-100x.out', 'circoviridae-deletion-100x-zoom', 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150407-15x-poisson.out', 'circoviridae-deletion-15x-poisson-zoom', 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150414-10x-poisson-50del.out', 'circoviridae-deletion-15x-poisson', variation_label='Deletion', features=('Evidence',), show_features=True )#, 900, 1200, variation_label='Deletion' )
plot_depth( 'out/read-depth-deletion-150414-10x-poisson-50del-err05.out', 'circoviridae-deletion-15x-poisson-err05', variation_label='Deletion', features=('Evidence',), show_features=True )#, 900, 1200, variation_label='Deletion' )
plot_depth( 'out/read-depth-deletion-150414-10x-poisson-50del-err05.out', 'circoviridae-deletion-15x-poisson-err05-coverage', variation_label='Deletion', features=('Evidence',), show_features=False )#, 900, 1200, variation_label='Deletion' )

#plot_heatmap( 'out/circoviridae-insertion-heatmap-150309a.out', 'read_length', 'max_insertion_len', 'vcf_f1', 'circoviridae-insertion-heatmap', 'Insertion length' )
#plot_heatmap( 'out/ecoli-insertion-heatmap-150309a.out', 'read_length', 'max_insertion_len', 'vcf_f1', 'ecoli-insertion-heatmap', 'Insertion length' )
#plot_heatmap( 'out/circoviridae-deletion-heatmap-150309a.out', 'read_length', 'max_deletion_len', 'vcf_f1', 'circoviridae-deletion-heatmap', 'Deletion length' )

#mutation_hiv('out/pipeline_batch_mutation_hiv_20140813.out') # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_hiv_snp_f1('out/pipeline_batch_mutation_hiv_20140813.out' ) # mutation-f1-snp-hiv.pdf
#mutation_ecoli_snp_f1( 'out/pipeline_batch_mutation_ecoli_20140813.out' ) # ecoli-mutations-snps-2.pdf
#plot_read_length_vs_alignment_ecoli_detailed()
#plot_coverage_vs_alignment_ecoli_high_errors()
#coverage_ecoli('out/pipeline_batch_coverage_ecoli_20140813.out' ) # ecoli-coverage-snp.pdf
#snp_vs_map_ecoli( 'out/pipeline_batch_read_length_ecoli_20140813_inchworm.out' ) # ecoli-map_vs_snp-400.pdf plot of f1-scores of snps found and mapping accuracy, vs mutation rate, for e-coli
#plot_insertion_vs_readlength( "out/insert_readlength_ecoli_detailed_150304.out", "ecoli-insertion-vs-readlength")
#plot_insertion_vs_readlength( "out/delete_readlength_ecoli_detailed_150304.out", "ecoli-deletion-vs-readlength")
#plot_insertion_vs_variation_distance( 'out/ecoli_variation_distance150305.out', 'min_variation_distance', 'Variation Distance', 'Accuracy vs Variation Distance', 'ecoli-insert-variation-distance' )
#plot_insertion_vs_variation_distance( 'out/ecoli_variation_distance_deletion_150305.out', 'min_variation_distance', 'Variation Distance', 'Accuracy vs Variation Distance', 'ecoli-deletion-variation-distance' )
#plot_insertion_vs_alignment( "out/insert_experiment_141209.out", "Insertion length", "ecoli-insertion-vs-alignment", "max_insertion_len" )
#plot_insertion_vs_alignment( "out/delete_experiment_ecoli_141211.out", "Deletion length", "ecoli-deletion-vs-alignment", "max_deletion_len" )
#plot_insertion_vs_variation_distance( 'out/ecoli-cluster-150309.out', 'min_variation_distance', 'Variation Distance', 'Accuracy vs Variation Distance', 'ecoli-insert-cluster' )
#plot_insertion_vs_variation_distance( 'out/ecoli-cluster-delete-150309.out', 'min_variation_distance', 'Variation Distance', 'Accuracy vs Variation Distance', 'ecoli-delete-cluster' )
#plot_comparison( 'out/insert-mutations-150310.out', (0, 11, 22), ('Novel', 'Homopolymer', 'Duplication'), "error_prob", "Substitution Rate (%)", "vcf_f1", "F1 Accuracy (%)", 'vcf_precision', 'Precision', 'vcf_recall', 'Recall', 'insertion-vs-mutation-150311' )
#plot_comparison( 'out/ecoli-insert-homopolymer-compare-150311.out', (0, 11, 22), ('Novel', 'Homopolymer', 'Duplication'), "error_prob", "Homopolymer Insertion Rate (%)", "vcf_f1", "F1 Accuracy (%)", 'vcf_precision', 'Precision', 'vcf_recall', 'Recall', 'insertion-vs-homopolymer-150311', with_precision_recall=False )
#plot_insertion_vs_variation_distance( 'out/ecoli-insert-homopolymer-150311a.out', 'custom', 'Homopolymer Insertion Length', 'title', 'ecoli-insert-homopolymer-150311', legend="lower left" )

#plot_reference_bias_ecoli_example() # reference-bias-profile.pdf

#### in report - migrated from ipython
#error_unmapped_hiv()
#ecoli_mutations_snps_f1()
#ecoli_mutations_snps_unmapped()
#mutation_f1_hiv()
#error_f1_hiv()


# old - no longer used
# ecoli-mutations-snps-unmapped-2.pdf (not used)
# mutation-f1-hiv.pdf (not used)
# ecoli-mutations-snps-f1-2.pdf (not used)
#plot_read_length_vs_alignment_ecoli()

##### old
#plot_insertion_vs_variation_distance( 'out/insert-homopolymer-150311a.out', 'custom', 'Homopolymer Insertion Length', 'title', 'circoviridae-insert-homopolymer-150311' )
#plot_insertion_vs_variation_distance( 'out/insert-homopolymer-rate-150311.out', 'error_prob', 'Error Rate', 'title', 'circoviridae-homopolymer-error-rate', as_percent=True, legend='lower left' )

#mutation_ecoli_snp_histogram() # ecoli-mutations-snps-2.pdf
#read_length_ecoli() # ecoli-read-length-bowtie-snp.pdf
#read_length_ecoli_map() # ecoli-read-length-bowtie-mapping.pdf
#read_length_hiv() # hiv-read-length-vcf-snps.pdf (lots of noise)

##### experimental
#read_length_ecoli_map_low_mutation() # ecoli-read-length-bowtie-mapping.pdf
#plot_mutation_vs_bias_ecoli()
#plot_reference_bias_ecoli_example_with_errors() # 
#plot_insertion_vs_alignment_circoviridae()
#plot_deletion_vs_alignment_ecoli()
#plot_deletion_vs_alignment_circoviridae()
#plot_insertion_vs_readlength( 'out/insert_readlength_150108.out', 'circoviridae-insertion-vs-readlength' ) #plot_insertion_vs_readlength_circoviridae()
#plot_insertion_vs_variation_distance( 'out/insert_variation_distance_150226.out', 'min_variation_distance', 'Variation Distance', 'Accuracy vs Variation Distance', 'circoviridae-insert-variation-distance' )
#old plot_vcf_parent_vs_child()
#
#plot_vcf_parent_vs_child( 'Plasmodium_falciparum_3d_p_CHR02_recovered.vcf', 'Plasmodium_falciparum_3d7_1q_CHR02_recovered.vcf', 'malaria-3d7-p-vs-3d7-1q-vcf' )
#plot_vcf( 'hiv_efree.vcf', 'hiv_vcf_error_free' )
#plot_vcf( 'hiv_raw_bwa_150303a.vcf', 'hiv_vcf_raw_bwa_150303' )
#plot_vcf_parent_vs_child( 'hiv_raw_bwa_150303a.vcf', 'hiv_efree.vcf', 'hiv_efree_vs_raw', 'Raw Reads', 'Error Free Reads', highlight_differences=True )
#plot_vcf_parent_vs_child( 'hiv_efreeread_pg_150303.vcf', 'hiv_efree.vcf', 'hiv_efree_pg_vs_efree_provied', 'Generated Error Free', 'Provided Error Free', highlight_differences=True )
#plot_vcf_parent_vs_child( 'hiv_raw_bwa_150303a.vcf', 'hiv_raw_stripped_bwa_150303.vcf', 'hiv_raw_nostrip_vs_strip', 'not stripped', 'stripped', highlight_differences=True )
#plot_vcf_parent_vs_child( 'freebayes_hiv_efree_2015-02-18_normed.vcf', 'hiv_efree.vcf', 'hiv_vcf_vs_freebayes_error_free', 'FreeBayes', 'MGSA' )
#plot_vcf_parent_vs_child( 'hiv_efree_bowtie.vcf', 'hiv_efree.vcf', 'hiv_vcf_error_free_bowtie_vs_bwamem', 'Bowtie', 'BWA mem' )
#plot_vcf_parent_vs_child( 'hiv_efree_bwasw.vcf', 'hiv_efree.vcf', 'hiv_vcf_error_free_bwasw_vs_bwamem', 'BWASW', 'BWA mem' )
#plot_vcf_parent_vs_child( 'freebayes_hiv_raw_bowtie_normed_2015-02-25.vcf', 'hiv_raw_bowtie_filtered.vcf', 'hiv_vcf_vs_freebayes_raw', 'FreeBayes', 'MGSA' )
#plot_vcf_parent_vs_child( 'freebayes_hiv_raw_bowtie_10_normed.vcf', 'hiv_raw_bowtie_filtered.vcf', 'hiv_vcf_vs_freebayes_raw_20150225', 'FreeBayes', 'MGSA' )
#plot_vcf_parent_vs_child_chromosomes( 'sample_multi.vcf', 'sample_multi_recovered.vcf', 'sample_multi' )
#plot_vcf_parent_vs_child_chromosomes( 'Plasmodium_falciparum_p_resequenced_recovered.vcf', 'Plasmodium_falciparum_1q_resequenced_recovered.vcf', 'Plasmodium_falciparum_3d7_1q.vcf', 'malaria_vcf_1q_snp', legend_loc='upper right', legend_bbox=(1.,0.9) )
#plot_vcf_parent_vs_child_chromosomes( 'freebayes_3d_p_filtered_normed.vcf', 'freebayes_1q_filtered_normed.vcf', 'Plasmodium_falciparum_3d7_1q.vcf', 'freebayes_malaria_vcf_1q_snp', ymax=120, legend_loc='upper right', legend_bbox=(1.,0.9) )
#plot_vcf_parent_vs_child_chromosomes( 'Plasmodium_falciparum_3d_p_resequenced_recovered.vcf', 'Plasmodium_falciparum_1a_resequenced_recovered.vcf', 'Plasmodium_falciparum_3d7_1a.vcf', 'malaria_vcf_1a_snp', legend_loc='upper right', legend_bbox=(1.,0.9) )
#plot_vcf_parent_vs_child_chromosomes( 'freebayes_3d_p_filtered_normed.vcf', 'freebayes_1a_normed.vcf', 'Plasmodium_falciparum_3d7_1a.vcf', 'freebayes_malaria_vcf_1a_snp', ymax=50, legend_loc='upper right', legend_bbox=(1.,0.9) )
#plot_vcf_parent_vs_child_chromosomes( 'Plasmodium_falciparum_p_ref13_resequenced_recovered.vcf', 'Plasmodium_falciparum_1q_ref13_resequenced_recovered.vcf', 'Plasmodium_falciparum_1q_ref13.vcf', 'malaria_vcf_1q_ref13_snp', legend_loc='upper right', legend_bbox=(1.,0.9) )
#plot_vcf_parent_vs_child_chromosomes( 'freebayes_p_ref13_normed.vcf', 'freebayes_1q_ref13_normed.vcf', 'Plasmodium_falciparum_1q_ref13.vcf', 'freebayes_malaria_vcf_1q_ref13_snp', ymax=10, legend_loc='upper right', legend_bbox=(1.,0.9) )
#plot_vcf_parent_vs_child_chromosomes( 'Plasmodium_falciparum_p_ref9_resequenced_recovered.vcf', 'Plasmodium_falciparum_1q_ref9_resequenced_recovered.vcf', 'Plasmodium_falciparum_3d7_1q.vcf', 'malaria_vcf_1q_snp_ref9', legend_loc='upper right', legend_bbox=(1.,0.9) )#, chromosome_map=MALARIA_CHROMOSOME_MAP )
# ----- mappability -----
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
#plot_mappability( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88.out', 'HIV NC001722 with 1 SNP', 'hiv001722-bwasw-snp' )
#plot_mappability_hist( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88.out', 'HIV NC001722 with 1 SNP', 'hiv001722-bwasw-snp' )
#plot_bias( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88.out', 'HIV NC001722 with 1 SNP', 'hiv001722-bwasw-snp' )
#plot_bias_hist( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88.out', 'HIV NC001722 with 1 SNP', 'hiv001722-bwasw-snp' )
#plot_all_mappability( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88_ins2-2.out', 'HIV NC001722 with 2x2 inserts', 'hiv001722-bwasw-snp-ins22' )
#plot_all_mappability( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88_snp10.out', 'HIV NC001722 with 10 SNPs', 'hiv001722-bwasw-snp-snp10' )
#plot_all_mappability( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88_snp5.out', 'HIV NC001722 with 5 SNPs', 'hiv001722-bwasw-snp-snp5' )
#plot_all_mappability( 'out/mappability_hiv-NC_001722_bwasw_snp_rl88_ins5.out', 'HIV NC001722 with 5 inserts', 'hiv001722-bwasw-snp-ins5' )
#
#plot_all_mappability( 'out/mappability_hiv-NC_001802_bwasw_snp_rl88_ins2-2.out', 'HIV NC001802 with 2x2 inserts', 'hiv001802-bwasw-snp-ins22' )
#plot_all_mappability( 'out/mappability_hiv-NC_001802_bwasw_snp_rl88_snp10.out', 'HIV NC001802 with 10 SNPs', 'hiv001802-bwasw-snp-snp10' )
#plot_all_mappability( 'out/mappability_hiv-NC_001802_bwasw_snp_rl88_snp5.out', 'HIV NC001802 with 5 SNPs', 'hiv001802-bwasw-snp-snp5' )
#plot_all_mappability( 'out/mappability_hiv-NC_001802_bwasw_snp_rl88_ins5.out', 'HIV NC001802 with 5 inserts', 'hiv001802-bwasw-snp-ins5' )
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
#plot_reference_bias( 'out/ecoli_reference_bias_vs_repeats_150212.out', 'ecoli' )
