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

import scipy.stats
import numpy as np

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

columns = [ 'unmapped', 'incorrect', 'read_precision', 'read_recall', 'read_f1', 'vcf_tp', 'vcf_fp', 'vcf_fn', 'vcf_precision', 'vcf_recall', 'vcf_f1', 'vcf_covered', 'vcf_bucket_tp', 'vcf_bucket_fp', 'vcf_bucket_fn', 'reference_bias', 'error_bias', 'unmapped_variations', 'total_variations', 'mean_reference', 'mean_error' ]
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
  
def get_column_values( columns, fh, default=None ):
  while True:
    line = fh.readline()
    if line == "":
      print "get_value", column, "returning None"
      return None
    if line.startswith( '#' ):
      continue
    values = []
    for column in columns:
      value = find_column( line, column )
      if value is not None:
        if default is not None and default[0] == float(value):
          values.append( default[1] )
        else:
          values.append( value )
    return values

def get_value( column, fh, default=None ):
  while True:
    line = fh.readline()
    if line == "":
      print "get_value", column, "returning None"
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

from matplotlib.ticker import FuncFormatter
def to_millions(y, position):
  s = '%.1fM' % ( y / 1000000 )
  return s

def to_thousands(y, position):
  s = '%.1fk' % ( y / 1000 )
  return s

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

##########################################################################
# run job
def mutation_unmapped( out_file, target_file, has_base=True, x_dist=[], x_bwa=19, x_bowtie=20, include_incorrect=True ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_base = []
  y_bow = []
  y_bwa = []
  for x in xrange(len(x_dist)):
    if include_incorrect:
      values = get_column_values( ['unmapped', 'incorrect'], fh )
      y_bow.append( 100 - float( values[0] ) - float( values[1] ))
      values = get_column_values( ['unmapped', 'incorrect'], fh )
      y_bwa.append( 100 - float( values[0] ) - float( values[1] ))
    else:
      if has_base:
        y_base.append( 100 - float( get_value( 'unmapped', fh ) ) )
      y_bow.append( 100 - float( get_value( 'unmapped', fh ) ) )
      y_bwa.append( 100 - float( get_value( 'unmapped', fh ) ) )
  
  # draw graph
  #x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  if has_base:
    ax.plot(x_dist, y_base, label='Bowtie2 exhaustive', color='b', marker='s')
  ax.plot(x_dist, y_bow, label='Bowtie2 local', color='g', marker='^')
  ax.plot(x_dist, y_bwa, label='BWA', color='r', marker='o')

  ax.axvline(x=x_bwa, color='black', alpha=0.2)#, label='BWA k')
  ax.axvline(x=x_bowtie, color='black', alpha=0.2)#, label='Bowtie2 k')

  ax.set_ylabel('Mapped %')
  ax.set_xlabel('Distance between mutations')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  ax.text(x_bwa - 1.4, 3,'BWA\nseed',rotation=0, color='gray')
  ax.text(x_bowtie + 0.2, 3,'Bowtie2\nseed',rotation=0, color='gray')

  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def mutation_seed( out_file, target_file, x_dist=[] ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_bwa = []
  for x in xrange(len(x_dist)):
      y_bwa.append( 100 - float( get_value( 'unmapped', fh ) ) )
  
  # draw graph
  #x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x_dist, y_bwa, label='BWA', color='r', marker='o')

  ax.set_ylabel('Mapped %')
  ax.set_xlabel('Seed Length')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def mutation_seed_2( out_file, target_file, x_dist=[] ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_bwa = []
  y_bwa_2 = []
  for x in xrange(len(x_dist)):
      y_bwa.append( 100 - float( get_value( 'unmapped', fh ) ) )
      y_bwa_2.append( 100 - float( get_value( 'unmapped', fh ) ) )
  
  # draw graph
  #x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x_dist, y_bwa, label='Seed length 20', color='r', marker='o')
  ax.plot(x_dist, y_bwa_2, label='Seed length 8', color='b', marker='s')

  ax.set_ylabel('Mapped %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )

def mutation_covered_variation( out_file, target_file, x_dist ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_bwa = []
  for x in xrange(len(x_dist)):
      y_bwa.append( float( get_value( 'vcf_covered', fh ) ) )
  
  print y_bwa
  # draw graph
  #x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x_dist, y_bwa, label='BWA', color='r', marker='o')

  ax.set_ylabel('Covered %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
  bio.log_stderr( 'extracting values from %s: done' % out_file )



def mutations_unmapped_random( target_file, y_predict, y, mapper_name ):
  # e-coli vs unmapped
  # from pipeline_batch_mutation_ecoli-140805.out
  #x = np.linspace( 2, 16, 8 )
  x = np.linspace( 1, 10, 10 )
  print x
  print y_predict
  print y
  #ybase = [39.362935, 63.848565, 78.354124, 86.817800, 92.319459, 95.545825, 97.332601, 98.467877, 99.101603, 99.497479]
  #ybow = [ 100 - z for z in (0.982743, 4.099772, 9.130041, 15.653223, 23.007691, 31.118041, 39.108282, 46.678516, 53.948746, 60.810065) ]
  #ybow = [ 100 - z for z in () ]
  #0.000000 0.000000 0.000000 0.000000 1.072125 3.216374 5.360624 7.407407 12.670565 16.666667 20.760234

  # calculated using pipeline_batch.py batch/pipeline_batch_mutation_hiv.cfg out/pipeline_batch_mutation_hiv_150830.out
  #ybwa = [ 100 - z for z in (0.430347, 2.065667, 5.465788, 10.584389, 17.202473, 25.366899, 34.136181, 43.066400, 51.827064, 60.169230) ]
  #ybow_predict = [ 100 - 100. * z for z in ( 0.000307166262346, 0.00594893896531, 0.0277733871657, 0.0731128114326, 0.141654810427, 0.227446007244, 0.322403441667, 0.418934270243, 0.511219584631, 0.595502841535, ) ] # k=22
  #ybow_predict = [ 100 - 100. * z for z in ( 0.000307166262346, 0.00594893896531, 0.0277733871657, 0.0731128114326, 0.141654810427, 0.227446007244, 0.322403441667, 0.418934270243, 0.511219584631, 0.595502841535, ) ] # k=22
  #ybow_predict = [ 100 - 100. * z for z in ( 0.000200201348861, 0.00405744832828, 0.0197607458106, 0.0540955478139, 0.108650306444, 0.180287987205, 0.263309530458, 0.351499923261, 0.43942862253, 0.523021264566, ) ] # k=20

  # calculated using expected_unmapped
  #ybwa_predict = [ 100 - 100. * z for z in ( 2.75910996713e-05, 0.0010491635715, 0.00719589559447, 0.0246796454746, 0.0582705179043, 0.109216203632, 0.175332840799, 0.252272694499, 0.334955622291, 0.418639358939, ) ] # k=19
  #ybwa_predict = [ 100 - 100. * z for z in (0, 0.000158723357544, 0.0032913141409, 0.0163771142189, 0.0457385540111, 0.0935852486123, 0.157969252687, 0.234361386071, 0.317364640847, 0.4019352168, 0.484023796648, ) ] # k=19
  #ybwa = [ 100 - z for z in (14.814815, 14.814815, 14.717349, 14.814815, 15.399610, 16.666667, 16.861598, 18.518519, 22.027290, 25.146199, 34.600390,) ]
  #ybwa = [ 100 - z for z in (0.000000, 0.000000, 0.000000, 0.000000, 0.255102, 1.020408, 1.530612, 5.612245, 16.836735, 16.071429, 14.030612,) ]
  


  #bowfit = np.poly1d(np.polyfit(x, ybow, 3))(x)
  #bwafit = np.poly1d(np.polyfit(x, ybwa, 1))(x)
  #yave = [ (z[0]+z[1])/2. for z in zip(ybow,ybwa) ]
  #avefit = np.poly1d(np.polyfit(x, yave, 1))(x)
  
  fig = plt.figure()
  ax = fig.add_subplot(111)

  #ax.plot(x, ybase, label='Baseline', color='b', marker='s')
  #ax.plot(x, ybow_predict, label='Bowtie2 Prediction', color='b', marker='s')
  ax.plot(x, y_predict, label='%s Prediction' % mapper_name, color='b', marker='s')
  #ax.plot(x, ybow, label='Bowtie2 Actual', color='g', marker='^')
  ax.plot(x, y, label='%s Actual' % mapper_name, color='r', marker='o')
  #ax.plot(x, avefit, '--', label='Linear Fit', color='#a0a0a0')
  #plt.plot(x, linear, label='Linear', color='#cccccc')
  #ax.set_ylim(ymin=0.0)
  
  ax.set_ylabel('Mapped %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
 
def mutation_hiv( out_file, target_file, has_base=True ):
  #out_file = "out/pipeline_batch_mutation_hiv_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_mutation_hiv.cfg %s" % out_file )
  
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  y_base = []
  y_bow = []
  y_bwa = []
  for x in xrange(0, 11):
    if has_base:
      y_base.append( get_value( 'unmapped', fh ) )
    y_bow.append( get_value( 'unmapped', fh ) )
    y_bwa.append( get_value( 'unmapped', fh ) )
  
  # draw graph
  x = np.linspace( 0, 10, 11 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  if has_base:
    ax.plot(x, y_base, label='Baseline', color='b', marker='s')
  ax.plot(x, y_bow, label='Bowtie2', color='g', marker='^')
  ax.plot(x, y_bwa, label='BWA', color='r', marker='o')
  ax.set_ylabel('Unmapped %')
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)
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

def coverage_ecoli(out_file, target, measure, y_axis, depths=(5, 10, 20)):
  #out_file = "out/pipeline_batch_coverage_ecoli_%s.out" % datetime.datetime.now().strftime("%Y%m%d")
  #run( "python pipeline_batch.py batch/pipeline_batch_coverage_ecoli.cfg %s" % out_file )
  colors = ('b', 'g', 'r', 'orange')
  markers = ('s', '^', 'o', '+')
  
  bio.log_stderr( 'extracting values from %s...' % out_file )
  fh = open( out_file, 'r' )
  results = [ [] for _ in depths ] 
  print results
  for x in xrange(0, 11):
    for i, depth in enumerate(depths):
      value = get_value( measure, fh, default=(0,100) )
      results[i].append( value )

  # draw graph
  x = np.linspace( 1, 10, 10 )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  for i, depth in enumerate(depths):
    print x, results[i][1:]
    ax.plot(x, results[i][1:], label='%ix coverage' % depth, color=colors[i], marker=markers[i])
  #ax.plot(x, b10[1:], label='10x coverage', color='g', marker='^')
  #ax.plot(x, b20[1:], label='20x coverage', color='r', marker='o')
  ax.set_ylabel(y_axis)
  ax.set_xlabel('Mutation %')
  leg = ax.legend(prop={'size':12})
  leg = ax.legend(loc='lower left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target), format='pdf', dpi=1000)
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

def plot_read_length( out_file, target, show_covered=True, start_from=0 ):
  bio.log_stderr( 'extracting values from %s...' % out_file )
  x = []
  y_mapped = []
  y_covered = []
  # python find_repeats.py 50 < ../../data/hiv.fasta
  # 50, 100, 150, 200, 300, 400, 500, 750, 1000
  repeated = [ 1.352169, 1.140607, 1.019120, 0.928980, 0.792520, 0.684024, 0.592569, 0.388655, 0.226665 ]
  #hiv = [ 7.780674, 
  fh = open( out_file, 'r' )
  for line in fh:
    if line.startswith( '#' ):
      continue
    x.append( find_parameter( line, 'read_length' ) )
    y_mapped.append( 100 - float( find_column( line, 'unmapped' ) ) )
    #y_incorrect.append( find_column( line, 'incorrect' ) )
    y_covered.append( find_column( line, 'vcf_covered' ) )

  print "x", x
  #print "y", y
  # draw graph
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x[start_from:], y_mapped[start_from:], label='Mapped %', color='r', marker='o')
  if show_covered:
    ax.plot(x[start_from:], y_covered[start_from:], label='Covered %', color='b', marker='s')
  ax.plot(x[start_from:], [ 100 - x for x in repeated ][start_from:], label='Unrepeated Content', color='g', marker='+')
  ax.set_ylabel('%')
  ax.set_xlabel('Read Length')
  #ax.set_ylim(ymin=90)
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target), format='pdf', dpi=1000)
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

  
def plot_vcf( fn_list, out_fn, labels ):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  colors = ('r', 'g', 'b')

  xmax = 1
  for idx, fn in enumerate(fn_list):
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
      xmax = max(xmax, snp['pos'] )
  
    ax.scatter(x, y, edgecolor="none", color=colors[idx], label='%s (%i)' % ( labels[idx], len(vcf.snp_list) ) ) #, log=True)

  ax.set_ylabel('Quality')
  ax.set_xlabel('Position')
  ax.set_xlim( xmin=0, xmax=xmax )
  ax.set_ylim( ymin=0 )

  leg = ax.legend(loc='upper left', prop={'size':10})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

def plot_vcf_parent_vs_child( parent_fn, child_fn, out_fn, parent_label='parent', child_label='child', highlight_differences=False, legend='upper right' ):
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
  leg = ax.legend(loc=legend, prop={'size':8})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, out_fn ), format='pdf', dpi=1000)

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

def plot_depth_predictor( src, short_name, my_depths=0, my_breakpoints_lines=1, legend='upper left', my_max_depth=None, window=0, x_label='Depth', normalize=True ):
  '''
    just does 1 depth set and 1 breakpoint
  '''
  lines = open( src ).readlines()
  depths = []
  breakpoints = set()
  idx = 0
  for line in lines:
    if line.startswith( '#' ) or line.strip() == '':
      continue
    if idx == my_depths:
      depths.extend( [ int(round(float(x))) for x in line.split(',') ] )
      print "added", idx
    if idx == my_breakpoints_lines:
      breakpoints.update( [ int(round(float(x))) for x in line.split(',') ] )
      print "added breakpoint ", idx, len(breakpoints)
    idx += 1

  max_depth = max(depths)
  print "building hist with", max_depth, "depths"
  if my_max_depth is not None:
    max_depth = my_max_depth
  x = range(0, max_depth + 1 )
  hist_p = [0] * ( max_depth + 1 )
  hist_n = [0] * ( max_depth + 1 )
  hist_c = [0] * ( max_depth + 1 )
  for position in xrange(0, len(depths)):
    if position > 100 and position < len(depths) - 100:
      #window_depth = sum( depths[ position - window: position + window + 1 ] )
      #ave_depth = window_depth / ( 2 * window + 1 )
      window_depth = min( depths[ position - window: position + window + 1 ] )
      ave_depth = window_depth
      if ave_depth < max_depth:
        if position in breakpoints:
          hist_p[ave_depth] += 1
        #elif position + 1 in breakpoints or position - 1 in breakpoints or position + 2 in breakpoints or position - 2 in breakpoints or position + 3 in breakpoints or position - 3 in breakpoints or position + 4 in breakpoints or position - 4 in breakpoints or position + 5 in breakpoints or position - 5 in breakpoints:
        #  hist_c[depth] += 1
        else:
          hist_n[ave_depth] += 1

  print "hist_p", sum(hist_p), "hist_n", sum(hist_n)
  # normalize
  if normalize:
    hist_n = [ 1. * v / sum(hist_n) for v in hist_n ]
    hist_p = [ 1. * v / ( sum(hist_p) + sum(hist_c) ) for v in hist_p ]
#  hist_c = [ 1. * v / ( sum(hist_p) + sum(hist_c) ) for v in hist_c ]

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.bar(x, hist_n, width=1, color='b', label='No breakpoint', edgecolor='b', alpha=0.5)
  ax.bar(x, hist_p, width=1, color='r', label='Breakpoint', edgecolor='r', alpha=0.5)
#  ax.bar(x, hist_c, width=1, color='orange', label='Near breakpoint', edgecolor='orange', alpha=0.5, bottom=hist_p)

  ax.set_ylabel('Normalized Frequency')
  ax.set_xlabel(x_label)

  leg = ax.legend( loc=legend, prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/predict-%s.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)

  # roc from hists
  return bio.generate_roc( hist_n, hist_p )

def plot_rocs( rocs, labels, legend='lower right', short_name='test', xmax=1., colors=('b', 'r', 'g'), markers=('s', 'o', '+') ):
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_ylabel('True Positive Rate')
  ax.set_xlabel('False Positive Rate')

  for idx, roc in enumerate(rocs):
    ax.plot( roc['fpr'], roc['tpr'], color=colors[idx], label=labels[idx], marker=markers[idx], markersize=4 )
    ax.fill_between( roc['fpr'], roc['tpr'], 0, color=colors[idx], alpha=0.1)

  plt.xlim(xmax=xmax)
  leg = ax.legend( loc=legend, prop={'size':12})
  leg.get_frame().set_alpha(0.8)

  fig.savefig('%s/roc-%s.pdf' % (REPORT_DIRECTORY, short_name), format='pdf', dpi=1000)
  

def plot_depth( src, short_name, range_start=0, range_stop=-1, variation_label='Variation', features=(), show_features=False, show_breakpoints=True, labels=('Coverage',), my_depths=(0,), my_breakpoints=None, my_breakpoints_lines=None, my_max=0, legend='upper left', colors=('b', 'g', 'magenta', 'red', 'black'), my_horizontal=None ):
  '''
    plots either depth or evidence
    @my_breakpoints: list of lists
  '''
  lines = open( src ).readlines()
  lists = []
  breakpoint_lines = []
  idx = 0
  for line in lines:
    if line.startswith( '#' ) or line.strip() == '':
      continue
    if idx in my_depths:
      lists.append( [ float(x) for x in line.split(',') ] )
      print "added", idx
    if idx in my_breakpoints_lines:
      breakpoint_lines.append( [ float(x) for x in line.split(',') ] )
      print "added breakpoint ", idx
    idx += 1

  # lists[0] = depth, lists[1] = breakpoints
  last_depth = lists[-1]
  fig = plt.figure()
  ax = fig.add_subplot(111)
  if range_stop < 0:
    range_stop = len(last_depth)
  first = True
  breakpoints = None
  if show_breakpoints:
    if my_breakpoints:
      breakpoints = my_breakpoints
    elif my_breakpoints_lines is not None:
      breakpoints = breakpoint_lines[0] # only supports one breakpoint list
    else:
      breakpoints = lists[1 + len(features)]
  if breakpoints is not None:
    for breakpoint in breakpoints:
      if range_start <= breakpoint < range_stop:
        if first:
          first = False
          ax.axvline(x=breakpoint, color='r', alpha=0.1, label=variation_label)
        else:
          ax.axvline(x=breakpoint, color='r', alpha=0.1 )

  if not show_features: # do the depth
    max_depths = my_max
    for idx, depths in enumerate(lists):
      max_depths = max(depths[range_start:range_stop] + [ max_depths ] )
      x = range(range_start, min(range_stop, len(depths)))
      y = depths[range_start:min(range_stop, len(depths))]
      print "plotting", idx, "x", len(x), "y", len(y)
      ax.plot(x, y, label=labels[idx], color=colors[idx] ) #, s=2)
      ax.fill_between( x, y, 0, color=colors[idx], alpha=0.1)
    plt.ylim(ymax=max_depths+0.5, ymin=-0.5)

  # additional features
  else:
    for idx in xrange(0, len(features)):
      ax.plot(range(range_start, range_stop), lists[1 + idx][range_start:range_stop], label=features[idx] ) #, s=2)

  if my_horizontal is not None:
    xlims = ax.get_xlim()
    x = range(range_start, int(xlims[1]))
    ax.axhline( y=my_horizontal[0], color='r', alpha=0.3, label='Detector')
    ax.axhline( y=my_horizontal[1], color='r', alpha=0.3)
    ax.fill_between( x, my_horizontal[0], my_horizontal[1], color='r', alpha=0.3, label='Detector')

  ax.set_ylabel('Coverage')
  ax.set_xlabel('Position')

  leg = ax.legend( loc=legend, prop={'size':12})
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

def plot_multi_zero_depth( src, fn, labels, colors, show_index=None, max_depth=None, min_length=0 ):
  lines = open( src ).readlines()
  lists = []
  runs = []
  actual_max_depth = 0
  for line in lines:
    if line.startswith( '#' ) or line.strip() == '':
      continue
    lists.append( [ float(x) for x in line.split(',') ] ) # depths
    trimmed = bio.trim_edges( lists[-1] ) # without edges
    all_runs = bio.all_runs( trimmed, 0 )
    runs.append( all_runs ) # dictionary of run length and count of runs
    longest = bio.longest_run( trimmed, 0 )
    total_runs = sum( [ all_runs[k] for k in  all_runs ] )
    threshold_runs = sum( [ all_runs[k] for k in all_runs if k <= min_length ] )
    print "longest run: %i @ %i. total runs: %i below threshold %i" % ( longest[0], longest[1], total_runs, threshold_runs )
    if len(lists) == 1:
      base_total = total_runs
      base_threshold = threshold_runs
    else:
      threshold_inc = threshold_runs - base_threshold
      total_inc = total_runs - base_total
      print "fp above %i %.4f below %i %.4f" % ( total_inc - threshold_inc, 1. * ( total_inc - threshold_inc ) / total_inc, threshold_inc, 1 - 1. * ( total_inc - threshold_inc ) / total_inc )

  fig = plt.figure()
  hists = []
  hist_labels = []
  for idx, run in enumerate( runs ):
    if show_index is None or idx in show_index:
      ax = fig.add_subplot(111)
      x = []
      y = []
      z = [] # expanded out i.e. { 50: 3 } -> [ 50, 50, 50 ]
      total = 0
      for k, v in run.iteritems():
        if k > min_length:
          x.append( k )
          y.append( v )
          z.extend( [k] * v )
          total += v
      #ax.scatter(x, y, label=label, color=colors[idx], alpha=0.7, edgecolors='none')
      #bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
      #bins = range( 0, 1000, 10 )
      #n, bins, patches = ax.hist(z, bins=50, label=label, histtype='step', alpha=0.5 )#, normed=1, histtype='stepfilled')
      if len(z) > 0:
        hists.append( z )
        label='%s (%i)' % ( labels[idx], total )
        print '%s: %s' % ( label, sorted( run.keys() ) )
        hist_labels.append( label )

  if len(hists) > 0:
    n, bins, patches = ax.hist(hists, label=hist_labels, bins=30 )#, normed=1, histtype='stepfilled')
  #ax.set_xscale('log')
  #ax.set_yscale('log')
  ax.set_xlabel( 'Deletion length' )
  ax.set_ylabel( 'Count' )

  leg = ax.legend(loc='upper right', prop={'size':12})
  #leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, fn ), format='pdf', dpi=1000)


def plot_multi_depth_hist( src, fn, labels, colors, max_depth=None ):
  lines = open( src ).readlines()
  lists = []
  means = []
  actual_max_depth = 0
  for line in lines:
    if line.startswith( '#' ) or line.strip() == '':
      continue
    lists.append( [ float(x) for x in line.split(',') ] )
    means.append( sum( lists[-1] ) / len( lists[-1] ) )
    actual_max_depth = max( actual_max_depth, max( lists[-1] ) )
    trimmed = bio.trim_edges( lists[-1] )
    longest = bio.longest_run( trimmed, 0 )
    print 'longest zero run: %i at %i-%i of %i. distribution: %s' % ( longest[0], longest[1], longest[1] + longest[0], len(lists[-1]), bio.all_runs( trimmed ) )

  if max_depth is None:
    max_depth = actual_max_depth

  x = range( 0, int(max_depth) + 1 ) 
  max_height = 0
  fig = plt.figure()
  for idx, depth in enumerate( lists ):
    y = bio.bucket( depth, x)
    max_height = max( max_height, max(y))
    ax = fig.add_subplot(111)
    ax.plot(x, y, label=labels[idx], color=colors[idx])
    ax.fill_between( x, y, 0, color=colors[idx], alpha=0.1)
    ax.axvline(x=means[idx], color=colors[idx], alpha=1)

  if max_height > 1e6:
    ax.get_yaxis().set_major_formatter(FuncFormatter(to_millions))
  #ax1.set_xlabel('my label ' + '$10^{{{0:d}}}$'.format(scale_pow))

  ax.set_ylabel('Frequency')
  ax.set_xlabel('Depth of Coverage')

  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, fn ), format='pdf', dpi=1000)

def plot_depth_hist_from_bedtools( srcs, fn, labels=None, max_depth=None, with_poisson=False, with_nbinom=False, nbinom_var=3 ):
  colors = ('b', 'g', 'r' )
  nbinom_colors = ('green', 'cyan')
  depths_dict = []
  means = []
  true_max_depth = 0
  max_count = 0
  pmfs = []
  for idx, src in enumerate(srcs):
    total_depths = 0
    total_count = 0
    depths_dict.append(collections.defaultdict(int))
    for line in open( src ):
      if line.startswith( 'genome' ):
        fields = line.split()
        depth = int(fields[1])
        count = int(fields[2])
        total_count += count
        total_depths += count * depth
        depths_dict[idx][depth] += count
        true_max_depth = max( true_max_depth, depth )
        max_count = max( max_count, count )
    mean = 1. * total_depths / total_count
    means.append( mean )
    pmfs.append( total_count )

  if max_depth is None:
    max_depth = true_max_depth

  x = range( 0, max_depth + 1 )
  # y = bio.bucket( depths, x )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  for i in xrange(len(depths_dict)):
    y = [ depths_dict[i][depth] for depth in x ]
    ax.plot(x, y, label=labels[i], color=colors[i] )#, marker='s')
    ax.fill_between( x, y, 0, color=colors[i], alpha=0.1)
    if with_poisson or with_nbinom:
      ax.axvline(x=means[i], color='black', alpha=0.3)
    else:
      ax.axvline(x=means[i], color=colors[i], alpha=1)
    if i == 0:
      with open( 'out.txt', 'w' ) as fh:
        for yi in y:
          fh.write( '%i\n' % int(yi) )

  if with_poisson:
    dist = scipy.stats.poisson(means[0])
    y = pmfs[0] * dist.pmf(x)
    ax.plot( x, y, label='Poisson', color='r' )
    ax.fill_between( x, y, 0, color='r', alpha=0.1)
        
  if with_nbinom:
    for idx, v in enumerate(nbinom_var):
      var = means[0] * v # np.std( [ depths_dict[0][depth] for depth in x ] )
      #var = depths_sd * depths_sd
      #n = 200 #0.4 #( var * var - means[0] ) / ( means[0] * means[0] )
      #p = 0.4 #n / ( means[0] + n )
      n = ( means[0] * means[0] ) / ( var - means[0] )
      p = n / ( means[0] + n )
      dist = scipy.stats.nbinom(n, p) #means[0], 3)
      y = pmfs[0] * dist.pmf(x)
      ax.plot( x, y, label='Negative binomial %ix' % v, color=nbinom_colors[idx] )
      ax.fill_between( x, y, 0, color=nbinom_colors[idx], alpha=0.1)
        
  if max_count > 1e5:
    ax.get_yaxis().set_major_formatter(FuncFormatter(to_millions))
  elif max_count > 1e4:
    ax.get_yaxis().set_major_formatter(FuncFormatter(to_thousands))
  print "max_count", max_count

  ax.set_ylabel('Frequency')
  ax.set_xlabel('Depth of Coverage')
  leg = ax.legend(loc='lower right', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, fn ), format='pdf', dpi=1000)


def plot_depth_hist( src, fn ):
  lines = open( src ).readlines()
  lists = []
  means = []
  for line in lines:
    if line.startswith( '#' ) or line.strip() == '':
      continue
    lists.append( [ float(x) for x in line.split(',') ] )
    means.append( sum( lists[-1] ) / len( lists[-1] ) )

  # lists[0] = depth, lists[1] = breakpoints
  depths = lists[0]
  x = range( 0, int(max(depths)) + 1 ) 
  y = bio.bucket( depths, x )
  print "first 10 coverages:", y[:10]
  print "total genome:", len(depths)
  print "total coverage:", sum(depths)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, y, label='Baseline', color='b' )#, marker='s')
  ax.axvline(x=means[0], color='b', alpha=1)
  ax.set_ylabel('Frequency')
  ax.set_xlabel('Depth of Coverage')
  fig.savefig('%s/%s.pdf' % ( REPORT_DIRECTORY, fn ), format='pdf', dpi=1000)

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

def plot_tandem( bwa, bowtie, bowtie_ecoli, target ):
  x = np.linspace( 1, len(bwa), len(bwa) )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  btm = [a * 100./44 for a in bowtie]
  #btme = [a * 100./44 for a in bowtie_ecoli]
  bwam = [a * 100./60 for a in bwa]
  ax.plot(x, btm, label='Bowtie2', marker='^')
  #ax.plot(x, btme, label='Bowtie2 Ecoli', marker='x')
  ax.plot(x, bwam, label='BWA', marker='o')
  
  ax.set_ylabel('Mapping Quality (%)')
  ax.set_xlabel('Repeat Period')
  leg = ax.legend(prop={'size':12}, loc='lower right')
  leg.get_frame().set_alpha(0.8)
  
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target), format='pdf', dpi=1000)


def plot_mutation_unmapped( bwa, bowtie, target ):
  x = np.linspace( 0, len(bwa)-1, len(bwa) )
  fig = plt.figure()
  ax = fig.add_subplot(111)
  btm = [100-a for a in bowtie]
  bwam = [100-a for a in bwa]
  print x, btm, bwam
  ax.plot(x, btm, label='Bowtie2', marker='^')
  ax.plot(x, bwam, label='BWA', marker='o')
  
  ax.set_ylabel('Unmapped %')
  ax.set_xlabel('Read Distance')
  leg = ax.legend(prop={'size':12}, loc='lower right')
  leg.get_frame().set_alpha(0.8)
  
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target), format='pdf', dpi=1000)


def plot_heatmap( out_file, x_field_name, y_field_name, z_field_name, target_file, y_label, x_label='Read Length', z_field=False ):
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
    if z_field:
      z_field = find_parameter( line, z_field_name )
    else:
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
  ax.set_xlabel(x_label)

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

def plot_indels( target_file ):
  x = (100, 200, 300, 400, 500, 750, 1000)
  bwa_ins = (28,  59,  94,  122, 163, 267, 365 )
  bowtie_ins = ( 23,   48,  73, 98, 123, 182, 245 )
  bwa_del = ( 45,  94,  130, 132, 156, 211, 273 )
  bowtie_del = ( 30 , 50 , 70 , 90 , 107, 159, 195)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(x, bwa_ins, label='BWA Insertion Length', color='r', alpha=0.8, marker='s' )
  ax.plot(x, bowtie_ins, label='Bowtie2 Insertion Length', color='g', alpha=0.8, marker='o' )
  ax.plot(x, bwa_del, label='BWA Deletion Length', color='b', alpha=0.8, marker='^' )
  ax.plot(x, bowtie_del, label='Bowtie2 Deletion Length', color='orange', alpha=0.8, marker='+' )
  leg = ax.legend(loc='upper left', prop={'size':12})
  leg.get_frame().set_alpha(0.8)
  fig.savefig('%s/%s.pdf' % (REPORT_DIRECTORY, target_file), format='pdf', dpi=1000)

#### experimental
#plot_vcf_parent_vs_child( 'e-coli-bw2952.vcf', 'e-coli-mg1655.vcf', 'e-coli-bw2952-mg1655', 'K12', 'MG1655', highlight_differences=True, legend='lower left' )
#plot_vcf_parent_vs_child( 'e-coli-iai39.vcf', 'e-coli-mg1655.vcf', 'e-coli-iai39-mg1655', 'IAI39', 'MG1655', highlight_differences=True, legend='lower left' )
#plot_vcf( ('e-coli-iai39.vcf','e-coli-bw2952.vcf', 'e-coli-mg1655.vcf'), 'ecoli-mapping-vcf-comparison', labels=('IAI39', 'K12', 'MG1655') )

# novel insertion rd 10
#plot_depth( 'out/read-depth-novel-150430-coverage.out', 'cicroviridae-depth-unique-insertion-200-150430', 941, 1059, show_breakpoints=True, labels=('50bp', '100bp', '500bp', '1000bp'), my_depths=(0, 2, 4, 6), legend='lower right', my_breakpoints_lines=(1,) )
#plot_depth( 'out/read-depth-novel-150430-coverage-poisson.out', 'cicroviridae-depth-unique-insertion-200-150430-poisson', 941, 1059, show_breakpoints=False, labels=('50bp', '100bp', '500bp', '1000bp'), my_depths=(0, 2, 4, 6), legend='lower right', my_breakpoints_lines=(1,) )
#plot_depth( 'out/read-depth-novel-150430-coverage-poisson.out', 'cicroviridae-depth-unique-insertion-200-150430-poisson-all', show_breakpoints=True, labels=('100bp',), my_depths=(2,), my_breakpoints_lines=(3,), legend='lower right' )
#
# novel insertion rd 100
#plot_depth( 'out/read-depth-novel-150503-coverage-short.out', 'cicroviridae-depth-unique-insertion-rd100-150503-short', 941, 1059, show_breakpoints=True, labels=('10bp', '20bp', '30bp', '100bp'), my_depths=(0, 2, 4, 8), legend='lower right', my_breakpoints_lines=(1,) )
#plot_depth( 'out/read-depth-novel-150430-coverage-rd100.out', 'cicroviridae-depth-unique-insertion-rd100-150430', 941, 1059, show_breakpoints=True, labels=('50bp', '100bp', '500bp', '1000bp'), my_depths=(0, 2, 4, 6), legend='lower right', my_breakpoints_lines=(1,) )
#plot_depth( 'out/read-depth-novel-150430-coverage-poisson-rd100.out', 'cicroviridae-depth-unique-insertion-rd100-150430-poisson', 941, 1059, show_breakpoints=False, labels=('50bp', '100bp', '500bp', '1000bp'), my_depths=(0, 2, 4, 6), legend='lower right', my_breakpoints_lines=(1,) )
#plot_depth( 'out/read-depth-novel-150430-coverage-poisson-rd100.out', 'cicroviridae-depth-unique-insertion-rd100-150430-poisson-all', show_breakpoints=True, labels=('100bp',), my_depths=(2,), my_breakpoints_lines=(3,), legend='lower left', my_horizontal=(50, 70,) )
# don't do this plot_depth( 'out/read-depth-novel-150430-coverage-poisson-ecoli.out', 'ecoli-depth-unique-insertion-rd100-150430-poisson-all', show_breakpoints=True, labels=('100bp',), my_depths=(0,), my_breakpoints_lines=(1,), legend='lower right', my_horizontal=(50, 70,) )
#
#bwa_roc = plot_depth_predictor( 'out/read-depth-novel-150430-coverage-poisson-ecoli.out', 'ecoli-depth-unique-insertion-rd100-150430-poisson-window5', my_depths=0, my_breakpoints_lines=1, legend='upper right', my_max_depth=150, window=8, x_label='Min depth (16bp window)')
#bowtie_roc = plot_depth_predictor( 'out/read-depth-novel-150506-coverage-poisson-ecoli-bowtie.out', 'ecoli-depth-unique-insertion-rd100-150430-poisson-window5-bowtie', my_depths=0, my_breakpoints_lines=1, legend='upper right', my_max_depth=150, window=8, x_label='Min depth (16bp window)')

#lr_fpr = [ 0.078876679, 0.050061050, 0.039316239, 0.033455433, 0.029304029, 0.022466422, 0.021245421, 0.018559219, 0.016849817, 0.015628816, 0.015140415, 0.013675214, 0.012454212, 0.011233211, 0.010744811, 0.009279609, 0.008302808, 0.006593407, 0.005372405, ]
#lr_tpr = [ 0.9801105, 0.9779006, 0.9767956, 0.9745856, 0.9701657, 0.9679558, 0.9624309, 0.9502762, 0.9469613, 0.9425414, 0.9303867, 0.9138122, 0.9027624, 0.8773481, 0.8541436, 0.8176796, 0.7690608, 0.6950276, 0.5624309, ]
#lr_roc = { 'fpr': lr_fpr, 'tpr': lr_tpr }
#plot_rocs( ( bwa_roc, bowtie_roc, lr_roc ), labels=('BWA', 'Bowtie', 'Regression'), short_name='bwa-vs-bowtie-vs-regression', xmax=0.5 )

###
#bwa_roc = plot_depth_predictor( 'out/read-depth-novel-150430-coverage-poisson-ecoli.out', 'ecoli-depth-unique-insertion-rd100-150430-poisson-window5-norm', my_depths=0, my_breakpoints_lines=1, legend='upper right', my_max_depth=150, window=8, x_label='Min depth (16bp window)', normalize=False)
#bowtie_roc = plot_depth_predictor( 'out/read-depth-novel-150506-coverage-poisson-ecoli-bowtie.out', 'ecoli-depth-unique-insertion-rd100-150430-poisson-window5-bowtie-norm', my_depths=0, my_breakpoints_lines=1, legend='upper right', my_max_depth=150, window=8, x_label='Min depth (16bp window)', normalize=False)
#plot_rocs( ( bwa_roc, bowtie_roc ), labels=('BWA', 'Bowtie'), short_name='bwa-vs-bowtie-norm', xmax=0.5 )
#
#plot_depth_predictor( 'out/read-depth-novel-150430-coverage-poisson-ecoli.out', 'ecoli-depth-unique-insertion-rd100-150430-poisson', my_depths=0, my_breakpoints_lines=1, legend='upper right', my_max_depth=200, window=0)
#plot_depth_predictor( 'out/read-depth-novel-150430-coverage-poisson-rd100.out', 'circoviridae-depth-unique-insertion-rd100-150430-poisson', my_depths=0, my_breakpoints_lines=1, legend='upper right', my_max_depth=200, window=8)

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

#### poisson coverage deletions bwa ####
# no deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats.out', 'e-coli-depth-vs-repeats-150424-empty', show_breakpoints=False, range_start=299000, range_stop=301000, skip=0 )
# unique deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats.out', 'e-coli-depth-vs-repeats-150424-unique', show_breakpoints=False, range_start=299000, range_stop=301000, skip=1, my_breakpoints=( range(299999, 300200), range(3300198, 3300190) ), variation_label='Deletion' )
# no deletion across repeat
#plot_depth( 'out/ecoli-150424-depth-vs-repeats.out', 'e-coli-depth-vs-repeats-150424-empty-repeat', show_breakpoints=False, range_start=257000, range_stop=259000, skip=0 )
# deletion across repeat
#plot_depth( 'out/ecoli-150424-depth-vs-repeats.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat', show_breakpoints=False, range_start=257000, range_stop=259000, skip=3, my_breakpoints=( range(257999, 258200) ), variation_label='Deletion' )
# corresponding repeat with no deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats.out', 'e-coli-depth-vs-repeats-150424-empty-repeat-other-loc', show_breakpoints=False, range_start=19000, range_stop=21000, skip=0 )
# corresponding repeat with deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-other-loc', show_breakpoints=False, range_start=19000, range_stop=21000, skip=3, my_breakpoints=( range(19887, 20088) ), variation_label='Pseudo Deletion' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out

#plot_depth_hist_from_bedtools( ('out/Plasmodium_falciparum_1q_orig_dedup.hist',), 'malaria-depth-hist-150422-dedup', max_depth=500, labels=('Observed',), with_poisson=True )#, 900, 1200, variation_label='Deletion' )
#### uniform coverage deletions bwa ####
# no deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-empty-uniform', show_breakpoints=False, range_start=299000, range_stop=301000, skip=0 )
# unique deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-unique-uniform', show_breakpoints=False, range_start=299000, range_stop=301000, skip=1, my_breakpoints=( range(299999, 300200), range(3300198, 3300190) ), variation_label='Deletion' )
# no deletion across repeat
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-empty-repeat-uniform', show_breakpoints=False, range_start=257000, range_stop=259000, skip=0 )
# deletion across repeat
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-uniform', show_breakpoints=False, range_start=257000, range_stop=259000, skip=3, my_breakpoints=( range(257999, 258200) ), variation_label='Deletion' )
# corresponding repeat with no deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-empty-repeat-other-loc-uniform', show_breakpoints=False, range_start=19000, range_stop=21000, skip=0, my_max=50 )
# corresponding repeat with deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-other-loc-uniform', show_breakpoints=False, range_start=19000, range_stop=21000, skip=3, my_breakpoints=( range(19887, 20088) ), variation_label='Pseudo Deletion', my_max=50 ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out
#### uniform coverage deletions bowtie ####
# no deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-empty-uniform-bowtie', show_breakpoints=False, range_start=299000, range_stop=301000, skip=0 )
# unique deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-unique-uniform-bowtie', show_breakpoints=False, range_start=299000, range_stop=301000, skip=1, my_breakpoints=( range(299999, 300200), range(3300198, 3300190) ), variation_label='Deletion' )
# no deletion across repeat
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-empty-repeat-uniform-bowtie', show_breakpoints=False, range_start=257000, range_stop=259000, skip=0 )
# deletion across repeat
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-uniform-bowtie', show_breakpoints=False, range_start=257000, range_stop=259000, skip=3, my_breakpoints=( range(257999, 258200) ), variation_label='Deletion' ) # bps via sed -n 11p
# corresponding repeat with no deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-empty-repeat-other-loc-uniform-bowtie', show_breakpoints=False, range_start=19000, range_stop=21000, skip=0 )
# corresponding repeat with deletion
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-other-loc-uniform-bowtie', show_breakpoints=False, range_start=19000, range_stop=21000, skip=3, my_breakpoints=( range(19887, 20088) ), variation_label='Pseudo Deletion' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out
#plot_depth_hist( 'out/malaria-depth-150422.out', 'malaria-depth-hist-150422' )#, 900, 1200, variation_label='Deletion' )
#plot_depth_hist_from_bedtools( ['out/Plasmodium_falciparum_1q_orig.hist', 'out/Plasmodium_falciparum_1q_orig_dedup.hist'], 'malaria-orig-depth-hist-150422-bedtools', labels=['Original', 'Deduped'], max_depth=400 )#, 900, 1200, variation_label='Deletion' )
#plot_depth_hist( 'out/malaria-depth-dedup-150422.out', 'malaria-depth-hist-150422-dedup-fixed')#, max_depth=500 )#, 900, 1200, variation_label='Deletion' )

#plot_multi_zero_depth( 'out/malaria-depth-150422.out', 'malaria-depth-150422.out', ('Real Readset',), colors=('blue',), show_index=set( (0,) ), max_depth=19, min_length=50 )
#plot_multi_zero_depth( 'out/malaria-depth-dedup-150422.out', 'malaria-depth-dedup-150422.out', ('Real Readset',), colors=('blue',), show_index=set( (0,) ), max_depth=19, min_length=50 )

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

#plot_depth_hist( 'out/read-depth-deletion-150414-10x-poisson-50del-err05.out', 'circoviridae-deletion-15x-poisson-err05-depth-hist' )#, 900, 1200, variation_label='Deletion' )
#plot_multi_depth_hist( 'out/read-depth-deletion-150416-depth-vs-errors.out', 'circoviridae-depth-vs-errors', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta') )#, 900, 1200, variation_label='Deletion' )
#plot_multi_depth_hist( 'out/ecoli-150416-depth-vs-errors.out', 'ecoli-depth-vs-errors', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta'), max_depth=19 )#, 900, 1200, variation_label='Deletion' )
#plot_multi_depth_hist( 'out/ecoli-150416-depth-vs-errors-bowtie.out', 'ecoli-depth-vs-errors-bowtie', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta'), max_depth=19 )#, 900, 1200, variation_label='Deletion' )
#plot_multi_depth_hist( 'out/ecoli-150416-depth-vs-errors-hom5.out', 'ecoli-depth-vs-errors-hom5', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta'), max_depth=19 )#, 900, 1200, variation_label='Deletion' )
#plot_multi_depth_hist( 'out/ecoli-150416-depth-vs-errors-hom5-bowtie.out', 'ecoli-depth-vs-errors-hom5-bowtie', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta'), max_depth=19 )#, 900, 1200, variation_label='Deletion' )
#plot_multi_depth_hist( 'out/ecoli-150416-depth-vs-errors-snp.out', 'ecoli-depth-vs-errors-snp', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta'), max_depth=19 )#, 900, 1200, variation_label='Deletion' )

#plot_multi_zero_depth( 'out/read-depth-deletion-150416-depth-vs-errors.out', 'circoviridae-depth-vs-errors-zeros', ('1%', '2%', '4%', '8%', '16%'), colors=('blue', 'green', 'red', 'cyan', 'magenta') )#, 900, 1200, variation_label='Deletion' )

#### in report
# reference repeats etc
#plot_heatmap( 'out/result-cand-100samples.fix', 'x', 'y', 'z', 'reference-read-candidate-distance', 'Candidate Distance', 'Read Distance', z_field=True )
#plot_mutation_unmapped( bwa = [ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 98.0, 95.0, 82.0, 79.0, 75.0, ], bowtie = [ 100.0, 100.0, 100.0, 100.0, 100.0, 98.0, 99.0, 92.0, 78.0, 76.0, 58.0, 53.0, ], target = 'reference-unmapped-distance' )
# this is from result-tandem.out, which is circoviridae with 100 samples and from the 7 candidate column
#plot_tandem( bwa = [ 1.7, 2.0, 6.9, 7.2, 8.7, 8.3, 10.5, 13.7, 15.2, 15.7, 18.5, 18.9, 19.9, 21.6, 24.0, 25.3, 25.3, 27.0, 28.0, 28.3, 29.5, 30.6, 31.7, 33.9, 34.5, 36.0, 36.4, 38.0, 39.3, 40.5, 41.4, 42.5, 43.9, 45.3, 45.7, 48.0, 48.1, 49.8, 50.5, ], bowtie = [ 6.5, 9.3, 14.3, 13.4, 15.4, 16.4, 13.3, 14.7, 14.9, 15.7, 16.3, 16.9, 16.4, 16.7, 16.4, 21.8, 20.6, 21.1, 20.7, 22.9, 20.0, 19.5, 20.5, 22.9, 21.2, 23.0, 23.4, 23.4, 19.9, 20.1, 20.5, 21.2, 20.3, 20.3, 21.3, 21.5, 21.2, 20.9, 22.9, ], bowtie_ecoli = [8.5, 5.9, 13.4, 15.5, 14.4, 17.8, 13.4, 13.8, 14.7, 14.8, 14.9, 15.9, 18.8, 17.2, 17.9, 22.6, 19.7, 20.4, 20.7, 27.9, 19.0, 19.4, 20.2, 24.0, 24.9, 26.1, 23.7, 24.4, 20.3, 18.0, 20.3, 24.1, 20.1, 19.8, 20.7, 21.1, 21.0, 21.0, 23.0], target='reference-tandem-quality' )

# this is from result-tandem.out, which is circoviridae with 100 samples
#plot_tandem( bwa = [ ], bowtie = [ ], target='reference-tandem-quality-ecoli' )

# combined unique with empty
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-deleted-unique-other-loc-uniform-combined', show_breakpoints=False, range_start=299000, range_stop=301000, my_depths=(0,1), my_breakpoints=( range(299999, 300200) ), variation_label='Deletion location', my_max=25, labels=('Without deletion', 'With deletion'), legend='upper right' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out
# combined deletion with empty
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-other-loc-uniform-combined', show_breakpoints=False, range_start=257000, range_stop=259000, my_depths=(0,3), my_breakpoints=( range(257999, 258200) ), variation_label='Deletion location', my_max=25, labels=('Without deletion', 'With deletion'), legend='upper right' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out
# combined pseudo with empty
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform.out', 'e-coli-depth-vs-repeats-150424-pseudo-deleted-repeat-other-loc-uniform-combined', show_breakpoints=False, range_start=19000, range_stop=21000, my_depths=(0,3), my_breakpoints=( range(19887, 20088) ), variation_label='Pseudo deletion location', my_max=70, labels=('Without deletion', 'With deletion'), legend='upper right' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out

#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-deleted-unique-other-loc-uniform-bowtie-combined', show_breakpoints=False, range_start=299000, range_stop=301000, my_depths=(0,1), my_breakpoints=( range(299999, 300200) ), variation_label='Deletion location', my_max=25, labels=('Without deletion', 'With deletion'), legend='upper right' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-deleted-repeat-other-loc-uniform-bowtie-combined', show_breakpoints=False, range_start=257000, range_stop=259000, my_depths=(0,3), my_breakpoints=( range(257999, 258200) ), variation_label='Deletion location', my_max=25, labels=('Without deletion', 'With deletion'), legend='upper right' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out
#plot_depth( 'out/ecoli-150424-depth-vs-repeats-uniform-bowtie.out', 'e-coli-depth-vs-repeats-150424-deleted-pseudo-other-loc-uniform-bowtie-combined', show_breakpoints=False, range_start=19000, range_stop=21000, my_depths=(0,3), my_breakpoints=( range(19887, 20088) ), variation_label='Pseudo deletion location', my_max=25, labels=('Without deletion', 'With deletion'), legend='upper right' ) # bps via sed -n 8p out/ecoli-150424-depth-vs-repeats.out

#plot_depth_hist_from_bedtools( ('out/Plasmodium_falciparum_1q_orig_dedup.hist',), 'malaria-depth-hist-150422-dedup-poisson', max_depth=500, labels=('Observed',), with_poisson=True, with_nbinom=True, nbinom_var=(3, 6) )#, 900, 1200, variation_label='Deletion' )
#plot_depth_hist_from_bedtools( ('out/ERR588535-bwa-basic.hist',), 'e-coli-mg1655-ERR588535-bwa-basic', max_depth=500, labels=('Observed',), with_poisson=True, with_nbinom=True, nbinom_var=(3, 6) )#, 900, 1200, variation_label='Deletion' )

#plot_depth( 'out/read-depth-deletion-150407-100x.out', 'circoviridae-deletion-100x-zoom', 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150407-15x-poisson.out', 'circoviridae-deletion-15x-poisson-zoom', 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150414-10x-poisson-50del.out', 'circoviridae-deletion-15x-poisson', variation_label='Deletion', features=('Evidence',), show_features=True )#, 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150414-10x-poisson-50del-err05.out', 'circoviridae-deletion-15x-poisson-err05', variation_label='Deletion', features=('Evidence',), show_features=True )#, 900, 1200, variation_label='Deletion' )
#plot_depth( 'out/read-depth-deletion-150414-10x-poisson-50del-err05.out', 'circoviridae-deletion-15x-poisson-err05-coverage', variation_label='Deletion', features=('Evidence',), show_features=False )#, 900, 1200, variation_label='Deletion' )

#plot_heatmap( 'out/circoviridae-insertion-heatmap-150309a.out', 'read_length', 'max_insertion_len', 'vcf_f1', 'circoviridae-insertion-heatmap', 'Insertion length' )
#plot_heatmap( 'out/ecoli-insertion-heatmap-150309a.out', 'read_length', 'max_insertion_len', 'vcf_f1', 'ecoli-insertion-heatmap', 'Insertion length' )
#plot_heatmap( 'out/circoviridae-deletion-heatmap-150309a.out', 'read_length', 'max_deletion_len', 'vcf_f1', 'circoviridae-deletion-heatmap', 'Deletion length' )

#mutation_hiv('out/pipeline_batch_mutation_hiv_20140813.out', 'mutation-hiv') # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_sensitivity('out/pipeline_batch_mutation_hiv_uniform_1.out', 'mutation-hiv-sensitivity-1', False, [35,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15]) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_unmapped('out/pipeline_batch_mutation_hiv_uniform_nomapq.out', 'mutation-hiv-sensitivity-nomapq', False, [35,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15]) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_unmapped('out/pipeline_batch_mutation_hiv_poisson.out', 'mutation-hiv-unmapped-poisson', False, [35,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15]) # mutation-hiv.pdf (% unmapped as mutation rate increases)

#mutation_unmapped('out/pipeline_batch_mutation_ecoli_uniform.out', 'mutation-ecoli-unmapped-nomapq', False, [35,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15]) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_unmapped('out/pipeline_batch_mutation_ecoli_uniform_rl88.out', 'mutation-ecoli-unmapped-nomapq-rl88', False, [35,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15]) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_unmapped('out/pipeline_batch_mutation_hiv_uniform_mapq_rl44.out', 'mutation-hiv-44-mapq', False, [35,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15], include_incorrect=False) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_unmapped('out/pipeline_batch_mutation_hiv_uniform_mapq_rl100.out', 'mutation-hiv-mapq-rl100', False, [30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15], include_incorrect=False) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_unmapped('out/pipeline_batch_mutation_hiv_uniform_150830.out', 'mutation-hiv-uniform-150830', True, [30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15], include_incorrect=False) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutation_seed( 'out/pipeline_batch_mutation_hiv_random_seed_150830.out', 'mutation-hiv-random-seed-150830', [19, 18, 17, 15, 14, 13, 12, 11, 10, 9])
#mutation_seed_2( 'out/pipeline_batch_mutation_ecoli_bwa_seed_compare.out', 'mutation-ecoli-random-seed-compare-150830', [2,4,6,8,10,12,14,16])
#mutation_covered_variation( 'out/pipeline_coverage_ecoli_poisson.out', 'mutation-covered-variation-ecoli-poisson', [2, 4, 6, 8, 10, 12, 14, 16] )
#plot_read_length( 'out/pipeline_batch_many_read_lengths_ecoli_150830.out', 'read-length-ecoli-150830', start_from=1 )
#plot_read_length( 'out/pipeline_batch_many_read_lengths_ecoli_150831.out', 'read-length-ecoli-150831', start_from=1 )
#plot_read_length( 'out/pipeline_batch_many_read_lengths_ecoli_error_free_150830.out', 'read-length-ecoli-error-free-150830', show_covered=False )
#mutations_unmapped_random( 'ecoli-mutations-unmapped-random-bowtie2.pdf' ) # in report
#mutations_unmapped_random( 'hiv-mutations-unmapped-random-bwa.pdf' ) # in report
#mutation_unmapped('out/pipeline_batch_mutation_hiv_uniform_150830.out', 'mutation-hiv-mapq-rl100', True, [30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15], include_incorrect=False) # mutation-hiv.pdf (% unmapped as mutation rate increases)
#mutations_unmapped_random( 'ecoli-mutations-unmapped-random-bwa.pdf', \
#  [ 100 * z for z in ( 1.0, 1.0, 1.0, 0.998, 0.997, 0.983, 0.967, 0.927, 0.876, 0.826,) ], \
#  [ 100 - z for z in ( 0.000431, 0.010772, 0.087255, 0.376814, 1.109974, 2.603872, 4.867555, 8.536803, 13.033793, 18.754443,) ], \
#  'BWA' )

#mutation_hiv_snp_f1('out/pipeline_batch_mutation_hiv_20140813.out' ) # mutation-f1-snp-hiv.pdf
#mutation_hiv_snp_f1('out/pipeline_batch_mutation_hiv_20140813.out' ) # mutation-f1-snp-hiv.pdf
#mutation_ecoli_snp_f1( 'out/pipeline_batch_mutation_ecoli_20140813.out' ) # ecoli-mutations-snps-2.pdf
#plot_read_length_vs_alignment_ecoli_detailed()
#plot_coverage_vs_alignment_ecoli_high_errors()
#coverage_ecoli('out/pipeline_batch_coverage_ecoli_20140813.out', 'ecoli-coverage-snp-sensitivity', 'vcf_recall', 'Sensitivity' ) # ecoli-coverage-snp.pdf
#coverage_ecoli('out/pipeline_batch_coverage_ecoli_150905.out', 'ecoli-coverage-snp-sensitivity-poisson-150906', 'vcf_recall', 'Sensitivity', depths=(2,5,50,100) ) # ecoli-coverage-snp.pdf
#coverage_ecoli('out/pipeline_batch_coverage_ecoli_uniform_150906.out', 'ecoli-coverage-snp-sensitivity-uniform-150906', 'vcf_recall', 'Sensitivity', depths=(2,5,50,100) ) # ecoli-coverage-snp.pdf
#coverage_ecoli('out/pipeline_batch_coverage_ecoli_20140813.out', 'ecoli-coverage-snp-specificity', 'vcf_precision', 'Precision' ) # ecoli-coverage-snp.pdf
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

#plot_multi_zero_depth( 'out/ecoli-150416-depth-vs-errors-snp.out', 'ecoli-depth-vs-errors-snp-zeros', ('No errors', '1% error', '2% error', '4% error', '8% error', '16% error'), colors=('black', 'blue', 'green', 'red', 'cyan', 'magenta'), show_index=set( (0, 3, 4) ), max_depth=19, min_length=50 )#, 900, 1200, variation_label='Deletion' )
#plot_multi_zero_depth( 'out/ecoli-150419-depth-vs-errors-hom5.out', 'ecoli-depth-vs-errors-hom5-zeros', ('No errors', '1% error', '2% error', '4% error', '8% error', '16% error'), colors=('black', 'blue', 'green', 'red', 'cyan', 'magenta'), show_index=set( (0, 2, 3) ), max_depth=19, min_length=50 )#, 900, 1200, variation_label='Deletion' )
#plot_multi_zero_depth( 'out/ecoli-150419-depth-vs-errors-hom5-bowtie.out', 'ecoli-depth-vs-errors-hom5-zeros-bowtie', ('No errors', '1% error', '2% error', '4% error', '8% error', '16% error'), colors=('black', 'blue', 'green', 'red', 'cyan', 'magenta'), show_index=set( (0, 2, 3) ), max_depth=19, min_length=50 )#, 900, 1200, variation_label='Deletion' )

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
#plot_vcf( ('hiv_efree.vcf',), 'hiv_vcf_error_free', labels=('Error Free',) ) # broken
#plot_vcf( ('hiv_raw_bwa_150303a.vcf',), 'hiv_vcf_raw_bwa_150303', labels=('Raw',) ) # broken
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

plot_indels( 'indel-max' )
