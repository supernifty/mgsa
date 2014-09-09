
import sys

import matplotlib.pyplot as plt

import helpers

def plot_series():
  xs, ys = helpers.series_from_pipeline_batch( fh=open(sys.argv[1]), x='insert_prob', y='incorrect', bias_report=False )
  xs = [ x * 100 for x in xs ]
  plt.plot( xs, ys, label='BWA' )
  plt.ylabel('Incorrect %')
  plt.xlabel('Mutation %')
  plt.legend()
  plt.show()

def plot_error_bias():
  #xs, ys = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True )
  xs, ys = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='error_bias', bias_report=True )
  xs = [ x * 100 for x in xs ]
  print xs, ys
  # weird crap
  ys = [ y + 1e-7 for y in ys ]
  width = 100. / len(ys) * 0.9
  plt.bar( left=xs, height=ys, label='Bias', width=width )
  plt.ylabel('SNPs')
  plt.xlabel('Error Bias (%)')
  #plt.xlabel('Reference Bias (%)')
  
  plt.legend()
  plt.show()

def plot_reference_bias():
  xs, ys = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True )
  xs = [ x * 100 for x in xs ]
  print xs, ys
  # weird crap
  ys = [ y + 1e-7 for y in ys ]
  width = 100. / len(ys) * 0.9
  plt.bar( left=xs, height=ys, label='Bias', width=width )
  plt.ylabel('SNPs')
  plt.xlabel('Reference Bias (%)')
  
  plt.legend()
  plt.show()

def plot_bias():
  xs, yrs = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True )
  _, yes = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='error_bias', bias_report=True )
  xs = [ x * 100 for x in xs ]
  print xs, yrs, yes
  # weird crap
  yrs = [ y * 100 + 1e-7 for y in yrs ]
  yes = [ y * 100 + 1e-7 for y in yes ]
  width = 100. / len(yrs) * 0.4
  #plt.bar( left=xs, height=yrs, label='Reference Bias', width=width, color='b', log=True )
  #plt.bar( left=[x+width for x in xs], height=yes, label='Error Bias', width=width, color='y', log=True )
  plt.bar( left=[x-width/2 for x in xs][1:], height=yrs[1:], label='Reference Bias', width=width, color='b', log=False, align='center')
  plt.bar( left=[x+width/2 for x in xs][1:], height=yes[1:], label='Error Bias', width=width, color='g', log=False, align='center')
  plt.ylabel('SNVs (%)')
  plt.xlabel('Bias (%)')

  # determine averages
  expected_ref = 0.
  for i in xrange(len(xs)):
    expected_ref += xs[i] * yrs[i] / 100.
  expected_err = 0.
  for i in xrange(len(xs)):
    expected_err += xs[i] * yes[i] / 100.
  plt.text( 108, 2.5, "Average reference bias: %.1f%%\nAverage error bias: %.1f%%" % ( expected_ref, expected_err ), ha='right' )
  plt.xlim( xmin=5, xmax=110 )
  plt.legend()
  plt.show()

def plot_bias_compare():
  xs, yrs = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True )
  _, yes = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True, item=2 )
  xs = [ x * 100 for x in xs ]
  print xs, yrs, yes
  # weird crap
  yrs = [ y * 100 + 1e-7 for y in yrs ]
  yes = [ y * 100 + 1e-7 for y in yes ]
  width = 100. / len(yrs) * 0.35
  #plt.bar( left=xs, height=yrs, label='Reference Bias', width=width, color='b', log=True )
  #plt.bar( left=[x+width for x in xs], height=yes, label='Error Bias', width=width, color='y', log=True )
  plt.bar( left=[x-width/2 for x in xs][1:], height=yrs[1:], label='Read length 50bp', width=width, color='b', log=False, align='center')
  plt.bar( left=[x+width/2 for x in xs][1:], height=yes[1:], label='Read length 100bp', width=width, color='g', log=False, align='center')
  plt.ylabel('SNVs (%)')
  plt.xlabel('Bias (%)')

  # determine averages
  expected_ref = 0.
  for i in xrange(len(xs)):
    expected_ref += xs[i] * yrs[i] / 100.
  expected_err = 0.
  for i in xrange(len(xs)):
    expected_err += xs[i] * yes[i] / 100.
  plt.text( 108, 2.5, "Average bias (50bp): %.1f%%\nAverage bias (100bp): %.1f%%" % ( expected_ref, expected_err ), ha='right' )
  plt.xlim( xmin=5, xmax=110 )
  plt.legend()
  plt.show()

if __name__ == '__main__':
  #plot_series()
  #plot_error_bias()
  #plot_reference_bias()
  #plot_bias()
  plot_bias_compare()
