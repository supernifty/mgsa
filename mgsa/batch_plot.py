
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



if __name__ == '__main__':
  #plot_series()
  plot_error_bias()
  plot_reference_bias()
