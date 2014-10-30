
import sys

import matplotlib.pyplot as plt
import numpy as np

import helpers

def plot_series():
  xs, ys = helpers.series_from_pipeline_batch( fh=open(sys.argv[1]), x='mult_snp_prob', y='mean_reference', bias_report=True )
  xs = [ x * 100 for x in xs ]
  ys = [ y * 100 for y in ys ]
  plt.plot( xs, ys, label='Bias' )
  plt.ylabel('Reference Bias (%)')
  plt.xlabel('Reference Mutation')
  #xs, ys = helpers.series_from_pipeline_batch( fh=open(sys.argv[1]), x='mult_snp_prob', y='unmapped', bias_report=True )
  #xs = [ x * 100 for x in xs ]
  #ys = [ y * 100 for y in ys ]
  #plt.plot( xs, ys, label='Unmapped' )
  #plt.ylabel('Unmapped (%)')
  #plt.xlabel('Reference Mutation')
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

def plot_reference_bias(include_zero=False):
  xs, ys = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True )
  xs = [ x * 100 for x in xs ]
  print xs, ys

  if include_zero:
    plt.xlim( xmin=-5, xmax=110 )
  else:
    plt.xlim( xmin=5, xmax=110 )

  # weird crap
  ys = [ y * 100 + 1e-7 for y in ys ]
  width = 100. / len(ys) * 0.75
  #plt.bar( left=xs, height=ys, label='Bias', width=width )
  if include_zero:
    plt.bar( left=[x for x in xs], height=ys, label='SNPs affected (%)', width=width, color='b', log=False, align='center')
  else:
    plt.bar( left=[x for x in xs][1:], height=ys[1:], label='SNPs affected (%)', width=width, color='b', log=False, align='center')
  plt.ylabel('SNPs (%)')
  plt.xlabel('Reference Bias (%)')
  
  plt.legend()
  plt.show()

def plot_bias( include_zero=False, include_unmapped=False ):
  xs, yrs = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True, item=1 )
  _, yes = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='error_bias', bias_report=True, item=1 )
  xs = [ x * 100 for x in xs ]

  # determine averages
  expected_err = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='mean_error', bias_report=True, item=1 ) * 100.
  expected_ref = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='mean_reference', bias_report=True, item=1 ) * 100.
  #expected_ref = 0.
  #for i in xrange(len(xs)):
  #  expected_ref += xs[i] * yrs[i]
  #expected_err = 0.
  #for i in xrange(len(xs)):
  #  expected_err += xs[i] * yes[i]
  plt.text( 108, 2.5, "Average reference bias: %.2f%%\nAverage error bias: %.2f%%" % ( expected_ref, expected_err ), ha='right' )

  if include_zero:
    plt.xlim( xmin=-5, xmax=110 )
  else:
    plt.xlim( xmin=5, xmax=110 )

  if include_unmapped:
    plt.xlim( xmax=115 )
    xs.append( 110. ) # fix later
    unmapped = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='unmapped_variations', bias_report=True )
    total = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='total_variations', bias_report=True )
    yrs.append( 1.0 * unmapped / total )
    yes.append( 0 ) # dummy
    plt.xticks( xs[1:], [ '%i' % int(x) for x in xs[1:-1] ] + [ 'Unmapped' ], rotation=-45 )
  print xs, yrs, yes

  # weird crap
  yrs = [ y * 100 + 1e-7 for y in yrs ]
  yes = [ y * 100 + 1e-7 for y in yes ]
  width = 100. / len(yrs) * 0.4
  #plt.bar( left=xs, height=yrs, label='Reference Bias', width=width, color='b', log=True )
  #plt.bar( left=[x+width for x in xs], height=yes, label='Error Bias', width=width, color='y', log=True )
  if include_zero:
    plt.bar( left=[x-width/2 for x in xs], height=yrs, label='Reference Bias', width=width, color='b', log=False, align='center')
    plt.bar( left=[x+width/2 for x in xs], height=yes, label='Error Bias', width=width, color='g', log=False, align='center')
  else:
    plt.bar( left=[x-width/2 for x in xs][1:], height=yrs[1:], label='Reference Bias', width=width, color='b', log=False, align='center')
    plt.bar( left=[x+width/2 for x in xs][1:], height=yes[1:], label='Error Bias', width=width, color='g', log=False, align='center')
  plt.ylabel('SNVs (%)')
  plt.xlabel('Bias (%)')

  plt.legend()
  plt.show()

def plot_bias_compare( include_unmapped=False, column_offset=0 ):
  xs, yrs = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True, column_offset=column_offset )
  _, yes = helpers.series_from_pipeline_result( fh=open(sys.argv[1]), y='reference_bias', bias_report=True, item=2, column_offset=column_offset )
  xs = [ x * 100 for x in xs ]
  plt.xlim( xmin=5, xmax=110 )

  if include_unmapped:
    plt.xlim( xmax=115 )
    xs.append( 110. ) # fix later
    unmapped = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='unmapped_variations', bias_report=True, column_offset=column_offset )
    total = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='total_variations', bias_report=True, column_offset=column_offset )
    yrs.append( 1.0 * unmapped / total )

    unmapped = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='unmapped_variations', bias_report=True, item=2, column_offset=column_offset )
    total = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='total_variations', bias_report=True, item=2, column_offset=column_offset )
    yes.append( 1.0 * unmapped / total )

    plt.xticks( xs[1:], [ '%i' % int(x) for x in xs[1:-1] ] + [ 'Unmapped' ], rotation=-45 )

  print xs, yrs, yes
  # weird crap
  yrs = [ y * 100 + 1e-7 for y in yrs ]
  yes = [ y * 100 + 1e-7 for y in yes ]
  width = 100. / len(yrs) * 0.35
  #plt.bar( left=xs, height=yrs, label='Reference Bias', width=width, color='b', log=True )
  #plt.bar( left=[x+width for x in xs], height=yes, label='Error Bias', width=width, color='y', log=True )
  plt.bar( left=[x-width/2 for x in xs][1:], height=yrs[1:], label='50 bp', width=width, color='b', log=False, align='center')
  plt.bar( left=[x+width/2 for x in xs][1:], height=yes[1:], label='100 bp', width=width, color='g', log=False, align='center')
  plt.ylabel('SNVs (%)')
  plt.xlabel('Bias (%)')

  # determine averages
  expected_ref_1 = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='mean_reference', bias_report=True, item=1, column_offset=column_offset ) * 100.
  expected_ref_2 = helpers.item_from_pipeline_result( fh=open(sys.argv[1]), y='mean_reference', bias_report=True, item=2, column_offset=column_offset ) * 100.
  plt.text( 108, 2.5, "Average bias (50bp): %.1f%%\nAverage bias (100bp): %.1f%%" % ( expected_ref_1, expected_ref_2 ), ha='right' )
  if include_unmapped:
    plt.legend(loc='upper left')
  else:
    plt.legend(loc='upper right')
  #plt.legend()
  plt.show()

if __name__ == '__main__':
  from matplotlib import rcParams
  rcParams.update({'figure.autolayout': True})
  #plot_series()
  #plot_error_bias()
  plot_reference_bias(include_zero=False)
  #plot_bias(include_zero=False, include_unmapped=True)
  #plot_bias_compare( include_unmapped=True, column_offset=0 )
  #plot_bias_compare( include_unmapped=False, column_offset=0 )
