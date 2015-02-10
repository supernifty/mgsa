
import bio
import vcf 

def multi_chromosome_comparison( parent, child, truth ):
  '''
    params are all MultiChromosomeVCFs
  '''
  offset = 0

  chromosome_offsets = []
  chromosome_names = []
  x_tp, x_fp, x_fn, y_tp, y_fp, y_fn, x_t, y_t = [], [], [], [], [], [], [], []
  sorted_chromosomes = child.vcfs.keys()
  sorted_chromosomes.sort()

  count = 0
  for chromosome in sorted_chromosomes:
    diff = vcf.VCFDiff( vcf_correct=parent.vcfs[chromosome], vcf_candidate=child.vcfs[chromosome], generate_positions=True, log=bio.log_quiet )
    x_tp_local, y_tp_local = diff.confidence_map( 'tp' ) # parent and child
    x_fp_local, y_fp_local = diff.confidence_map( 'fp' ) # child
    x_fn_local, y_fn_local = diff.confidence_map( 'fn' ) # parent
    x_t_local = []
    if chromosome in truth.vcfs:
      for i in xrange(len(x_fp_local)):
        if x_fp_local[i] in truth.vcfs[chromosome].snp_map:
          x_t_local.append( x_fp_local[i] )
          y_t.append( y_fp_local[i] )
    x_tp.extend( [ x + offset for x in x_tp_local ] )
    x_fp.extend( [ x + offset for x in x_fp_local ] )
    x_fn.extend( [ x + offset for x in x_fn_local ] )
    x_t.extend( [ x + offset for x in x_t_local ] )
    y_tp.extend( y_tp_local )
    y_fp.extend( y_fp_local )
    y_fn.extend( y_fn_local )
    count += 1
    chromosome_names.append( ( offset, count % 2, chromosome.split('|')[-1] ) )
    offset += diff.length
    chromosome_offsets.append( offset )
    s = []
    for i in xrange(len(x_fp_local)):
      if y_fp_local[i] >= 7:
        s.append( '%i: %i' % ( x_fp_local[i], y_fp_local[i] ) )

  return { 'x_tp': x_tp, 'x_fp': x_fp, 'x_fn': x_fn, 'y_tp': y_tp, 'y_fp': y_fp, 'y_fn': y_fn, 'x_t': x_t, 'y_t': y_t, 'chromosome_offsets': chromosome_offsets, 'chromosome_names': chromosome_names, 'offset': offset }

