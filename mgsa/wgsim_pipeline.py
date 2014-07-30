
import os
import sys

if len(sys.argv) < 2:
  print "Usage: %s multiplier (snp|ins|del|none) (bowtie2|bowtie2_e2e|bwa|bwasw)" % sys.argv[0]
  sys.exit(0)

mult = sys.argv[1]

if len(sys.argv) > 2:
  vartype = sys.argv[2]
else:
  vartype = 'wgsim'

if len(sys.argv) > 3:
  mapper = sys.argv[3]
else:
  mapper = 'bowtie2'

percent_var = 0.01
if vartype == 'snp':
  percent_vars = '%f 0 0' % percent_var
if vartype == 'ins':
  percent_vars = '0 %f 0' % percent_var
if vartype == 'del':
  percent_vars = '0 0 %f' % percent_var
if vartype == 'none':
  percent_vars = '0 0 0'

os.system( "wgsim ../../data/circoviridae.fasta -N 1000 ../../data/circoviridae_%s_x%s.fastq ../../circoviridae_wgsim_p2.fastq" % ( vartype, mult ) )
os.system( "python mapper_selector.py %s ../../data/circoviridae_x%s.fasta ../../data/circoviridae_%s_x%s.fastq ../../data/circoviridae_%s_x%s_%s.sam" % ( mapper, mult, vartype, mult, vartype, mult, mapper ) )
os.system( "wgsim_eval ../../data/circoviridae_%s_x%s_%s.sam" % ( vartype, mult, mapper ) )
