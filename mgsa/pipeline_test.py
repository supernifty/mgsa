
import os
import sys

if len(sys.argv) < 2:
  print "Usage: %s multiplier (snp|ins|del)" % sys.argv[0]
  sys.exit(0)

mult = sys.argv[1]

if len(sys.argv) > 2:
  vartype = sys.argv[2]
else:
  vartype = 'snp'

percent_var = 0.01
if vartype == 'snp':
  percent_vars = '%f 0 0' % percent_var
if vartype == 'ins':
  percent_vars = '0 %f 0' % percent_var
if vartype == 'del':
  percent_vars = '0 0 %f' % percent_var

#os.system( "python build_repeated_fasta.py ../../data/circoviridae.fasta %s > ../../data/circoviridae_x%s.fasta" % ( mult, mult ) )
#os.system( "python generate_mutation.py ../../data/circoviridae_x%s.fasta ../../data/circoviridae_%s_x%s.vcf %s > ../../data/circoviridae_%s_x%s.fasta" % ( mult, vartype, mult, percent_vars, vartype, mult ) )
os.system( "python generate_reads.py 50 10 ../../data/circoviridae_%s_x%s.vcf < ../../data/circoviridae_%s_x%s.fasta > ../../data/circoviridae_%s_x%s.fastq" % ( vartype, mult, vartype, mult, vartype, mult ) )
os.system( "python assemble.py ../../data/circoviridae_x%s.fasta ../../data/circoviridae_%s_x%s.fastq ../../data/circoviridae_%s_x%s.sam" % ( mult, vartype, mult, vartype, mult ) )
os.system( "python evaluate_reads.py < ../../data/circoviridae_%s_x%s.sam > ../../data/circoviridae_%s_x%s_evaluation.txt" % ( vartype, mult, vartype, mult ) )

