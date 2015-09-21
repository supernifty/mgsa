
import collections
import os
import random
import sys

import numpy as np

import bio

# usage
# cmd fasta
# config

# human
fasta_start = 20000000
fasta_finish = 40000000
step = 1000
prefix = 'human_chr21'

# pms2
#fasta_start = 6012800
#fasta_finish = 6049000
#step = 100
#prefix = 'human_pms2'

# drosophila
#fasta_start = 10000
#fasta_finish = 20000000
#step = 1000
#prefix = 'drosophila'

# malaria
#fasta_start = 10000
#fasta_finish = 900000
#step = 50

# ecoli
#fasta_start = 100000
#fasta_finish = 4000000
#step = 100

# hiv
#fasta_start = 1000
#fasta_finish = 10000
#step = 1
read_length = 100

fasta_file = sys.argv[1]
variation = sys.argv[2]
idx = random.randint( 0, 1e9 )
def run( cmd ):
  '''
    run a system command
  '''
  bio.log_stderr( cmd )
  os.system( cmd )

def read_fasta( fh ):
  result = []
  for line in fh:
    if line.startswith( '>' ):
      continue
    result.append( line.strip() )
  return ''.join( result )

print 'fasta_file:', fasta_file

# index it
#run( 'bwa index %s' % fasta_file )

with open( fasta_file, 'r' ) as fh:
  fasta = read_fasta( fh )

# make a variation free fastq
fastq_file = 'tmp%i.fq' % idx
with open( fastq_file, 'w' ) as fastq_fh:
  for fasta_pos in xrange( fasta_start, fasta_finish, step ):
    for offset in xrange( -read_length + 1, 1 ):
      read = fasta[fasta_pos + offset: fasta_pos + offset + read_length]
      fastq_fh.write( '@read-%i-%i\n%s\n+\n%s\n' % ( fasta_pos, fasta_pos + offset, read, '~' * len(read)  ) )

# do alignment
sam_file = 'tmp%i.sam' % idx
run( 'bwa mem -t 8 %s %s > %s' % ( fasta_file, fastq_file, sam_file ) )

# save original
#run( "samtools view -bS tmp%i.sam > original.bam" % idx )
#run( "samtools sort original.bam originalsorted" )
#run( "samtools index originalsorted.bam" )
#sys.exit(1)

# evaluate % mapped correctly
base_tp = collections.defaultdict(int)
base_fp = collections.defaultdict(int)
base_fn = collections.defaultdict(int)

with open( sam_file, 'r' ) as sam_fh:
  for line in sam_fh:
    if line.startswith( '@' ):
      continue
    fields = line.split('\t')
    start = int(fields[0].split('-')[1])
    correct = int(fields[0].split('-')[2])
    actual = int(fields[3]) - 1
    if correct == actual:
      base_tp[start] += 1
    elif actual == -1:
      base_fn[start] += 1
    else:
      base_fp[start] += 1

# make fastq with variation
fastq_file = 'tmp%i.fq' % idx
variation_length = int(variation[:-1])
variation_type = variation[-1]
with open( fastq_file, 'w' ) as fastq_fh:
  for fasta_pos in xrange( fasta_start, fasta_finish, step ):
    for offset in xrange( -read_length + 1, 1 ):
      if variation_type == 's':
        substitution = 'C' if fasta[fasta_pos] == 'A' else 'A'
        read = ''.join( ( fasta[fasta_pos + offset: fasta_pos], substitution, fasta[fasta_pos + 1: fasta_pos + offset + read_length] ) )
      elif variation_type == 'i':
        insertion = 'A' * variation_length
        read = ''.join( ( fasta[fasta_pos + offset: fasta_pos], insertion, fasta[fasta_pos: fasta_pos + read_length - variation_length + offset] ) )
      elif variation_type == 'd':
        read = ''.join( ( fasta[fasta_pos + offset: fasta_pos], fasta[fasta_pos + variation_length: fasta_pos + read_length + offset + variation_length] ) ) 
      else:
        read = fasta[fasta_pos + offset: fasta_pos + offset + read_length]

      fastq_fh.write( '@read-%i-%i\n%s\n+\n%s\n' % ( fasta_pos, fasta_pos + offset, read, '~' * len(read)  ) )

# do alignment
sam_file = 'tmp%i.sam' % idx
run( 'bwa mem %s %s > %s' % ( fasta_file, fastq_file, sam_file ) )

# evaluate % mapped correctly
var_tp = collections.defaultdict(int)
var_fp = collections.defaultdict(int)
var_fn = collections.defaultdict(int)

with open( sam_file, 'r' ) as sam_fh:
  for line in sam_fh:
    if line.startswith( '@' ):
      continue
    fields = line.split('\t')
    start = int(fields[0].split('-')[1])
    correct = int(fields[0].split('-')[2])
    actual = int(fields[3]) - 1
    if correct == actual: # right position
      if variation_type == 's' or variation.upper() in fields[5]:
        #print 'match', fields
        var_tp[start] += 1
      else: # variation not in cigar
        var_fn[start] += 1
    elif actual == -1: # unmapped
      var_fn[start] += 1
    else: # wrong position
      var_fp[start] += 1

#print 'base tp: %s fp: %s fn: %s' % ( base_tp, base_fp, base_fn )
#print 'var tp: %s fp: %s fn: %s' % ( var_tp, var_fp, var_fn )
def keep_highest( candidate, candidate_val, current, current_vals, count=5 ):
  if len(current) < count:
    current.append( candidate )
    current_vals.append( candidate_val )
  else:
    for i, v in enumerate(current):
      if candidate > v:
        current[i] = candidate
        current_vals[i] = candidate_val
        break

base_unmapped = 0
var_unmapped = 0
total = 0
base_worst_tp = 100
var_worst_tp = 100
var_worst_pos = -1
base_worst_pos = -1

change_worst_tp = []
change_worst_pos = []

change_tp = []

for fasta_pos in xrange( fasta_start, fasta_finish, step ):
  if base_tp[fasta_pos] == 0:
    base_unmapped += 1
  if var_tp[fasta_pos] == 0:
    var_unmapped += 1
  if base_tp[fasta_pos] < base_worst_tp:
    base_worst_tp = base_tp[fasta_pos]
    base_worst_pos = fasta_pos
  if var_tp[fasta_pos] < var_worst_tp:
    var_worst_tp = var_tp[fasta_pos]
    var_worst_pos = fasta_pos
  change = base_tp[fasta_pos] - var_tp[fasta_pos]
  change_tp.append( change )
  keep_highest( change_tp[-1], fasta_pos, change_worst_tp, change_worst_pos )
  total += 1

# add result to file
with open( '%s.result' % prefix, 'a' ) as fh:
  fh.write( '%s: %s\n' % ( variation, ','.join( [ str(s) for s in change_tp ] ) ) )

print "unmapped before %i/%i after %i/%i" % ( base_unmapped, total, var_unmapped, total )
print "worst tp before %i at %i after %i at %i biggest change %s at %s" % ( base_worst_tp, base_worst_pos, var_worst_tp, var_worst_pos, change_worst_tp, change_worst_pos )
print "change_tp max %i min %i ave %.1f sd %.1f" % ( max( change_tp ), min( change_tp ), np.mean( change_tp ), np.std( change_tp ) )

for i in xrange(0, len(change_worst_tp) ):
  # make minibam, and ready for viewing with igv
  run( "grep '\\(@\\|read-%i\\)' < tmp%i.sam > mini.sam" % ( change_worst_pos[i], idx ) )
  run( "samtools view -bS mini.sam > mini.bam" )
  run( "samtools sort mini.bam %s-%s-%i-%i" % ( prefix, variation, change_worst_pos[i], read_length ) )
  run( "samtools index %s-%s-%i-%i.bam" % ( prefix, variation, change_worst_pos[i], read_length ) )

run( "rm tmp%i.fq" % idx )
run( "rm tmp%i.sam" % idx )
