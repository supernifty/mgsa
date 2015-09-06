
import config
import datetime
import os
import sys

import bio

# set these
insertion = False
mapper = 'bwa_indel' # bowtie2
fasta = 'e-coli-mg1655'

# all good
if insertion:
  tmpl = 'batch/find_longest_indel.cfg'
  start = 20
else:
  tmpl = 'batch/find_longest_deletion.cfg'
  start = 20
cfg = 'tmp.cfg'
target = 'tmp.out'
results = {}
recall_field = 21

#start = 80

def run( cmd ):
  '''
    run a system command
  '''
  bio.log_stderr( cmd )
  os.system( cmd )

for rl in (100, 200, 300, 400, 500, 750, 1000):
#for rl in (200, 300, 400, 500, ):
  while True:
    run( "sed 's/read_length 100/read_length %i/g' < %s > %s" % ( rl, tmpl, cfg ) )
    #start += 1
    start += 5
    run( "sed -i.bck 's/mapper bowtie2/mapper %s/g' %s" % ( mapper, cfg ) )
    run( "sed -i.bck 's/min_variation_distance 100/min_variation_distance %i/g' %s" % ( rl, cfg ) )
    run( "sed -i.bck 's/fasta hiv/fasta %s/g' %s" % ( fasta, cfg ) )
    if insertion:
      run( "sed -i.bck 's/min_insertion_len 25/min_insertion_len %i/g' %s" % ( start, cfg ) )
      run( "sed -i.bck 's/max_insertion_len 25/max_insertion_len %i/g' %s" % ( start, cfg ) )
    else: # deletion
      run( "sed -i.bck 's/min_deletion_len 25/min_deletion_len %i/g' %s" % ( start, cfg ) )
      run( "sed -i.bck 's/max_deletion_len 25/max_deletion_len %i/g' %s" % ( start, cfg ) )
    run( "python pipeline_batch.py %s %s" % ( cfg, target ) )
    lines = open( target, 'r' ).readlines()
    recall = float(lines[2].split(',')[recall_field])
    if recall < 0.01:
      results[rl] = start
      start -= 5
      break

print results
