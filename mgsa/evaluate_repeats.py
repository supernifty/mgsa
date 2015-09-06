
import glob
import os
import random
import re
import sys

import numpy as np

import bio
import config

def align( read, siteA, siteB, mapper_name, rm=True, alt_pos=1000, correct_pos=3000, single_insertion=False ):
  prefix = 'tmp/tmp%i' % random.randint(1, 1e9)
  # make fastq
  fastq = '%s.fastq' % prefix
  with open( fastq, 'w' ) as fq:
    data = [ '@A', read, '+', '~' * len( read ) ]
    fq.write( '\n'.join( data ) )

  # fasta - both variations
  #commands = 'i.1000.%s,i.3000.%s' % ( siteA, siteB )
  commands = 'i.%i.%s,i.%i.%s' % ( alt_pos, siteB, correct_pos, siteA )
  fasta = '%s.fa' % prefix
  bio.fasta_edit( open( sys.argv[1], 'r' ), open( fasta, 'w' ), commands )
  mapper = bio.get_mapper( config, mapper_name, fasta )

  # align
  sam = '%s.sam' % prefix
  mapper( fastq, sam )
  result_both = bio.SamFeatures( open( sam, 'r' ), None, classes="A", exclude_unmapped=False ).result[0]

  if single_insertion:
    # fasta - main variation only
    commands = 'i.1000.%s' % ( siteA )
    fasta = '%s.fa' % prefix
    bio.fasta_edit( open( sys.argv[1], 'r' ), open( fasta, 'w' ), commands )
    mapper = bio.get_mapper( config, mapper_name, fasta )

    # align
    sam = '%s.sam' % prefix
    mapper( fastq, sam )
    result_single = bio.SamFeatures( open( sam, 'r' ), None, classes="A", exclude_unmapped=False ).result[0]
  else:
    result_single = None

  #result = bio.SamFeatures( open( sam, 'r' ), None, classes="A" ).result[0]
  if rm:
    [ os.remove( x ) for x in glob.glob( '%s*' % prefix ) ]
  return result_both, result_single

# 
# experiment to see what position does with one mutation
# what is the effect of position?
#
def experiment_pos( result_fh ):
  #result_fh.write( "%s,%s,%s,%s,%s,%s,%s\n" % ( "read", "location", "alt", "mapper", "mapq_post", "mut_pos", "align_pos" ) )
  result_fh.write( "%s,%s,%s,%s\n" % ( "mapper", "mapq_post", "mut_pos", "align_pos" ) )
  read_len = 100
  max_cand_mutations = 10
  max_pos = 30 #100
  matrix = { 'bwa': {}, 'bowtie2': {} }
  sample = 100
  gc = 0.5
  for _ in xrange(sample):
    s1 = bio.SequenceGenerator( read_len, (0.5-gc/2, gc/2, gc/2, 0.5-gc/2) )
    s2 = bio.SequenceGenerator( 0 )
    s2.sequence = s1.sequence
    # 1 mutation
    for mapper_name in matrix.keys(): 
      for p in xrange(max_pos):
        s2.sequence = s1.sequence
        # just the 1 mutation to the read
        read = bio.SequenceGenerator( 0 )
        read.sequence = s1.sequence
        read.mutate( p )
        # apply mutations to candidate
        for cand_mutations in xrange(max_cand_mutations):
          cand_mutation_list = []
          s2.sequence = s1.sequence
          for m in xrange(cand_mutations):
            #m_pos = 10 + ( read_len - 10 ) / cand_mutations * m + 5
            while True:
              #m_pos = random.randint(10, read_len - 1)
              #m_pos = random.randint(0, 20)
              m_pos = random.randint(10, read_len - 10)
              if m_pos not in cand_mutation_list and m_pos != p:
                break
            #if m_pos == p:
            #  m_pos -= 2
            cand_mutation_list.append( m_pos )
            s2.mutate( cand_mutation_list[-1] )
          r = align( read.sequence, s1.sequence, s2.sequence, mapper_name )
          result_fh.write( "%s,%s,%i,%s,%f,%s,%s\n" % ( mapper_name, r[0]['mapq'], p, r[0]['pos'], gc, r[0]['cigar'], '-'.join([ str(x) for x in cand_mutation_list]) ) )
          key = '%i,%i' % ( cand_mutations, p )
          if key not in matrix[mapper_name]:
            matrix[mapper_name][key] = []
          matrix[mapper_name][key].append( int(r[0]['mapq']) )
          #matrix[mapper_name]['%i,%i' % ( cand_mutations, p ) ] = int(r[0]['mapq'])
  result_fh.write( 'x axis is candidate mutations, y axis is position of mutation on read\n' )
  show_dist_matrix( result_fh, matrix, max_cand_mutations, max_pos )

# 
# experiment to see what position does with one mutation
# what is the effect of position?
#
def experiment_genome_pos( result_fh ):
  #result_fh.write( "%s,%s,%s,%s,%s,%s,%s\n" % ( "read", "location", "alt", "mapper", "mapq_post", "mut_pos", "align_pos" ) )
  result_fh.write( "%s,%s,%s,%s\n" % ( "mapper", "mapq_post", "mut_pos", "align_pos" ) )
  read_len = 100
  offset_step = 25
  matrix = { 'bwa': {}, 'bowtie2': {} }
  max_cand_mutations = 6
  gc = 0.5
  samples = 100
  for _ in xrange(samples):
    s1 = bio.SequenceGenerator( read_len, (0.5-gc/2, gc/2, gc/2, 0.5-gc/2) )
    s2 = bio.SequenceGenerator( 0 )
    s2.sequence = s1.sequence
    p = 20
    # 1 mutation
    for mapper_name in matrix.keys(): 
      offset_count = 0
      for offset in xrange(0, 1000, offset_step):
        s2.sequence = s1.sequence
        # just the 1 mutation to the read
        read = bio.SequenceGenerator( 0 )
        read.sequence = s1.sequence
        read.mutate( p )
        # apply mutations to candidate
        for cand_mutations in xrange(max_cand_mutations):
          cand_mutation_list = []
          s2.sequence = s1.sequence
          for m in xrange(cand_mutations):
            #m_pos = 10 + ( read_len - 10 ) / cand_mutations * m + 5
            while True:
              #m_pos = random.randint(10, read_len - 1)
              #m_pos = random.randint(0, 20)
              m_pos = random.randint(10, read_len - 10)
              if m_pos not in cand_mutation_list and m_pos != p:
                break
              #if m_pos == p:
              #  m_pos -= 2
            cand_mutation_list.append( m_pos )
            s2.mutate( cand_mutation_list[-1] )
          r = align( read.sequence, s1.sequence, s2.sequence, mapper_name, alt_pos=1000+offset, correct_pos=3000-offset )
          result_fh.write( "%s,%s,%i,%s,%f,%s,%s\n" % ( mapper_name, r[0]['mapq'], p, r[0]['pos'], gc, r[0]['cigar'], '-'.join([ str(x) for x in cand_mutation_list]) ) )
          key = '%i,%i' % ( cand_mutations, offset_count )
          if key not in matrix[mapper_name]:
            matrix[mapper_name][key] = []
          matrix[mapper_name][key].append( int(r[0]['mapq']) )
        offset_count += 1
      #matrix[mapper_name]['%i,%i' % ( cand_mutations, p ) ] = int(r[0]['mapq'])
  result_fh.write( 'x axis is candidate mutations, y axis is offset\n' )
  show_dist_matrix( result_fh, matrix, max_cand_mutations, offset_count )

#
# what is the mapping quality as the distance between repeats and the read change?
#
def experiment_candidate_distance( result_fh ):
  result_fh.write( "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ( "read", "location", "alt", "mapper", "mapq_post", "mapq_pre", "pos_post", "pos_pre", "cand_mutations", "read_mutations", "post_cigar" ) )
  read_len = 100
  max_read_mutations = 12
  max_cand_mutations = 12
  samples = 100
  gc = 0.5
  #for gc in np.arange( 0., 1., 0.2 ):
  matrix = { 'bwa': {}, 'bowtie2': {} }
  for sample in xrange(samples):
    s1 = bio.SequenceGenerator( read_len, (0.5-gc/2, gc/2, gc/2, 0.5-gc/2) )
    s2 = bio.SequenceGenerator( 0 )
    s2.sequence = s1.sequence
    # 1 mutation
    for mapper_name in ('bwa', 'bowtie2'):
      #for p in xrange(0, read_len, 10):
      #for p in xrange(10):
      #  s2.sequence = s1.sequence
      #  s2.mutate( p )
      #  r = align( s1.sequence, s1.sequence, s2.sequence, mapper_name )
      #  #result_fh.write( "%s,%s,%s,%s,%s\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r ) )
      #  result_fh.write( "%s,%s,%s,%s,%i,%i,1\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r[0], r[1] ) )
      for read_mutations in xrange(max_read_mutations):
        read = bio.SequenceGenerator( 0 )
        read.sequence = s1.sequence
        # add mutations to read
        read_mutation_list = []
        for m in xrange(read_mutations):
          while True:
            m_pos = random.randint( 10, read_len - 10 )
            if m_pos not in read_mutation_list:
              break
          read_mutation_list.append( m_pos ) # don't mutate in the same place
          read.mutate( read_mutation_list[-1] ) 
        for cand_mutations in xrange(max_cand_mutations):
          # add mutations to alternative
          cand_mutation_list = []
          s2.sequence = s1.sequence
          for m in xrange(cand_mutations):
            while True:
              m_pos = random.randint( 10, read_len - 10 )
              if m_pos not in read_mutation_list and m_pos not in cand_mutation_list:
                break
            cand_mutation_list.append( m_pos )
            s2.mutate( cand_mutation_list[-1] )
            s2.mutate( cand_mutation_list[-1] ) # mutate twice so it is still different from cand
          # align and add result
          r = align( read.sequence, s1.sequence, s2.sequence, mapper_name ) # only populates r[0]
          #result_fh.write( "%s,%s,%s,%s,%s\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r ) )
          result_fh.write( "%i,%s,%i,%i,%s,%s,%s\n" % ( sample, mapper_name, int(r[0]['mapq']), int(r[0]['pos']), '-'.join( [ str(x) for x in cand_mutation_list ] ), '-'.join( [ str(x) for x in read_mutation_list ] ), r[0]['cigar'] ) )
          #result_fh.write( "%s,%s,%s,%s,%i,%i,%i,%i,%s,%s,%s\n" % ( read.sequence, s1.sequence, s2.sequence, mapper_name, int(r[0]['mapq']), int(r[1]['mapq']), int(r[0]['pos']), int(r[1]['pos']), cand_mutations, read_mutations, r[0]['cigar'] ) )
          key = '%i,%i' % ( cand_mutations, read_mutations ) 
          if key not in matrix[mapper_name]:
            matrix[mapper_name][key] = []
          matrix[mapper_name][key].append( int(r[0]['mapq']) )
  result_fh.write( 'x axis is candidate mutations, y axis is read mutations\n' )
  show_dist_matrix( result_fh, matrix, max_cand_mutations, max_read_mutations )

#
# what is the mapping quality with ultra low entropy
#
def experiment_entropy( result_fh ):
  result_fh.write( "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ( "read", "location", "alt", "mapper", "mapq_post", "mapq_pre", "pos_post", "pos_pre", "cand_mutations", "read_mutations", "post_cigar" ) )
  read_len = 100
  #max_read_mutations = 3
  max_cand_mutations = 8
  samples = 3 # 100
  entropy_max = 0.1
  #gc = 0.2
  #for gc in np.arange( 0., 1., 0.2 ):
  matrix = { 'bwa': {}, 'bowtie2': {} }
  for sample in xrange(samples):
    for entropy in np.arange( 0., entropy_max, 0.01 ):
      s1 = bio.SequenceGenerator( read_len, (1 - entropy, entropy / 3, entropy / 3, entropy / 3 ) )
      s2 = bio.SequenceGenerator( 0 )
      s2.sequence = s1.sequence
      # 1 mutation
      for mapper_name in matrix.keys():
        #for p in xrange(0, read_len, 10):
        #for p in xrange(10):
        #  s2.sequence = s1.sequence
        #  s2.mutate( p )
        #  r = align( s1.sequence, s1.sequence, s2.sequence, mapper_name )
        #  #result_fh.write( "%s,%s,%s,%s,%s\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r ) )
        #    result_fh.write( "%s,%s,%s,%s,%i,%i,1\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r[0], r[1] ) )
        # 1 read mutation
        read = bio.SequenceGenerator( 0 )
        read.sequence = s1.sequence
        # add mutations to read
        read_mutation_list = []
        while True:
          m_pos = random.randint( 10, read_len - 10 )
          if m_pos not in read_mutation_list:
            break
        read_mutation_list.append( m_pos ) # don't mutate in the same place
        read.mutate( read_mutation_list[-1] ) 
        for cand_mutations in xrange(max_cand_mutations):
          # add mutations to alternative
          cand_mutation_list = []
          s2.sequence = s1.sequence
          for m in xrange(cand_mutations):
            while True:
              m_pos = random.randint( 10, read_len - 10 )
              if m_pos not in read_mutation_list and m_pos not in cand_mutation_list:
                break
            cand_mutation_list.append( m_pos )
            s2.mutate( cand_mutation_list[-1] )
            s2.mutate( cand_mutation_list[-1] ) # mutate twice so it is still different from cand
          # align and add result
          r = align( read.sequence, s1.sequence, s2.sequence, mapper_name ) # only populates r[0]
          #result_fh.write( "%s,%s,%s,%s,%s\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r ) )
          result_fh.write( "%i,%s,%i,%i,%s,%s,%s,%s,%s,%s\n" % ( sample, mapper_name, int(r[0]['mapq']), int(r[0]['pos']), '-'.join( [ str(x) for x in cand_mutation_list ] ), '-'.join( [ str(x) for x in read_mutation_list ] ), r[0]['cigar'], read.sequence, s1.sequence, s2.sequence ) )
          #result_fh.write( "%s,%s,%s,%s,%i,%i,%i,%i,%s,%s,%s\n" % ( read.sequence, s1.sequence, s2.sequence, mapper_name, int(r[0]['mapq']), int(r[1]['mapq']), int(r[0]['pos']), int(r[1]['pos']), cand_mutations, read_mutations, r[0]['cigar'] ) )
          key = '%i,%i' % ( cand_mutations, int(entropy*100) ) 
          if key not in matrix[mapper_name]:
            matrix[mapper_name][key] = []
          matrix[mapper_name][key].append( int(r[0]['mapq']) )
  result_fh.write( 'x axis is candidate mutations, y axis is entropy\n' )
  show_dist_matrix( result_fh, matrix, max_cand_mutations, int(entropy_max*100) )

#
# what is the mapping quality with tandem repeats
#
def experiment_tandem_repeats( result_fh ):
  result_fh.write( "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % ( "read", "location", "alt", "mapper", "mapq_post", "mapq_pre", "pos_post", "pos_pre", "cand_mutations", "read_mutations", "post_cigar" ) )
  read_len = 100
  #max_read_mutations = 3
  max_cand_mutations = 8
  samples = 100
  repeat_max = 40
  gc = 0.5
  #for gc in np.arange( 0., 1., 0.2 ):
  matrix = { 'bwa': {}, 'bowtie2': {} }
  for sample in xrange(samples):
    for repeat in np.arange( 1, repeat_max ):
      s1 = bio.SequenceGenerator( read_len, (0.5-gc/2, gc/2, gc/2, 0.5-gc/2), repeat_length=repeat )
      s2 = bio.SequenceGenerator( 0 )
      s2.sequence = s1.sequence
      # 1 mutation
      for mapper_name in matrix.keys():
        #for p in xrange(0, read_len, 10):
        #for p in xrange(10):
        #  s2.sequence = s1.sequence
        #  s2.mutate( p )
        #  r = align( s1.sequence, s1.sequence, s2.sequence, mapper_name )
        #  #result_fh.write( "%s,%s,%s,%s,%s\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r ) )
        #    result_fh.write( "%s,%s,%s,%s,%i,%i,1\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r[0], r[1] ) )
        # 1 read mutation
        read = bio.SequenceGenerator( 0 )
        read.sequence = s1.sequence
        # add mutations to read
        read_mutation_list = []
        while True:
          m_pos = random.randint( 10, read_len - 10 )
          if m_pos not in read_mutation_list:
            break
        read_mutation_list.append( m_pos ) # don't mutate in the same place
        read.mutate( read_mutation_list[-1] ) 
        for cand_mutations in xrange(max_cand_mutations):
          # add mutations to alternative
          cand_mutation_list = []
          s2.sequence = s1.sequence
          for m in xrange(cand_mutations):
            while True:
              m_pos = random.randint( 10, read_len - 10 )
              if m_pos not in read_mutation_list and m_pos not in cand_mutation_list:
                break
            cand_mutation_list.append( m_pos )
            s2.mutate( cand_mutation_list[-1] )
            s2.mutate( cand_mutation_list[-1] ) # mutate twice so it is still different from cand
          # align and add result
          r = align( read.sequence, s1.sequence, s2.sequence, mapper_name ) # only populates r[0]
          #result_fh.write( "%s,%s,%s,%s,%s\n" % ( s1.sequence, s1.sequence, s2.sequence, mapper_name, r ) )
          result_fh.write( "%i,%s,%i,%i,%s,%s,%s,%s\n" % ( sample, mapper_name, int(r[0]['mapq']), int(r[0]['pos']), '-'.join( [ str(x) for x in cand_mutation_list ] ), '-'.join( [ str(x) for x in read_mutation_list ] ), r[0]['cigar'], read.sequence ) )
          #result_fh.write( "%s,%s,%s,%s,%i,%i,%i,%i,%s,%s,%s\n" % ( read.sequence, s1.sequence, s2.sequence, mapper_name, int(r[0]['mapq']), int(r[1]['mapq']), int(r[0]['pos']), int(r[1]['pos']), cand_mutations, read_mutations, r[0]['cigar'] ) )
          key = '%i,%i' % ( cand_mutations, repeat ) 
          if key not in matrix[mapper_name]:
            matrix[mapper_name][key] = []
          matrix[mapper_name][key].append( int(r[0]['mapq']) )
  result_fh.write( 'x axis is candidate mutations, y axis is repeat size\n' )
  show_dist_matrix( result_fh, matrix, max_cand_mutations, repeat_max, ymin=1 )

def count_unmapped( fh_in, fh_out ):
  matrix = { 'bwa': {}, 'bowtie2': {} }
  xmax = 0
  ymax = 0
  first = True
  for line in fh_in:
    if first:
      first = False
      continue
    fields = line.strip().split(',')
    x = len(filter( None, re.split('-', fields[4]) ) )
    y = len(filter( None, re.split('-', fields[5]) ) )
    xmax = max(xmax, x)
    ymax = max(ymax, y)
    is_unmapped = fields[6] == '*'
    key = '%i,%i' % (x, y)
    if key not in matrix[fields[1]]:
      matrix[fields[1]][key] = []
    matrix[fields[1]][key].append( 0 if is_unmapped else 100 )
    
  fh_out.write( 'x axis is candidate mutations, y axis is read mutations\n' )
  show_dist_matrix( fh_out, matrix, xmax+1, ymax+1)

def count_clipped( fh_in, fh_out ):
  matrix = { 'bwa': {}, 'bowtie2': {} }
  xmax = 0
  ymax = 0
  first = True
  for idx, line in enumerate(fh_in):
    if first:
      first = False
      continue
    fields = line.strip().split(',')
    if len(fields) > 5:
      x = len(filter( None, re.split('-', fields[6]) ) )
      y = int(fields[2])
      xmax = max(xmax, x)
      ymax = max(ymax, y)
      is_clipped = fields[5] != '100M'
      key = '%i,%i' % (x, y)
      if key not in matrix[fields[0]]:
        matrix[fields[0]][key] = []
      matrix[fields[0]][key].append( 100 if is_clipped else 0 )
      print "line", idx
    
  fh_out.write( 'x axis is candidate mutations, y axis is read mutations\n' )
  show_dist_matrix( fh_out, matrix, xmax+1, ymax+1)


def show_dist_matrix( fh, matrix, xmax, ymax, ymin=0 ):
  fh.write( 'Mean Quality' )
  for mapper_name in matrix:
    fh.write( '\n' + mapper_name )
    fh.write( '\t' + '\t'.join( [ str(x) for x in xrange( xmax ) ] ) + '\n' )
    for y in xrange(ymin, ymax): 
      fh.write( str(y) + '\t' )
      xs = []
      for x in xrange(xmax): 
        xs.append( np.mean( matrix[mapper_name]['%i,%i' % (x,y) ] ) )
      fh.write( '\t'.join( [ '%.1f' % x for x in xs ] ) + '\n' )
  fh.write( '\nSD' )
  for mapper_name in matrix:
    fh.write( '\n' + mapper_name )
    fh.write( '\t' + '\t'.join( [ str(x) for x in xrange( xmax ) ] ) + '\n' )
    for y in xrange(ymin, ymax): 
      fh.write( str(y) + '\t' )
      xs = []
      for x in xrange(xmax): 
        xs.append( np.std( matrix[mapper_name]['%i,%i' % (x,y) ] ) )
      fh.write( '\t'.join( [ '%.1f' % x for x in xs ] ) + '\n' )

# make fasta
with open(sys.argv[2], 'w') as result_fh:
  if sys.argv[3] == 'position':
    experiment_pos( result_fh )
  elif sys.argv[3] == 'genome':
    experiment_genome_pos( result_fh )
  elif sys.argv[3] == 'distance':
    experiment_candidate_distance( result_fh )
  elif sys.argv[3] == 'entropy':
    experiment_entropy( result_fh )
  elif sys.argv[3] == 'tandem':
    experiment_tandem_repeats( result_fh )
  elif sys.argv[3] == 'unmapped':
    count_unmapped( sys.stdin, result_fh )
  elif sys.argv[3] == 'clip':
    count_clipped( sys.stdin, result_fh )
