import collections
import itertools
import math
import random

import bio
import fasta
# before error_bias_no_variation
#y_map_old = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }

def fasta_stats( fh, log=bio.log_stderr ):
  return fasta.FastaStats( fh, log=log )

def fasta_set_width( fh, out, width=50, log=bio.log_stderr ):
  b = ''
  fr = fasta.FastaReader( fh, include_headers = True )
  for line in fr.items():
    if line.startswith( '>' ):
      if len(b) > 0:
        out.write( '%s\n',  b )
        b = ''
      out.write( '%s\n' % line )
    else:
      b += line
      while len(b) >= width:
        out.write( '%s\n' % b[:width] )
        b = b[width:]
  out.write( '%s\n' %  b )
    
def _update_range( line, start, command_list ):
  original_len = len(line)
  for i in xrange(len(command_list)-1, -1, -1): # assume sorted
    command = command_list[i]
    if start <= command[1] < start + original_len:
      # apply command
      if command[0] == 'i':
        line = line[:command[1]-start] + command[2] + line[command[1]-start:]
      elif command[0] == 'd':
        to_delete = int(command[2])
        if command[1] + to_delete > start + original_len:
          line = line[:command[1]-start]
          command_list[i] = ( 'd', start + original_len, to_delete - ( original_len - command[1] + start ) )
        else:
          line = line[:command[1]-start] + line[command[1] - start + to_delete:]
      else:
        print "unknown command", command
  return line
  
def fasta_edit( fh, out, commands, width=50 ):
  '''
    command is comma separated list of commands of the form [i|d].pos.[bases|count]
    e.g. i.50.A -> insert at pos 50, the base A
         d.80.20 -> delete at pos 80, 20 bases
  '''
  command_list = []
  for command in commands.split(','):
    new_command = command.split( '.' )
    new_command[1] = int(new_command[1]) # i|d, pos, bases/count
    command_list.append( new_command )

  b = ''
  fr = fasta.FastaReader( fh, include_headers = True )
  pos = 1
  for line in fr.items():
    if line.startswith( '>' ):
      if len(b) > 0:
        out.write( '%s\n' % b )
        b = ''
      out.write( '%s\n' % line )
    else:
      original_len = len(line)
      b += _update_range( line, pos, command_list )
      pos += original_len
      while len(b) >= width:
        out.write( '%s\n' % b[:width] )
        b = b[width:]
  out.write( '%s\n' %  b )
  
def fasta_filter( fh, out, commands ):
  '''
    command is comma separated list of locations of the form [i|d].pos.[bases|count]
    e.g. CHR1:123456-654321
    ignore CHR2 for now
  '''
  start, finish = commands.split('-')
  if ':' in start:
    start_chr, start_pos = start.split( ':' )
  else:
    start_chr = ''
    start_pos = start
  start_pos = int(start_pos) - 1
  finish_pos = int(finish) - 1

  fr = fasta.FastaReader( fh, include_headers = True )
  pos = 0
  chromosome = ''
  in_filter = False
  for line in fr.items():
    if line.startswith( '>' ):
      chromosome = line[1:].split(' ')[0]
      pos = 0
      if in_filter and start_chr != '' and chromosome != start_chr:
        in_filter = False
    else:
      line_len = len(line)
      if in_filter:
        if finish_pos + 1 >= pos + line_len:
          out.write( '%s' % line )
        else:
          out.write( '%s' % line[:finish_pos - pos + 1] )
          in_filter = False
      else: # not in_filter
        if ( start_chr == '' or chromosome == start_chr ) and start_pos >= pos and start_pos < pos + line_len:
          in_filter = True
          if finish_pos >= pos + line_len:
            out.write( '%s' % line[start_pos - pos:] )
          else:
            out.write( '%s' % line[start_pos - pos:finish_pos - pos + 1] )
            in_filter = False
      pos += line_len 
  out.write('\n')

def fasta_find( fh, query ):
  '''
    find query in fasta fh
    @return list of locations matching query
  '''
  fr = fasta.FastaReader( fh, include_headers=False )
  # for the moment, not memory efficient
  window = ''
  for line in fr.items():
    window += line
  results = []
  start = 0
  while True:
    cand = window.find( query, start )
    if cand >= 0:
      results.append( cand )
      start = cand + 1
    else:
      break
  return results

def poisson( lmbda ):
  '''
    return a poisson distribution
  '''
  l = math.exp( -lmbda )
  k = 0
  p = 1
  while p > l:
    k += 1
    p = p * random.random()
  
  return k - 1

def choose( n, k ):
  denom = math.factorial( k )
  num = 1
  for multiplier in xrange( n, n - k, -1 ):
    num *= multiplier
  return num / denom

def binomial_prob( k=0, m=1000, n=1000 ):
  '''
    what is the probability that, when there are m items to go in n slots, a slot will end up with k items
  '''
  p = 1. / n
  prob = p ** k * ( 1 - p ) ** ( m - k ) * choose( m, k )
  return prob

def binomial( m=1000, n=1000 ):
  '''
    how many items in a slot given m trials and n slots?
  '''
  r = random.random() # 0..1
  # prob of 0
  cumulative = 0.
  for k in xrange(0, m):
    cumulative += binomial_prob( k, m, n)
    if r <= cumulative:
      return k
  return m

def logistic( x ):
  return 1. / ( 1 + math.exp( -x ) )

def trim_edges( l ):
  start, end = find_edges( l )
  return l[start:end + 1]

def find_edges( l ):
  start = 0
  while l[start] == l[0]:
    start += 1
  end = -1
  while l[end] == l[-1]:
    end -= 1
  return start, end
  
def longest_run( l, x=0 ):
  '''
    longest run of x in list l
    e.g. ( [ 1, 0, 0, 1 ], 0 ) --> 2
  '''
  result = 0
  winning_pos = 0
  pos = 0
  for val, cand in itertools.groupby( l ):
    count = len(list(cand))
    if x == val:
      if count > result:
        result = count
        winning_pos = pos
    pos += count
  return result, winning_pos

def all_runs( l, x=0 ):
  '''
    find all runs of x in list l, 
    @return a dict of the form { 'length': 'count', ... }
  '''
  result = collections.defaultdict(int)
  for val, cand in itertools.groupby( l ):
    count = len(list(cand))
    if x == val:
      result[count] += 1
  return result
  
