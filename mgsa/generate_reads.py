#
# generate reads from a fasta file
# includes tags in each sequence comment of the correct position
#
# python generate_reads.py 50 10 mutations.vcf < donor.fasta > donor.fastq
#
import os
import random
import sys
import helpers

import bio

def write( pos, sequence, quality, variations=set(), extra='' ):
  if len(variations) == 0:
    sys.stdout.write( '@mgsa_seq_%i_%s\n' % ( (pos+1), extra ) ) # sam is 1 indexed
  else:
    sys.stdout.write( '@mgsa_seq_%i_variation_%s_%s\n' % ( (pos+1), ','.join([ v for v in variations ]), extra ) ) # sam is 1 indexed
  sys.stdout.write( sequence )
  sys.stdout.write( '\n+\n' )
  sys.stdout.write( quality ) # quality 
  sys.stdout.write( '\n' )


if len(sys.argv) < 3:
  print "Usage: %s readlength coverage [vcf]" % sys.argv[0]
  sys.exit(1)

read_length = int(sys.argv[1])
coverage = int(sys.argv[2])
if len(sys.argv) == 4:
  vcf = bio.VCF( reader=open(sys.argv[3], 'r') )
else:
  vcf = None

probabilistic = False

read_probability = 1. * coverage / read_length 
read_every = 1. / read_probability

# read candidate fasta
dna = ''
base_candidate_position = 0
time_to_next = read_every
max_pos = 0
for line in sys.stdin:
  if line.startswith( '>' ):
    continue
  line = line.strip()
  dna += line
  process = len(dna) - read_length + 1
  if process > 0:
    for i in xrange(0, process):
      time_to_next -= 1
      if probabilistic and random.random() < read_probability or \
         not probabilistic and time_to_next <= 0:
        candidate_position = base_candidate_position + i # reference == candidate
        if vcf is None:
          reference_position = candidate_position # reference == candidate
          variations = set()
        else:
          net = vcf.net_insertions( candidate_position ) # net_insertions expects a reference position
          reference_position = candidate_position - net # guess reference position
          reference_net = vcf.net_insertions( reference_position ) # get variation count
          move = 0
          if net != reference_net: # there was a variation between reference_position and candidate_position, e.g. net=4, ref=3
            #net = reference_net
            reference_net = vcf.net_insertions( reference_position + net - reference_net )
            #reference_position = reference_position + net - reference_net 
            reference_position = candidate_position - reference_net
            #print "net %i new_net %i new_position %i old_pos %i" % ( net, reference_net, reference_position, candidate_position - net )
            move += 1
          # on an insertion?
          if vcf.find_indel( reference_position ) == None:
            include_start = False
            #print "no indel"
          else: # either on or just before?
            net_position = candidate_position - reference_position 
            net_insertions = vcf.net_insertions( reference_position ) 
            include_start = ( net_insertions >= net_position )
            #print "net_pos %i net_insertions %i -> %s" % ( net_position, net_insertions, include_start )
            
          #print "net for %i: %i" % ( ( base_candidate_position + i ), net )
          variations = vcf.variations( reference_position, reference_position + read_length, include_start=include_start )
        write( reference_position, dna[i:i+read_length], '~' * read_length, variations=variations, extra='net%i_for%i_is%s' % (reference_net, reference_position, include_start ) )
        max_pos = max(max_pos, reference_position + read_length)
        time_to_next += read_every
    base_candidate_position += process
    dna = dna[process:]

# process final fragment
#if len(dna) > 0:
#  for i in xrange(0, len(dna)):
#    time_to_next -= 1
#    if probabilistic and random.random() < read_probability:
#      write( position + i, dna[i:i+read_length], '~' * read_length )
#      max_pos = max(max_pos, position + i + read_length)
#    if not probabilistic and time_to_next <= 0:
#      write( position + i, dna[i:i+read_length], '~' * read_length )
#      max_pos = max(max_pos, position + i + read_length)
#      time_to_next += read_every
#  position += len(dna)

#print max_pos, " basepairs processed.", len(dna), " unprocessed."
