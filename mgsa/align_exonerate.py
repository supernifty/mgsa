
# align unmapped reads in a bam file using exonerate

import argparse
import subprocess
import sys

import bio
import config

cmd = 'exonerate --model affine:local tmp.fa %s'

parser = argparse.ArgumentParser(description='Align exonerate')
parser.add_argument('--bam', dest='bam', required=True, help='bam file to align')
parser.add_argument('--ref', dest='ref', required=True, help='target reference')
parser.add_argument('--count', dest='count', type=int, default=1e9, help='stop after trying this many reads')
args = parser.parse_args()
sam_fh = bio.BamReaderExternal( config.BAM_TO_SAM, args.bam )
count = args.count
unmapped = 0
still_unmapped = 0
mapped = 0
scores = []
best = ( 0, [] )

for idx, line in enumerate( sam_fh ):
  fields = line.strip().split()
  if len(fields) > 9:
    # is unmapped?
    flag = int( fields[1] )
    mapq = int( fields[4] )
    if flag & 0x04 != 0: # unmapped
      unmapped += 1
      # try to map it
      # write query
      with open( 'tmp.fa', 'w' ) as q:
        q.write( '>tmp\n%s' % fields[9] )

      what = cmd % args.ref
      #print "executing", what
      p = subprocess.Popen(what, shell=True, bufsize=0, stdin=subprocess.PIPE, stdout=subprocess.PIPE) #, close_fds=True)
      found = False
      current = []
      for line in p.stdout:
        current.append( line )
        if line.find( 'Raw score' ) != -1:
          score = int( line.split(':')[1].strip() )
          scores.append( score )
          found = True
          if score > best[0]: # read the whole response
            best = ( score, current )
          else:
            break
      if found:
        mapped += 1
      else:
        still_unmapped += 1
      
      if unmapped % 100 == 0:
        print "processed %i lines: unmapped %i mapped %i still unmapped %i" % ( idx, unmapped, mapped, still_unmapped )
      count -= 1
      if count == 0:
        break      
  if idx % 100000 == 0:
    print "processed %i lines: unmapped %i mapped %i still unmapped %i" % ( idx, unmapped, mapped, still_unmapped )

print "==========================================================="
print "processed %i lines: unmapped %i mapped %i still unmapped %i" % ( idx, unmapped, mapped, still_unmapped )
if len(scores) > 0:
  print "mean %.2f max %i min %i hit %.2f%%" % ( 1. * sum(scores) / len(scores), max(scores), min(scores), 100. * mapped / unmapped )
  print "==========================================================="
  print "best alignment:\n", ''.join( best[1] )
