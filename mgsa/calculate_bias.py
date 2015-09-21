
# the entire process of aligning and calculating bias
# parameters
# - donor: fasta file of target
# - reference: fasta file of reference
# - fastq: fastq files

import argparse
import os
import random

import bio
import config

def run( cmd ):
  '''
    run a system command
  '''
  bio.log_stderr( cmd )
  os.system( cmd )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('--donor', help='donor fasta')
  parser.add_argument('--reference', help='reference fasta')
  parser.add_argument('--job', required=False, help='use to continue previous pipeline')
  parser.add_argument('--start', required=False, help='start from this stage')
  parser.add_argument('--tmpdir', required=False, help='where to write files')
  parser.add_argument('fastq', help='fastq files to align')
  args = parser.parse_args()
  # now do each stage...
  if args.job:
    idx = int(args.job)
  else:
    idx = random.randint(1, 1e6)
  if args.start:
    start = int(args.start)
  else:
    start = 0
  if args.tmpdir:
    tmpdir = args.tmpdir
  else:
    tmpdir = '/tmp'
  bio.log_stderr( 'Job ID: %i, Starting at stage %i' % (idx, start) )
  # TODO error correction (ec)

  # fasta indexes
  stage = 1
  if start <= stage:
    run( '%s index %s' % ( config.BWA_PATH, args.donor ) )
    run( '%s index %s' % ( config.BWA_PATH, args.reference ) )
    bio.log_stderr( 'Stage %i: Indexing completed' % stage )

  stage += 1 # 2
  # alignment (aln)
  if start <= stage:
    run( '%s mem -t 8 %s %s > %s/donor%i.sam' % ( config.BWA_PATH, args.donor, args.fastq, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Donor alignment completed' % stage )

  stage += 1 # 3
  if start <= stage:
    run( '%s mem -t 8 %s %s > %s/reference%i.sam' % ( config.BWA_PATH, args.reference, args.fastq, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Reference alignment completed' % stage )

  # genome alignment (mauve)
  stage += 1 # 4
  if start <= stage:
    run( '%s --output=%s/mauve%i %s %s' % ( config.MAUVE_PATH, tmpdir, idx, args.donor, args.reference ) )
    bio.log_stderr( 'Stage %i: Mauve completed' % stage )

  # realignment
  stage += 1 # 5
  if start <= stage:
    donor_accession = open( args.donor, 'r' ).readline().strip().split()[0][1:]
    run( 'python remap_bam.py --xmfa %s/mauve%i --origin 2 --target 1 --output_not_covered %s/notcovered%i.sam --output %s/remapped%i.sam %s/reference%i.sam --new_reference \'%s\' > %s/remap_bam%i.stats' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, donor_accession, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Remap completed' % stage )

  # converto to bam
  stage += 1 # 6
  if start <= stage:
    run( 'samtools view -bS %s/donor%i.sam > %s/donor%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools view -bS %s/reference%i.sam > %s/reference%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    # fix remapped
    with open( '%s/donor%i.sam' % ( tmpdir, idx ), 'r' ) as dfh:
      l = (dfh.readline(), dfh.readline())
    with open( '%s/remapped%i.head' % ( tmpdir, idx ), 'w' ) as rfh:
      rfh.write( l[0] )
      rfh.write( l[1] )
    run( 'cat %s/remapped%i.head %s/remapped%i.sam | samtools view -bS - > %s/remapped%i.bam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    with open( '%s/reference%i.sam' % ( tmpdir, idx ), 'r' ) as dfh:
      l = (dfh.readline(), dfh.readline())
    with open( '%s/notcovered%i.head' % ( tmpdir, idx ), 'w' ) as rfh:
      rfh.write( l[0] )
      rfh.write( l[1] )
    run( 'cat %s/notcovered%i.head %s/notcovered%i.sam | samtools view -bS - > %s/notcovered%i.bam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Convert to bam completed' % stage )

  stage += 1 # 7
  if start <= stage:
    # reads in donor
    run( 'samtools flagstat %s/donor%i.bam > %s/donorflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/reference%i.bam > %s/referenceflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/remapped%i.bam > %s/remappedflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/notcovered%i.bam > %s/notcoveredflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    # coverage
    
    run( 'bedtools genomecov -ibam %s/donor%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/donor%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/reference%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/reference%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/remapped%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/remapped%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Coverage analysis completed' % stage )

  stage += 1 # 8
  if start <= stage:
    run( 'python compare_bams.py --compare_position True --subset_detail True --mismatch_detail 1 --xmfa %s/mauve%i --origin 2 --target 1 %s/donor%i.bam %s/remapped%i.bam > %s/compare_bams%i.log' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    run( 'python extract_mismatched_reads.py --min_distance 50 %s/remapped%i.bam < %s/compare_bams%i.log > %s/mismatched%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    run( 'samtools view -bS %s/mismatched%i.sam > %s/mismatched%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Mismatch analysis completed' % stage )

  stage += 1 # 9
  if start <= stage:
    run( 'samtools flagstat %s/mismatched%i.bam > %s/mismatchedflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/mismatched%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/mismatched%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/notcovered%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/notcovered%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    
    bio.log_stderr( 'Stage %i: Mismatch coverage completed' % stage )

  stage += 1
  if start <= stage:
    print "===== Stats ====="
    # reads
    print "-- Reads --"
    for line in open( '%s/donorflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split(' ')
      if len(fields) > 4 and fields[4] == 'total':
        print "Donor reads total: %s" % fields[0]
        drt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Donor reads mapped: %s" % fields[0]
        drm = int(fields[0])
    print "Donor reads %% mapped: %.1f" % ( 100. * drm / drt )
    for line in open( '%s/referenceflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Reference reads total: %s" % fields[0]
        rrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Reference reads mapped: %s" % fields[0]
        rrm = int(fields[0])
    print "Reference reads %% mapped: %.1f" % ( 100. * rrm / rrt )
    for line in open( '%s/remappedflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Remapped reads total: %s" % fields[0]
        rrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Remapped reads mapped: %s" % fields[0]
        rrm = int(fields[0])
    print "Remapped reads %% mapped: %.1f" % ( 100. * rrm / rrt )
    for line in open( '%s/mismatchedflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Mismatched reads total: %s" % fields[0]
        rrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Mismatched reads mapped: %s" % fields[0]
        rrm = int(fields[0])
    print "Mismatched reads %% mapped: %.1f" % ( 100. * rrm / rrt )
    for line in open( '%s/notcoveredflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Notcovered reads total: %s" % fields[0]
        rrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Notcovered reads mapped: %s" % fields[0]
        rrm = int(fields[0])
    print "Notcovered reads %% mapped: %.1f" % ( 100. * rrm / rrt )

    # coverage
    print "\n-- Coverage --"
    df = open( '%s/donor%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Donor not covered: %s" % df[0]
    print "Donor gaps: %s" % df[5]
    print "Donor max gap: %s" % df[2]

    rf = open( '%s/reference%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Reference not covered: %s" % rf[0]
    print "Reference gaps: %s" % rf[5]
    print "Reference max gap: %s" % rf[2]

    mf = open( '%s/remapped%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Remapped not covered: %s" % mf[0]
    print "Remapped gaps: %s" % mf[5]
    print "Remapped max gap: %s" % mf[2]

    print "\n-- Remapped incorrectly > 50bp --"
    xf = open( '%s/mismatched%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Bases affected by mismatch: %s" % xf[5]
    print "Max mismatch coverage: %s" % xf[2]

    print "\n-- Off target (outside mappable region) --"
    nf = open( '%s/notcovered%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Off target bases: %s" % nf[5]
    print "Max coverage of off target: %s" % nf[2]

    print "\n-- Remapping --"
    for line in open( '%s/remap_bam%i.stats' % ( tmpdir, idx ), 'r' ):
      fields = line.strip().split(':')
      if len(fields) >1 and fields[0] == 'count':
        maplen = int(fields[1])
        print "Mapped bases: %i" % maplen
      if len(fields) >1 and fields[0] == 'blocks':
        print "Mapped blocks: %s" % fields[1].strip()

    print "\n-- Summary --"
    reflen = int( open( '%s/notcovered%i.head' % ( tmpdir, idx ), 'r' ).readline().strip().split(':')[-1] )
    donorlen = int( open( '%s/remapped%i.head' % ( tmpdir, idx ), 'r' ).readline().strip().split(':')[-1] )
    coverage_loss = donorlen - ( reflen - int(rf[0]) )
    print 'Best case loss from reference coverage: %i / %i: %.1f%%' % ( coverage_loss, donorlen, 100. * coverage_loss / donorlen )
    print 'Best case loss from remapping: %i / %i: %.1f%%' % ( donorlen - maplen, donorlen, 100. * ( donorlen - maplen ) / donorlen )
    print 'Loss after remap coverage: %i / %i: %.1f%%' % ( int(mf[0]), donorlen, 100. * int(mf[0]) / donorlen )
    print 'Loss due to remap: %i / %i: %.1f%%' % ( int(mf[0]) - coverage_loss, donorlen, 100. * ( int(mf[0]) - coverage_loss ) / donorlen )
    print 'Potential mismatch impact: %i / %i: %.1f%%' % ( int(xf[5]), donorlen, 100. * int(xf[5]) / donorlen )
    print 'Off target: %i / %i: %.1f%%' % ( int(nf[5]), donorlen, 100. * int(nf[5]) / donorlen )
    print "===== "

    bio.log_stderr( 'Stage %i: Finished' % stage )

  stage += 1
  if start <= stage:
    run( 'rm %s/reference%i.sam %s/mismatched%i.sam %s/remapped%i.sam %s/notcovered%i.sam %s/donor%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    bio.log_stderr( 'Stage %i: Cleanup finished' % stage )
  

