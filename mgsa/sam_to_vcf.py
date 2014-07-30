
import bio
import os
import sys

if len(sys.argv) != 3:
  print "Usage: %s ref bam" % sys.argv[0]
  sys.exit(0)

# ref sam -> vcf

def run( cmd ):
  bio.log_stderr( cmd )
  os.system( cmd )

reference_file = sys.argv[1]
sam_file = sys.argv[2]

# sam -> bam
#run( "samtools view -Sb %s > %s.bam" % ( sam_file, sam_file ) ) 
#run( "samtools sort %s.bam %s_sorted" % ( sam_file, sam_file ) ) 

# bam, ref -> vcf
run( "samtools mpileup -ugf %s %s_sorted.bam | bcftools view -bvcg - | bcftools view - | vcfutils.pl varFilter" % ( reference_file, sam_file ) )
#run( "samtools faidx %s" % reference_file )
#run( "samtools mpileup -uf %s.fai %s_sorted.bam | bcftools view -" % ( reference_file, sam_file ) )
