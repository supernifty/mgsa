#
# build a samfile with bowtie or other aligner
#
import os
import sys

def get_mapper( config, name, fasta_file, seed_length=None ):
  if name == 'baseline':
    return baseline_map( config, fasta_file )
  elif name == 'bowtie2':
    if seed_length is None or seed_length == '0':
      return bowtie_map( config, fasta_file, '--local' )
    else:
      return bowtie_map( config, fasta_file, '--local -D 15 -R 2 -N 0 -L %s -i S,1,0.75' % seed_length )
  elif name == 'bowtie2_exhaustive':
    return bowtie_map( config, fasta_file, '--local -D 1000 -R 2 -N 0 -L 20 -i S,1,0' )
  elif name == 'bowtie2_e2e':
    return bowtie_map( config, fasta_file, '--end-to-end' )
  elif name == 'bwa':
    if seed_length is None or seed_length == '0':
      return bwa_map( config, fasta_file, 'mem' )
    else:
      return bwa_map( config, fasta_file, 'mem -k %s' % seed_length )
  elif name == 'bwa_indel':
    return bwa_map( config, fasta_file, 'mem -w 1000' )
  elif name == 'bwasw':
    return bwa_map( config, fasta_file, 'bwasw' )
  elif name == 'soap':
    return soap_map( config, fasta_file )
  else:
    raise Exception( "Unrecognized mapper %s" % name )
 
def baseline_map( config, fasta_file, params=None ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  def result( fastq_file, sam_file ):
    os.system( 'python mapper.py %s %s %s' % ( fasta_file, fastq_file, sam_file ) )
  return result

def bowtie_map( config, fasta_file, params ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  def result( fastq_file, sam_file, index=True ):
    if index:
      os.system( '%s-build %s %s-bt2' % ( config.BOWTIE_PATH, fasta_file, fasta_file ) )
    os.system( '%s %s -p 16 -x %s-bt2 -U %s --quiet -S %s' % ( config.BOWTIE_PATH, params, fasta_file, fastq_file, sam_file ) ) # -t adds time
  return result

def bwa_map( config, fasta_file, params ):
  '''build candidate.sam from candidate.fastq and reference.fasta'''
  # evaluate bowtie - build bowtie sam
  def result( fastq_file, sam_file, index=True ):
    if index:
      os.system( '%s index -p %s-bwa %s' % ( config.BWA_PATH, fasta_file, fasta_file ) )
    os.system( '%s %s %s-bwa %s > %s' % ( config.BWA_PATH, params, fasta_file, fastq_file, sam_file ) ) # -v 2 for mem
  return result

def soap_map( config, fasta_file ):
  def result( fastq_file, sam_file, index=True ):
    if index:
      os.system( '%s/2bwt-builder %s' % ( config.SOAP_PATH, fasta_file ) )
    os.system( '%s/soap -a %s -D %s.index -o %s' % ( config.SOAP_PATH, fastq_file, fasta_file, sam_file ) )
  return result

