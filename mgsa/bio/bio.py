
import datetime
import random
import re
import sys

def log_stderr(msg):
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  sys.stderr.write( '%s: %s\n' % ( when, msg ) )

def log_quiet(msg):
  pass

complement = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N' }
def reverse_complement(s):
  return ''.join( map( lambda n: complement[n], s )[::-1] )

class Config(object):
  '''
    generates configuration dictionary from a line config spec
    e.g. snp_prob 0.03, insert_prob 0, delete_prob 0, fasta hiv, mapper baseline, read_length 50
    description of parameters:
    * chromosomes: true|false => true to use multi chromosome processors
    * insert_prob: probability of an insertion at any base in the sequence
    * delete_prob: probability of a deletion at any base in the sequence
    * snp_prob: probability of a mutation at any base in the sequence
    * min_insertion_len: if an insertion added, it will be at least this length
    * max_insertion_len: if an insertion added, it will be at most this length
    * min_deletion_len: if a deletion added, at least this many bases deleted
    * max_deletion_len: if a deletion added, at most this many bases deleted
    * min_variation_distance: minimum distance between variations
    * min_variation_start: minimum start for any variation
    * mult
    * mult_snp_prob
    * mapper
    * mutation_type
    * read_length
    * fasta
    * coverage
    * error_prob: prob of error at any position
    * error_type: e.g. snv, homopolymer 3
    * command: fasta, mutation, reads, vcf: only generate the vcf, don't evaluate
    * inversion_prob
    * paired_end_length
    * paired_end_sd
    * reports
    * bias_report_buckets
    * insertion_source: random (not implemented) random = randomly generated bases according to error profile, repeat = take source from genome
    * call_strategy: majority (not implemented) specifies the strategy for predicting a variation. majority = most reads predict the variation
    * sam_source: generated or filename; generated means the sam will be generated by running a mapping job with the aligner specified by mapper
    * fasta_reference: generated or filename
    * vcf_exclude: filename for existing variations that should not be considered
    * vcf_out: filename for recovered vcf
    * vcf_source: generated or filename; generated means the vcf has been created by an earlier "mutation" command. filename means there's an existing vcf of correct mutations
  '''
  floats = ( 'insert_prob', 'delete_prob', 'snp_prob', 'error_prob', 'inversion_prob', 'mult_snp_prob' )
  ints = ( 'min_insertion_len', 'max_insertion_len', 'min_deletion_len', 'max_deletion_len', 'min_variation_distance', 'min_variation_start', 'read_length', 'coverage', 'mult', 'bias_report_buckets')
  # insertion_source: novel (from another genome), repeated (from this genome), random (generated)
  default_config = { 'insert_prob': 0, 'delete_prob': 0, 'snp_prob': 0, 'min_insertion_len': 1, 'max_insertion_len': 1, 'min_deletion_len': 1, 'max_deletion_len': 1, 'min_variation_distance': 1, 'min_variation_start': 0, 'mult': 1, 'mult_snp_prob': 0, 'mapper': 'bowtie2', 'mutation_type': '', 'read_length': 100, 'fasta': 'circoviridae', 'coverage': 10, 'error_prob': 0, 'error_type': 'snv', 'command': '', 'inversion_prob': 0, 'paired_end_length': 0, 'paired_end_sd': 0, 'reports': '-', 'bias_report_buckets': 10, 'insertion_source': 'random', 'call_strategy': 'majority', 'sam_source': 'generated', 'fasta_reference': 'generated', 'vcf_source': 'generated', 'vcf_exclude': 'none', 'chromosomes': 'false', 'vcf_out': 'generated' }

  def __init__( self ):
    pass

  def read_config_line( self, line ):
    # defaults
    cfg = Config.default_config.copy()
    for option in line.split( ',' ):
      key, value = option.strip().split()
      if key in Config.floats:
        cfg[key] = float(value)
      elif key in Config.ints:
        cfg[key] = int(value)
      else:
        cfg[key] = value
      # generate brief description for filename
      if key in ( 'insert_prob', 'delete_prob', 'snp_prob' ):
        cfg['mutation_type'] = cfg['mutation_type'] + '_' + key[:3]
      if key in ( 'max_insertion_len', 'max_deletion_len' ):
        cfg['mutation_type'] = cfg['mutation_type'] + '_' + value
    return cfg

  def read_config_file( self, fh ):
    # defaults
    cfg = Config.default_config.copy()
    for line in fh:
      if line.startswith( '#' ):
        continue
      #print "processing config line", line
      key, value = line.strip().split()
      if key in Config.floats:
        cfg[key] = float(value)
      elif key in Config.ints:
        cfg[key] = int(value)
      else:
        cfg[key] = value
    return cfg

class ErrorGenerator(object):
  def __init__( self, error_profile ):
    self.error_profile = error_profile

  def apply_errors( self, dna ):
    result = []
    for x in dna:
      result.append( self.error_profile( x ) )
    return ''.join( result )

  @staticmethod
  def create_uniform_error_profile( error_prob ):
    transitions = { 'A': 'TGC', 'T': 'GCA', 'G': 'ACT', 'C': 'AGT', 'N': 'N' }
    def uniform_error_profile( bp ):
      if random.random() < error_prob:
        return transitions[bp][random.randint(0, len(transitions[bp]) - 1)]
      else:
        return bp
    return uniform_error_profile

if __name__ == "__main__":
  #s = SamToFasta( sys.stdin, log )
  #log( 'size is %i' % s.fasta.size )

  #log( s.fasta.prior )
  #for i in xrange(1, 10):
  #  print i, s.fasta.count(i), s.fasta.confidence(i)

  #pos = 1
  #while pos < s.fasta.size:
  #  print s.fasta.consensus( pos, pos + 70 )
  #  pos += 70

  #d = FastaDiff( FastaReader( open( '../data/test/test1k.fasta', 'r' ) ), s.fasta, log )
  #log( d.error_total )
  #log( d.errors )

  m = FastaMutate( FastaReader( open( '../data/test/test1k.fasta', 'r' ) ), log, open( '../data/test/test1ksnp.vcf', 'w' ), 0.1 )
  o = open( '../data/test/test1ksnp.fasta', 'w' )
  o.write( '>mutated-snp-test1k\n')
  for line in m.items():
    o.write( '%s\n' % line )
  o.close()
  log( m.mutations )
