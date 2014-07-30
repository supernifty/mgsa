
import datetime
import random
import re
import sys

def log_stderr(msg):
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  sys.stderr.write( '%s: %s\n' % ( when, msg ) )

class Config(object):
  '''
    generates configuration dictionary from a line config spec
    e.g. snp_prob 0.03, insert_prob 0, delete_prob 0, fasta hiv, mapper baseline, read_length 50
  '''
  floats = ( 'insert_prob', 'delete_prob', 'snp_prob', 'error_prob' )
  ints = ( 'max_insertion_len', 'max_deletion_len', 'read_length', 'coverage', 'mult')

  def __init__( self ):
    pass

  def read_config_line( self, line ):
    cfg = { 'insert_prob': 0, 'delete_prob': 0, 'snp_prob': 0, 'max_insertion_len': 1, 'max_deletion_len': 1, 'mult': 1, 'mapper': 'bowtie2', 'mutation_type': '', 'read_length': 50, 'fasta': 'circoviridae', 'coverage': 10, 'error_prob': 0 }
    for option in line.split( ',' ):
      key, value = option.strip().split()
      if key in Config.floats:
        cfg[key] = float(value)
      elif key in Config.ints:
        cfg[key] = int(value)
      else:
        cfg[key] = value
      # generate brief description
      if key in ( 'insert_prob', 'delete_prob', 'snp_prob' ):
        cfg['mutation_type'] = cfg['mutation_type'] + '_' + key[:3]
      if key in ( 'max_insertion_len', 'max_deletion_len' ):
        cfg['mutation_type'] = cfg['mutation_type'] + '_' + value
    return cfg

  def read_config_file( self, fh ):
    cfg = { 'insert_prob': 0, 'delete_prob': 0, 'snp_prob': 0, 'max_insertion_len': 1, 'max_deletion_len': 1, 'mult': 1, 'mapper': 'bowtie2', 'mutation_type': '', 'read_length': 50, 'fasta': 'circoviridae', 'coverage': 10, 'error_prob': 0 }
    for line in fh:
      if line.startswith( '#' ):
        continue
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
