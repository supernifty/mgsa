
import bio
import fasta
# before error_bias_no_variation
#y_map_old = { 'unmapped': -20, 'incorrect': -19, 'read_precision': -18, 'read_recall': -17, 'read_f1': -16, 'vcf_tp': -15, 'vcf_fp': -14, 'vcf_fn': -13, 'vcf_precision': -12, 'vcf_recall': -11, 'vcf_f1': -10, 'vcf_bucket_tp': -9, 'vcf_bucket_fp': -8, 'vcf_bucket_fn': -7, 'reference_bias': -6, 'error_bias': -5, 'unmapped_variations': -4, 'total_variations': -3, 'mean_reference': -2, 'mean_error': -1 }

def fasta_stats( fh, log=bio.log_stderr ):
  return fasta.FastaStats( fh, log )

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
  return line
  
def fasta_edit( fh, out, commands, width=50 ):
  '''
    command is comma separated list of commands of the form [i|d].pos.[bases|count]
    e.g. i.50.A
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
  
