
import bio

# parses mauve output files

class MauveMap( object ):
  def __init__( self, fh, src_strand=1, target_strand=2, log=bio.log_stderr ):
    '''
      @fh: xmfa file like object of Mauve output
    '''
    self.log = log
    self.coverage = {}
    current_sequence = [ '', '' ] # src, target
    src_range = [ 0, 0 ] # src, target
    target_range = [ 0, 0 ] # src, target
    current = -1
    xmax = 0
    ymax = 0
    for line_num, line in enumerate(fh):
      line = line.strip()
      #print "processing", line
      if line.startswith( '>' ):
        items = line.split( ' ' )
        strand, sequence = items[1].split( ':' )
        strand = int(strand)
        if strand == src_strand:
          src_range = [ int(x) for x in sequence.split( '-' ) ]
          current = 0
          current_sequence[current] = ''
        elif strand == target_strand:
          target_range = [ int(x) for x in sequence.split( '-' ) ]
          current = 1
          current_sequence[current] = ''
        else:
          current = -1
      elif line == '=':
        target_pos = 0
        src_pos = 0
        if len( current_sequence[0] ) == len( current_sequence[1] ):
          self.log( "adding coverage for %s -> %s" % ( src_range, target_range ) )
          added = 0
          for pos, base in enumerate( current_sequence[0] ): # src
            src_base = current_sequence[0][pos]
            target_base = current_sequence[1][pos]
            if src_base != '-' and target_base != '-': # no mapping
              self.coverage[ src_pos + src_range[0] ] = target_pos + target_range[0]
              added += 1
              xmax = max( xmax, src_pos + src_range[0] )
              ymax = max( ymax, target_pos + target_range[0] )
            if src_base != '-':
              src_pos += 1
            if target_base != '-':
              target_pos += 1
          self.log( "added %i mappings" % added )
        else:
          self.log( 'skipping unmatched sequence of length %i at %i' % ( len( current_sequence[0] ), src_range[0] ) )
        current_sequence = [ '', '' ] # src, target
      else: # sequence
        #print current, "sequence", line
        if current >= 0:
          current_sequence[current] += line
      if line_num < 5 or line_num % 10000 == 0:
        self.log( '%i lines processed; %i mappings; src max: %i; target max: %i' % ( line_num, len(self.coverage), xmax, ymax ) )
    self.log( '%i lines processed; %i mappings; src max: %i; target max: %i' % ( line_num, len(self.coverage), xmax, ymax ) )

  def remap( self, sam_fh, output ):
    '''
      writes a sam file with the locations remapped
    '''
    self.stats = { 'total': 0, 'unmapped': 0, 'mapped': 0, 'reads_covered': 0, 'reads_notcovered': 0, 'reads_partial': 0 }
    for pos, line in enumerate(sam_fh):
      line = line.strip()
      fields = line.split()
      if len(fields) > 9:
        #print fields
        flag = int( fields[1] )
        self.stats['total'] += 1
        if flag & 0x04 != 0:
          self.stats['unmapped'] += 1
          output.write( '%s\n' % line ) # write unmapped reads verbatim
        else:
          self.stats['mapped'] += 1
          # assess coverage
          start_pos = int(fields[3])
          read_len = len(fields[9]) # TODO not always correct, not considering cigar
          if all( [ x in self.coverage for x in xrange( start_pos, start_pos + read_len ) ] ):
            self.stats['reads_covered'] += 1 # read is completely mappable, TODO fully mapped is not necessarily contiguous
            # convert to target position
            fields[3] = str( self.coverage[ start_pos ] )
            output.write( '\t'.join( fields ) )
            output.write( '\n' )
          elif any( [ x in self.coverage for x in xrange( start_pos, start_pos + read_len ) ] ):
            self.stats['reads_partial'] += 1 # read is partially mappable
            # convert to unmapped read TODO could still map this
            fields[1] = str( flag & 0x04 )
            fields[2] = '*'
            fields[3] = '0'
            fields[4] = '0'
            fields[5] = '*'
            fields[6] = '*'
            output.write( '\t'.join( fields ) )
            output.write( '\n' )
          else:
            self.stats['reads_notcovered'] += 1 # no part of the read is mappable
            # convert to unmapped read
            fields[1] = str( flag & 0x04 )
            fields[2] = '*'
            fields[3] = '0'
            fields[4] = '0'
            fields[5] = '*'
            fields[6] = '*'
            output.write( '\t'.join( fields ) )
            output.write( '\n' )
      else: # write non-reads verbatim
        output.write( '%s\n' % line )
      if pos < 10 or pos % 100000 == 0:
        self.log( "processed %i lines: %s..." % ( pos, self.stats ) )
 
