
import re
import itertools

import bio

# parses mauve output files

class MauveMap( object ):
  def __init__( self, fh, src_strand=1, target_strand=2, new_reference=None, log=bio.log_stderr ):
    '''
      @fh: xmfa file like object of Mauve output
    '''
    self.log = log
    self.coverage = {}
    self.target_map = None
    self.genome_stats = { 'xmin': 1e9, 'xmax': 0, 'ymin': 1e9, 'ymax': 0, 'blocks': 0 }
    self.new_reference = new_reference
    current_sequence = [ '', '' ] # src, target
    src_range = [ 0, 0 ] # src, target
    target_range = [ 0, 0 ] # src, target
    current = -1
    direction = [ '+', '+' ]
    xmax = 0
    ymax = 0
    for line_num, line in enumerate(fh):
      line = line.strip()
      #print "processing", line
      if line.startswith( '>' ): # e.g. > 1:213657-501434 - e-coli-mg1655.fasta
        items = line.split( ' ' )
        strand, sequence = items[1].split( ':' )
        strand = int(strand)
        if strand == src_strand:
          src_range = [ int(x) for x in sequence.split( '-' ) ]
          current = 0
          current_sequence[current] = ''
          direction[current] = items[2] # + or -
        elif strand == target_strand:
          target_range = [ int(x) for x in sequence.split( '-' ) ]
          current = 1
          current_sequence[current] = ''
          direction[current] = items[2] # + or -
        else:
          current = -1
      elif line == '=':
        self.genome_stats['blocks'] += 1
        if len( current_sequence[0] ) == len( current_sequence[1] ):
          self.log( "adding coverage for %s -> %s" % ( src_range, target_range ) )
          target_pos = 0 
          src_pos = 0 
          added = 0
          for forward_pos, base in enumerate( current_sequence[0] ): # src; forward_pos = 0..len
            reverse_pos = len( current_sequence[0] ) - forward_pos - 1
            src_base = self._current_base( current_sequence[0], forward_pos, reverse_pos, direction[0] )
            target_base = self._current_base( current_sequence[1], forward_pos, reverse_pos, direction[1] )
            if src_base != '-' and target_base != '-': # no mapping
              actual_src_pos = self._get_actual_pos( src_pos, src_range, direction[0] )
              actual_target_pos = self._get_actual_pos( target_pos, target_range, direction[1] )
              self.coverage[ actual_src_pos ] = actual_target_pos
              added += 1
              self.genome_stats['xmax'] = max( self.genome_stats['xmax'], actual_src_pos )
              self.genome_stats['ymax'] = max( self.genome_stats['ymax'], actual_target_pos )
              self.genome_stats['xmin'] = min( self.genome_stats['xmin'], actual_src_pos )
              self.genome_stats['ymin'] = min( self.genome_stats['ymin'], actual_target_pos )
            if src_base != '-':
              src_pos += 1
            if target_base != '-':
              target_pos += 1
          self.log( "added %i mappings" % added )
        else:
          self.log( 'skipping unmatched sequence of length %i at %i' % ( len( current_sequence[0] ), src_range[0] ) )
        current_sequence = [ '', '' ] # src, target
        direction = [ '+', '+' ]
      else: # sequence
        #print current, "sequence", line
        if current >= 0:
          current_sequence[current] += line
      if line_num < 5 or line_num % 10000 == 0:
        self.log( '%i lines processed; %i mappings; range: %s' % ( line_num, len(self.coverage), self.genome_stats ) )
    self.genome_stats['count'] = len(self.coverage)
    self.log( '%i lines processed; stats: %s' % ( line_num, self.genome_stats ) )

  def _get_actual_pos( self, src_pos, src_range, direction ):
    if direction == '+': # normal
      return src_pos + src_range[0]
    else: # reverse
      return src_range[1] - src_pos

  def _current_base( self, sequence, forward_pos, reverse_pos, direction ):
    if direction == '+':
      return sequence[forward_pos]
    else:
      return sequence[reverse_pos]

  def find_nearest_target( self, target_base, max_dist=100000 ):
    if self.target_map is None:
      self.target_map = set()
      for source, dest in self.coverage.items():
        self.target_map.add( dest )
        
    diff = 0
    while diff < max_dist:
      if target_base + diff in self.target_map:
        return diff
      if target_base - diff in self.target_map:
        return -diff
      diff += 1
    return None
 
  def find_nearest( self, src_base, max_dist=100000 ):
    diff = 0
    while diff < max_dist:
      if src_base + diff in self.coverage:
        return diff
      if src_base - diff in self.coverage:
        return -diff
      diff += 1
    return None
    
  def remap( self, sam_fh, output, remap_cigar=False, not_covered_output=None ):
    '''
      writes a sam file with the locations remapped
    '''
    self.stats = { 'total': 0, 'unmapped': 0, 'mapped': 0, 'reads_covered': 0, 'reads_notcovered': 0, 'reads_partial': 0 }
    for pos, line in enumerate(sam_fh):
      line = line.strip()
      if self.new_reference is not None and line.startswith( '@SQ' ) and line.find( 'SN:') != -1:
        line = re.sub(r"SN:[^\s]*", "SN:%s" % self.new_reference, line )
        line = re.sub(r"LN:[^\s]*", "LN:%i" % self.genome_stats['ymax'], line ) # also sub in ymax
        output.write( '%s\n' % line )
        continue
      if line.startswith('@'):
        output.write( '%s\n' % line )
        continue
      fields = line.split()
      if len(fields) > 9:
        #print fields
        if self.new_reference is not None:
          fields[2] = self.new_reference
        flag = int( fields[1] )
        self.stats['total'] += 1
        if flag & 0x04 != 0:
          self.stats['unmapped'] += 1
          #output.write( '%s\n' % line ) # write unmapped reads verbatim
          fields[2] = '*' # make sure it's unmapped!
          output.write( '\t'.join( fields ) )
          output.write( '\n' )
        else:
          self.stats['mapped'] += 1
          # assess coverage
          start_pos = int(fields[3])
          read_len = len(fields[9]) # TODO not always correct, not considering cigar
          if all( [ x in self.coverage for x in xrange( start_pos, start_pos + read_len ) ] ):
            self.stats['reads_covered'] += 1 # read is completely mappable, TODO fully mapped is not necessarily contiguous, and what about deletions?
            #self.log( 'read at %i -> %i covered: %s' % ( start_pos, start_pos + read_len, start_pos in self.coverage ) )
            # convert to target position
            if remap_cigar:
              fields[5] = self.calculate_cigar( int(fields[3]), len(fields[9]) )
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
            # write if desiredd
            if not_covered_output is not None:
                not_covered_output.write( '%s\n' % line ) # write unmapped reads verbatim
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
 
  @staticmethod
  def run_length_encode( cigar ):
    return ''.join( [ '%i%s' % (len(list(g)), k) for k,g in itertools.groupby(cigar)] )

  def calculate_cigar( self, source_pos, read_len ):
    '''
      look at the mauve mapping of each individual base to come up with a new cigar string for the read
    '''
    source_end_pos = source_pos + read_len
    cigar = []
    #self.log( 'looking for coverage at %i' % source_pos )
    expected_target_pos = self.coverage[ source_pos ] # assume 1st position covered
    while source_pos < source_end_pos:
      if source_pos in self.coverage:
        actual_target_pos = self.coverage[ source_pos ] # where it maps to
        diff = expected_target_pos - actual_target_pos
        if diff == 0: # match
          cigar.append( 'M' )
          expected_target_pos += 1
        elif diff > 0 and diff < source_end_pos - source_pos: # it's moved along the read (deletion)
            cigar.extend( 'D' * diff )
            expected_target_pos += diff + 1
        else: # it's moved off the read (insertion)
            cigar.append( 'I' )
            # keep on same target pos
      else: # not in coverage (insertion)
        cigar.append( 'I' )
        # keep on same target pos
      source_pos += 1

    # compress into cigar
    return MauveMap.run_length_encode( cigar )
