import sys

def depths_to_instances( fh_in, fh_out, depth_line, var_line, sample ):
  for count, line in enumerate( fh_in ):
    if count == depth_line:
      depths = line.split( ',' )
    if count == var_line:
      variations = set( [ int(x) for x in line.split( ',' ) ] )

  print len(depths), "depths", len(variations), "variations"

  # header
  fh_out.write( ','.join( [ 'f%i' % x for x in xrange(-100, 100) ] ) )
  fh_out.write( ',class\n' )

  # instances
  count = sample
  for pos in xrange( 1000, len( depths ) - 1000 ):
    if pos in variations:
      fh_out.write( ','.join( depths[pos - 100: pos + 100 ] ) )
      fh_out.write( ',1\n' )
    else:
      count -=1
      if count == 0:
        count = sample
        fh_out.write( ','.join( depths[pos - 100: pos + 100 ] ) )
        fh_out.write( ',0\n' )

depths_to_instances( sys.stdin, sys.stdout, 3, 4, 110 )

