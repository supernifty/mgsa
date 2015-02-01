
import sys
import bio

import mapper_selector

if len(sys.argv) < 2:
  print "Usage: %s fasta_file read_length mappability_results" % sys.argv[0]
  sys.exit(1)

pos_generator = bio.FastqPosGenerator( open( sys.argv[1], 'r' ) )
#mapper_function = mapper_selector.bowtie_map( sys.argv[1], '--local' )
#mapper_function = mapper_selector.bwa_map( sys.argv[1], 'mem' )
mapper_function = mapper_selector.bwa_map( sys.argv[1], 'bwasw' )
#mapper_function = mapper_selector.soap_map( sys.argv[1] ) # doesn't work
read_length = int(sys.argv[2])

#mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation='snp', min_pos = 50, max_pos = 51, remove_files=False )
mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation='insert 1' )#, min_pos = 0, max_pos = 1 ) #, min_pos = 50, max_pos = 51 )
#mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation='delete 1' )#, min_pos = 50, max_pos = 51, remove_files=False )

with open( sys.argv[3], 'w' ) as output:
  output.write( 'results for %s with read length %i\n' % ( sys.argv[1], read_length ) )
  output.write( 'accuracy across genome: %s\n' % mappability.accuracy )
  output.write( 'summary: %s\n' % mappability.summary )
