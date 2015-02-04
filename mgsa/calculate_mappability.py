
import sys
import bio

import config

if len(sys.argv) < 2:
  print "Usage: %s fasta_file read_length variation error mapper mappability_results" % sys.argv[0]
  # e.g. python calculate_mappability.py ../../data/circoviridae_snp_x2.fasta 50 "insert 1" "snp 1" bwasw out/mappability_circoviridae_bwa_sw_insert_1.out
  sys.exit(1)

pos_generator = bio.FastqPosGenerator( open( sys.argv[1], 'r' ) )
#mapper_function = mapper_selector.bowtie_map( sys.argv[1], '--local' )
mapper_function = bio.get_mapper( config, sys.argv[5], sys.argv[1] )
#mapper_function = mapper_selector.bwa_map( sys.argv[1], 'mem' )
#mapper_function = mapper_selector.bwa_map( sys.argv[1], 'bwasw' )
#mapper_function = mapper_selector.soap_map( sys.argv[1] ) # doesn't work
read_length = int(sys.argv[2])

#mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation='snp', min_pos = 50, max_pos = 51, remove_files=False )
#mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation='insert 1' )#, min_pos = 0, max_pos = 1 ) #, min_pos = 50, max_pos = 51 )
#mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation='delete 1' )#, min_pos = 50, max_pos = 51, remove_files=False )
mappability = bio.Mappability( pos_generator, mapper_function, read_length, variation=sys.argv[3], error=sys.argv[4], remove_files=True )#, min_pos=100, max_pos=101, remove_files=False )#, min_pos = 0, max_pos = 1 ) #, min_pos = 50, max_pos = 51 )

with open( sys.argv[6], 'w' ) as output:
  output.write( 'results for %s with read length %i, variation %s, error %s, mapper %s\n' % ( sys.argv[1], read_length, sys.argv[3], sys.argv[4], sys.argv[5] ) )
  output.write( 'accuracy across genome: %s\n' % mappability.accuracy )
  output.write( 'summary: %s\n' % mappability.summary )
