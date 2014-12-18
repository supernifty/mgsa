
import os
import sys

import config

def velvet_assemble( fastq, directory ):
  os.system( '%sh %s 31 -short -fastq %s' % ( config.VELVET_PATH, directory, fastq ) )
  os.system( '%sg %s' % ( config.VELVET_PATH, directory ) )

if len( sys.argv ) != 4:
  print "Usage: %s velvet fastq directory" % sys.argv[0]
  sys.exit()

if sys.argv[1] == 'velvet':
  velvet_assemble( sys.argv[2], sys.argv[3] )
else:
  print "Unrecognized assembler: %s" % sys.argv[1]
