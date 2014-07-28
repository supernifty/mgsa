#############################################
# evaluates a sam that has labelled reads
# optional configuration file on the command line
#############################################

import sys
import bio

READ_LENGTH = 50

if len(sys.argv) > 1:
  config_file = sys.argv[1]
  print "Configuration file: %s" % config_file
  for line in open( config_file, 'r' ):
    if line.startswith( '#' ):
      continue
    print "Configuration Option: %s" % line.strip()

if len(sys.argv) > 2:
  variation_map = open( sys.argv[2], 'r' )
else:
  variation_map = None

ev = bio.SamAccuracyEvaluator( sam=sys.stdin, verbose=True, variation_map=variation_map )
for stat in ev.stats:
  print stat, ev.stats[stat]

for incorrect in ev.incorrect:
  print incorrect

for incorrect in ev.incorrect_diff:
  print 'diff:', incorrect, ' frequency:', ev.incorrect_diff[incorrect]

# calculate some stuff
print "%mapped_unmapped", ( ev.stats['unmapped'] * 100. / ( ev.stats['unmapped'] + ev.stats['incorrect'] + ev.stats['correct'] ) )
print "%mapped_recall", ( ev.stats['correct'] * 100. / ( ev.stats['unmapped'] + ev.stats['incorrect'] + ev.stats['correct'] ) )
if ev.stats['incorrect'] + ev.stats['correct'] > 0:
  print "%mapped_incorrectly_(1-precision)", ( ev.stats['incorrect'] * 100. / ( ev.stats['incorrect'] + ev.stats['correct'] ) )
  print "%clipped", ( ( ev.stats['soft_clipping'] + ev.stats['hard_clipping'] ) * 100. / ( ev.stats['incorrect'] + ev.stats['correct'] ) / READ_LENGTH )
if 'incorrect_S' in ev.stats and 'correct_S' in ev.stats:
  print "%snp_mapped_incorrectly", ( ev.stats['incorrect_S'] * 100. / ( ev.stats['incorrect_S'] + ev.stats['correct_S'] ) )
if 'incorrect_I' in ev.stats and 'correct_I' in ev.stats:
  print "%ins_mapped_incorrectly", ( ev.stats['incorrect_I'] * 100. / ( ev.stats['incorrect_I'] + ev.stats['correct_I'] ) )
if 'incorrect_D' in ev.stats and 'correct_D' in ev.stats:
  print "%del_mapped_incorrectly", ( ev.stats['incorrect_D'] * 100. / ( ev.stats['incorrect_D'] + ev.stats['correct_D'] ) )
