
import sys

import bio

ev = bio.SamAccuracyEvaluator( sam=sys.stdin )
for stat in ev.stats:
  print stat, ev.stats[stat]

for incorrect in ev.incorrect:
  print incorrect

for incorrect in ev.incorrect_diff:
  print 'diff:', incorrect, ' frequency:', ev.incorrect_diff[incorrect]

# calculate some stuff
print "%mapped incorrectly", ( ev.stats['incorrect'] * 100. / ( ev.stats['incorrect'] + ev.stats['correct'] ) )
if 'incorrect_S' in ev.stats:
  print "%snp mapped incorrectly", ( ev.stats['incorrect_S'] * 100. / ( ev.stats['incorrect_S'] + ev.stats['correct_S'] ) )
if 'incorrect_I' in ev.stats:
  print "%variation mapped incorrectly", ( ev.stats['incorrect_I'] * 100. / ( ev.stats['incorrect_I'] + ev.stats['correct_I'] ) )
if 'incorrect_D' in ev.stats:
  print "%variation mapped incorrectly", ( ev.stats['incorrect_D'] * 100. / ( ev.stats['incorrect_D'] + ev.stats['correct_D'] ) )
