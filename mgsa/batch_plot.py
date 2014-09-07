
import sys

import matplotlib.pyplot as plt

import helpers

xs, ys = helpers.series_from_pipeline_batch( fh=open(sys.argv[1]), x='insert_prob', y='incorrect' )

print "xs", xs
print "ys", ys
xs = [ x * 100 for x in xs ]
plt.plot( xs, ys, label='BWA' )
plt.ylabel('Incorrect %')
plt.xlabel('Mutation %')
plt.legend()
plt.show()
