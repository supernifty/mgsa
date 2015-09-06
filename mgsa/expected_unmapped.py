
import numpy

import bio

print "bowtie2"
for e in numpy.arange(0, 0.11, 0.01):
  print e, bio.expected_unmapped( e, 20, 100 )

print "bwa"
for e in numpy.arange(0, 0.11, 0.01):
  print e, bio.expected_unmapped( e, 19, 100, segments=11 )

print "perfect"
for e in numpy.arange(0, 0.11, 0.01):
  print bio.run_probability( e, 20, 100 )
#for e in numpy.arange(0, 0.11, 0.01):
#  print e, bio.expected_unmapped( e, 19, 100, segments=11 )

print "e,c,p"
print "0.01,1,", bio.false_positive_rate(0.01, 0.01, 1), 100. / bio.false_positive_rate(0.01, 0.01, 1)
print "0.02,1,", bio.false_positive_rate(0.02, 0.01, 1), 100. / bio.false_positive_rate(0.02, 0.01, 1)
print "0.10,1,", bio.false_positive_rate(0.10, 0.01, 1), 100. / bio.false_positive_rate(0.10, 0.01, 1)
print "0.01,2,", bio.false_positive_rate(0.01, 0.01, 2), 100. / bio.false_positive_rate(0.01, 0.01, 2)
print "0.02,2,", bio.false_positive_rate(0.02, 0.01, 2), 100. / bio.false_positive_rate(0.02, 0.01, 2)
print "0.10,2,", bio.false_positive_rate(0.10, 0.01, 2), 100. / bio.false_positive_rate(0.10, 0.01, 2)
print "0.01,3,", bio.false_positive_rate(0.01, 0.01, 3), 100. / bio.false_positive_rate(0.01, 0.01, 3)
print "0.02,3,", bio.false_positive_rate(0.02, 0.01, 3), 100. / bio.false_positive_rate(0.02, 0.01, 3)
print "0.10,3,", bio.false_positive_rate(0.10, 0.01, 3), 100. / bio.false_positive_rate(0.10, 0.01, 3)
print "0.01,10,", bio.false_positive_rate(0.01, 0.01, 10), 100. / bio.false_positive_rate(0.01, 0.01, 10)
print "0.02,10,", bio.false_positive_rate(0.02, 0.01, 10), 100. / bio.false_positive_rate(0.02, 0.01, 10)
print "0.10,10,", bio.false_positive_rate(0.10, 0.01, 10), 100. / bio.false_positive_rate(0.10, 0.01, 10)
