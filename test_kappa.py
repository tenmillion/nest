# -*- coding: utf-8 -*-
#
# test_kappa.py

import pylab
import numpy
from numpy import ndarray

s1 = numpy.array([100+n*20 for n in range(20)])
s2 = numpy.array([105+n*20 for n in range(19)])
s2 = numpy.append(s2, 105+19*20-4)
spikes = [s1,s2]
print spikes

binwidth = 2.
xymat = numpy.zeros((len(spikes),len(spikes)))
id1 = 0
for st1 in spikes:
	id2 = 0
	for st2 in spikes:
#		intervals = numpy.array([])
		overlaps = 0
		for k in xrange(0, st1.shape[0]):
#			intervals = numpy.append(intervals, st1[k] - st2[numpy.nonzero((st1[k] - binwidth/2. < st2) & (st2 <= st1[k] + binwidth/2.))])
			overlaps += numpy.count_nonzero((st1[k] - binwidth/2. <= st2) & (st2 <= st1[k] + binwidth/2.))
#		xy = numpy.append(xy, intervals.shape[0])
		print overlaps, "overlaps"
		print "spike trains:", id1, "("+str(st1.shape[0])+" spikes) and", id2, "("+str(st2.shape[0])+" spikes)"
		xymat[id1][id2] = overlaps
		id2 +=1
	id1 +=1

print xymat