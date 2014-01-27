# -*- coding: utf-8 -*-
#
# test_table.py

import pylab
import numpy
from numpy import ndarray
from datetime import datetime
import sys
import os

spiketrain = numpy.loadtxt("sample.txt")  
print spiketrain.shape, "size of spiketrain"

delta_t = 1
recordtime = 300.
N_record = 50
starttime = 1000.

numbins = numpy.int(numpy.floor(recordtime / delta_t))
print numpy.max(spiketrain[:,0]), "max cell #" # 0th col
print numpy.max(spiketrain[:,1]), "max time" # 1st col
cells = numpy.int(numpy.max(spiketrain[:,0]))
print cells, "number of cells to count from"

# Create binned spike vector
tally = numpy.zeros(N_record)
for cell in range(cells):
 es = spiketrain[numpy.where(spiketrain[:,0] == cell+1),:] # Find entries where cell num eq cell
# print es, "es"
# print es.shape, "es.shape"
 for e in range(es.shape[1]-1): #Assuming es is sorted in time,
   if numpy.floor(es[0,e,1] / delta_t) <> numpy.floor(es[0,e+1,1] / delta_t):
    tally[cell]+=1
   else:
    print es[0,e,1], es[0,e,1], "cell", cell+1
 tally[cell] += 1

#a = numpy.array([])
#for cell in range(cells):
#  a=numpy.append(a,numpy.size(numpy.where(spiketrain[:,0] == cell+1)))
#print "for comparison, pure spike counts without the binning\n", a
#print tally < a
#print tally - a

# Create coincidence matrix
coinc = numpy.zeros((N_record,N_record))
numpy.fill_diagonal(coinc,1) # Count each 
for b in range(numbins):
 mask1 = numpy.array(numpy.where(starttime+delta_t*(b) <= spiketrain[:,1]))
 mask2 = numpy.array(numpy.where(spiketrain[:,1] < starttime+delta_t*(b+1)))
 mask = numpy.intersect1d(mask1.flatten(),mask2.flatten())
 if mask.any(): # If overlap is nonempty, i.e. there were spikes in that bin
  es = spiketrain[mask,:]
  #print mask.shape[0],"spikes in bin", b+1, "(", starttime+delta_t*(b), "to", starttime+delta_t*(b+1), ")"
  # Find all entries within a bin
  for e1 in range(es.shape[1]):
   for e2 in es[e1:]:
    if (numpy.floor(es[e1,0]) <> numpy.floor(e2[0])) or (es[e1,1] <> e2[1]): # If the cell IDs or exact spike times are different
     #print "C", es[e1,0], e2[0], "e1, e2", "bin", b
     coinc[es[e1,0]-1,e2[0]-1]+=1 # index within es starts from 1 
    else:
     print es[e1,0], e2[0], es[e1,1], e2[1], "bin", b
 else:
  print "no spikes in bin", b+1, "(", starttime+delta_t*(b), "to", starttime+delta_t*(b+1), ")"

print "total # spikes:", numpy.sum(tally)
print "total # coincs:", numpy.sum(coinc)
print numpy.mean(coinc)
print numpy.diagonal(coinc)

# Show coincidence matrix
pylab.pcolor(coinc)
pylab.colorbar()
pylab.show()

# Show scatter plot
#pylab.scatter(spiketrain[:,1],spiketrain[:,0])
#pylab.show()
