# -*- coding: utf-8 -*-
#
# test_cc.py

import pylab
import numpy
from numpy import ndarray
import sys
from NeuroTools import analysis, signals

cells = 50 # number of cells recorded from
starttime = 10.
recordtime = 300.
delta_t = 2.0

if len(sys.argv) == 1: # test mode
 print "test mode"
 times = numpy.arange(10,301,10)
 ss = (times.shape)[0]
 test = numpy.zeros((2,cells*ss))
 print ss
 for n in range(cells):
  for s in range(ss):
   print n*cells+s
   test[:,s*cells+n] = [n,times[s]]git 
 numpy.savetxt("./test.txt",numpy.transpose(test))
 spiketrain = numpy.loadtxt("./test.txt")

else:
 spiketrain = numpy.loadtxt(sys.argv[1])

print spiketrain.shape, "size of spiketrain"
pylab.scatter(spiketrain[:,1],spiketrain[:,0])
pylab.show()

print numpy.max(spiketrain[:,0]), "max cell #" # 0th col
print numpy.max(spiketrain[:,1]), "max time" # 1st col
print "number of cells to count from:", cells

#Spike Separation
indiv_trains = []
for cell in range(cells):
  locs = (numpy.where(spiketrain[:,0] == cell+1))[0].tolist() #where returns tuple
  indiv_trains.append(spiketrain[locs].tolist())
print indiv_trains

print "crosscorrelation"
#Pairwise kappa
ccmat = numpy.zeros((cells,cells))
for cell1 in range(cells):
 for cell2 in range(cells):
  #print numpy.array(indiv_trains[cell1])
  st1 = signals.SpikeTrain(numpy.array(indiv_trains[cell1]).flatten())
  st2 = signals.SpikeTrain(numpy.array(indiv_trains[cell2]).flatten())
  #print st1, st2
  #analysis.crosscorrelate(st1,st2,display=True,lag=1,kwargs={'bins':1000}) #plot requires display=True
  #pylab.show()
  out=analysis.crosscorrelate(st1,st2,lag=delta_t/2.,display=False) #return requires display=False
  #print out[0]
  ccmat[cell1,cell2] = (out[0].shape)[0]

#Check that it performed correctly
cell1=9
cell2=9
print len(indiv_trains[1]), len(indiv_trains[cell2]), ccmat[cell1,cell2], cell1, cell2

#Get kappa_ii
kappa_ii = ccmat.diagonal()
print kappa_ii

kappamat = numpy.zeros((cells,cells))
for cell1 in range(cells):
 for cell2 in range(cells):
  kappamat[cell1,cell2] = ccmat[cell1,cell2]/numpy.sqrt(kappa_ii[cell1]*kappa_ii[cell2])

pylab.pcolor(ccmat)
pylab.colorbar()
pylab.show()

pylab.pcolor(kappamat)
pylab.colorbar()
pylab.show()

print numpy.max(ccmat), numpy.min(ccmat)
print numpy.max(kappamat), numpy.min(kappamat)
print numpy.mean(kappamat)
