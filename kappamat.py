# -*- coding: utf-8 -*-
# kappamat.py
# Function to make kappa matrix
# Yoriko Yamamura Jan 29 2014

# Usage:
# from kappamat import *
# s = numpy.loadtxt("test.txt")
# ss=sort_spikes(s)
# kmat=kappamat(ss,2.)
# Or if you just want the average kappa,
# k = kappa(s,2.)

import numpy

def sort_spikes(spikes):
	spikes_tr=numpy.transpose(spikes)
	mincell=numpy.int(min(spikes_tr[0]))
	maxcell=numpy.int(max(spikes_tr[0])+1)
	sorted_spiketrains = []
	for cell in range(mincell,maxcell):
	  locs = (numpy.where(spikes_tr[0] == cell+1))[0] #where returns tuple
	  sorted_spiketrains.append(spikes_tr[1,locs])
	return sorted_spiketrains

def pairwise_kappa(s1, s2, binwidth): # expects numpy arrays
	overlaps = 0
	for k in xrange(0, s1.shape[0]):
		overlaps += numpy.count_nonzero((s1[k] - binwidth/2. <= s2) & (s2 <= s1[k] + binwidth/2.))
#		print overlaps, "overlaps"
#		print "spike trains:", id1, "("+str(st1.shape[0])+" spikes) and", id2, "("+str(st2.shape[0])+" spikes)"
	return overlaps
	
# Somehow when I import functions the array[i][j] syntax doesn't work even for numpy arrays
def kappamat(sorted, binwidth):
	xymat = []
	for s1 in range(len(sorted)):
		row = []
		st1=sorted[s1]
		for s2 in range(len(sorted)):
			st2=sorted[s2]
			if ((len(st1)!=0) and (len(st2)!=0)):
				row.append(pairwise_kappa(st1, st2, binwidth)) # not normalized
			else:
				row.append(0.)
		xymat.append(row)
	return xymat
	
def normalize_kappa(kmat):
	normalized=numpy.zeros(numpy.shape(kmat))
	kmatnp=numpy.array(kmat)
	dia=kmatnp.diagonal()
	for row in range(len(kmat)):
		for elem in range(len(kmat)):
			if ((dia[row]!=0) and (dia[elem]!=0)):
				normalized[row,elem]=kmat[row][elem]/numpy.sqrt(dia[row]*dia[elem])
			else:
				normalized[row,elem]=0.
				#Should be the same as /numpy.sqrt(len(st1)*len(st2)) for small bins
	return normalized.tolist()
	
def kappa(spikes, binwidth): # Returns normalized kappa for entire set of spike trains
	sorted = sort_spikes(spikes)
	kmat=kappamat(sorted, binwidth)
	nkmat = normalize_kappa(kmat)
#	k=numpy.mean(nkmat) # Count empty spike trains too
	npnkmat = numpy.array(nkmat)
	k=numpy.mean(npnkmat[numpy.nonzero(npnkmat)])
	return k