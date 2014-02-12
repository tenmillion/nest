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
	if len(spikes_tr.shape) > 1:
		flag = 0
		mincell=numpy.int(min(spikes_tr[0]))
		maxcell=numpy.int(max(spikes_tr[0])+1)
		sorted_spiketrains = []
		for cell in range(mincell,maxcell):
		  locs = (numpy.where(spikes_tr[0] == cell+1))[0] #where returns tuple
		  sorted_spiketrains.append(spikes_tr[1,locs])
		return sorted_spiketrains, flag
	else:
		flag = 1
		return spikes_tr, flag

def pairwise_kappa(s1, s2, binwidth, tstart, tstop): # expects numpy arrays
	w1 = s1[numpy.where((s1 >= tstart) & (s1 <= tstop))] # Windowed s1
	w2 = s2[numpy.where((s2 >= tstart) & (s2 <= tstop))]
	overlaps = 0
	for k in xrange(0, w1.shape[0]):
		overlaps += numpy.count_nonzero((w1[k] - binwidth/2. <= w2) & (w2 <= w1[k] + binwidth/2.))
#		print overlaps, "overlaps"
#		print "spike trains:", id1, "("+str(st1.shape[0])+" spikes) and", id2, "("+str(st2.shape[0])+" spikes)"
	return overlaps
	
# Somehow when I import functions the array[i][j] syntax doesn't work even for numpy arrays
def kappamat(sortedsps, binwidth, tstart, tstop):
	xymat = []
	for s1 in range(len(sortedsps)):
		row = []
		st1=sortedsps[s1]
		for s2 in range(len(sortedsps)):
			st2=sortedsps[s2]
			if ((len(st1)!=0) and (len(st2)!=0)):
				row.append(pairwise_kappa(st1, st2, binwidth, tstart, tstop)) # not normalized
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
	
def kappa(spikes, binwidth, tstart, tstop): # Returns normalized kappa for entire set of spike trains
	sortedsps, flag = sort_spikes(spikes)
	if flag == 0:
		kmat=kappamat(sortedsps, binwidth, tstart, tstop)
		nkmat = normalize_kappa(kmat)
	#	k=numpy.mean(nkmat) # Count empty spike trains too
		npnkmat = numpy.array(nkmat)
		k=numpy.mean(npnkmat[numpy.nonzero(npnkmat)])
	else:
		k=-1
	return k
