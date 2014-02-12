import nest.raster_plot
import nest
import os
import pylab as plt
import numpy as np

# Load the spikes by giving the list of neurons and the data file object
spikes = np.loadtxt("BGspike_detector-201-0.txt")

def sort_spikes(sts,ids):
	spikes_tr=np.transpose(sts)
	sorted_spiketrains = []
	for cellid in ids:
	  locs = (np.where(spikes_tr[0] == cellid+1))[0] #where returns tuple
	  sorted_spiketrains.append(spikes_tr[1,locs])
	return sorted_spiketrains

spikes = sort_spikes(spikes,np.arange(0,200,1))
print spikes[0:30]
#print spikes[0][-1]
rates = []
for st in spikes:
	rates.append(np.shape(np.where(st>0.)[0])[0]/2.)
print rates
print len(rates)
plt.show()
plt.plot(np.arange(0,200,1),rates)
plt.show()
