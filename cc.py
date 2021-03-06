import nest
import nest.raster_plot
import pylab

import numpy
from numpy import exp

from NeuroTools import signals

sd_filename = "output/tau60.0g7.5in9521.0_0124-0111-14_brunel-iaf-E-252-0"
spikes = signals.load_spikelist(sd_filename+".gdf", dims = 1, id_list=range(0,50))

ebin = 0.01 / spikes.mean_rate() * 1e3

cc = spikes.pairwise_cc(5000, time_bin=ebin, averaged=False, display=True)

pylab.savefig(sd_filename+'-cc.eps')
pylab.show()

#print "Pearson's CC	(ex) : %.2f" % cc % numpy.ndarray
