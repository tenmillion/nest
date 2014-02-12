import nest.raster_plot
import nest.voltage_trace
import nest
import os
import pylab as plt
import numpy as np

# Boilerplate
nest.ResetKernel()
nest.SetKernelStatus({"print_time": True, "local_num_threads": 1})
nest.SetKernelStatus({"resolution": .01})
nest.SetKernelStatus({"overwrite_files": False})
nest.SetKernelStatus({"data_path": "./", "data_prefix": "BG"})

nest.SetDefaults("hh_cond_exp_traub", # Using Borg-Graham inhibitory
                 {"g_Na": 8000.,
				  "g_K": 9000.,
				  "g_L": 15.,
				  "E_Na": 45.,
				  "E_K": -90.,
				  "E_L": -61.,
				  "C_m": 0.221,
				  "V_T": -25., # var
				  "t_ref": 10. # Change t_ref according to phi
})

# n= 20 to get 20 neurons
hh_neurons = nest.Create("hh_cond_exp_traub", n=100)
nest.SetStatus(hh_neurons, {'V_m': -65.})

for k in range(100):
	# note that nest.SetStatus's first input has to be a list
	# even if it is length 1
	nest.SetStatus([hh_neurons[k]],{"I_e": k * 0.000005+4.5965})
sd = nest.Create('spike_detector')
nest.SetStatus(sd, {'to_file': True})

# Using ConvergentConnect to project all neurons to one detector
nest.ConvergentConnect(hh_neurons, sd)
vm = nest.Create("voltmeter")
nest.SetStatus(vm,[{"to_file": True, "withtime": True, "withgid": True, "start":0.}])
nest.Connect(vm, hh_neurons[22:23])

vm2 = nest.Create("voltmeter")
nest.SetStatus(vm2,[{"to_file": True, "withtime": True, "withgid": True, "start":0.}])
nest.Connect(vm2, hh_neurons[23:24])

nest.Simulate(1500.)

plt.figure()
nest.voltage_trace.from_device(vm)
plt.show()
nest.voltage_trace.from_device(vm2)
plt.show()
nest.raster_plot.from_device(sd, hist=False, title='raster plot')
plt.show()

# Find out the file name
sd_filename = 'BGspike_detector-' + str(sd[0]) + '-0.gdf'
print sd_filename
os.rename(sd_filename, sd_filename[:-3]+"txt")

# Load the spikes by giving the list of neurons and the data file object
spikes = np.loadtxt(sd_filename[:-3]+"txt")

def sort_spikes(sts,ids):
	spikes_tr=np.transpose(sts)
	sorted_spiketrains = []
	for cellid in ids:
	  locs = (np.where(spikes_tr[0] == cellid+1))[0] #where returns tuple
	  sorted_spiketrains.append(spikes_tr[1,locs])
	return sorted_spiketrains

spikes = sort_spikes(spikes,np.arange(0,100,1))
print spikes[0:30]
#print spikes[0][-1]
rates = []
for st in spikes:
	rates.append(np.shape(np.where(st>500.)[0])[0]/1.)
print rates
print len(rates)
plt.show()
plt.plot(np.arange(0,100,1)*0.000005+4.5965,rates)
plt.show()
