# -*- coding: utf-8 -*-
#
# my_brunel2000_rand_hh.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

#import NeuroTools.signals # Needs to load before nest, or python core dumps. Not sure why.
#import NeuroTools.io # Needs latest version of NeuroTools downloaded thru svn
import nest
import nest.raster_plot
import nest.voltage_trace
import pylab
import numpy
from numpy import ndarray
from datetime import datetime
import sys
import os

# Network parameters. Some of these are given in Brunel (2000) J.Comp.Neuro.
phi	   = float(sys.argv[1])    # Default 1.
g      = float(sys.argv[2])    # Ratio of IPSP to EPSP amplitude: J_I/J_E
I_ext  = float(sys.argv[3])*100.  # Applied current in pA (default 1 uA/cm^2 = 100 pA)
Is_ext = float(sys.argv[4])*100.  # SD of I_ext in pA (scaled by 100 from original uA/cm^2)
J_E    = float(sys.argv[5]) # Should be 0.1 mS/cm^2. Weight has units in pS in NEST? So... 
J_range = 0.	# Range of synaptic weight (0 to 1)

N_E = 0
N_I = 100

M_syn_EE = 60. # Average number of ex-ex syns
M_syn_EI = 30.
M_syn_IE = 30.
M_syn_II = 30.

plotdistribs = True
plotresults = True

V_init	= -60.	# Initial membrane potential
V_range = 20.	# Range of initial membrane potential
delay   = 1.5   # synaptic delay in ms
d_range = 0.	# Range of synaptic delay (0 to 1)

N_neurons = N_E+N_I

if N_E > 0:
 p_conn_EE = M_syn_EE/float(N_E) # Probability of a synapse existing between ex-ex
 p_conn_EI = M_syn_EI/float(N_E)
else:
 p_conn_EE = 0.
 p_conn_EI = 0.
if N_I > 0:
 p_conn_IE = M_syn_IE/float(N_I)
 p_conn_II = M_syn_II/float(N_I)
else:
 p_conn_IE = 0.
 p_conn_II = 0.

J_I  = -g*J_E

# Set parameters of the NEST simulation kernel
nest.SetKernelStatus({"print_time": True,
                      "local_num_threads": 1})
tdatetime = datetime.now()
tstr = tdatetime.strftime('_%m%d-%H%M-%S_')
fnprefix = str('phi%.1f'%phi)+str('g%.1f'%g)+str('in%.1f'%I_ext)+\
		   str('J%.0f'%J_E)+str('+-%.1f'%J_range)+\
		   str('E%d'%N_E)+str('I%d'%N_I)+str('MsynII%d'%M_syn_II)+tstr
nest.SetKernelStatus({"data_path": "output", "data_prefix": fnprefix})

# Create and seed RNGs
ms = 1000 # master seed
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
pyrngs = [numpy.random.RandomState(s) for s in range(ms, ms+n_vp)]
nest.SetKernelStatus({'grng_seed': ms+n_vp,
                      'rng_seeds': range(ms+n_vp+1, ms+1+2*n_vp)})

# IAF mode for debugging
#nest.SetDefaults("iaf_psc_delta", {"C_m": 1.0, "tau_m":20*phi})
#nodes   = nest.Create("iaf_psc_delta",N_neurons)

nest.SetDefaults("hh_cond_exp_traub", # Using Wang-Buzsaki
                 {"phi_t": phi, # var
				  "I_e": I_ext, # var
				  "g_Na": 3500.,
				  "g_K": 900.,
				  "g_L": 10.,
				  "E_Na": 55.,
				  "E_K": -90.,
				  "E_L": -65.,
				  "E_ex": 0.,
				  "E_in": -75., # var
				  "tau_syn_ex": 10., # var
				  "tau_syn_in": 10., # var
				  "V_T": 0.,
                  "V_1": 35.0, 
                  "V_2": 60.0,
                  "V_3": 58.0,
                  "V_4": 28.0,
                  "V_5": 34.0,
                  "V_6": 44.0})
# These values are already set as defaults in the models.
# Only the variable ones need to be changed.

nodes = nest.Create("hh_cond_exp_traub",N_neurons)
nodes_E= nodes[:N_E]
nodes_I= nodes[N_E:]

# randomize membrane potential
if N_E > 0:
 node_E_info   = nest.GetStatus(nodes_E, ['global_id','vp','local'])
 local_nodes_E = [(gid,vp) for gid,vp,islocal in node_E_info if islocal]
 for gid,vp in local_nodes_E: 
   nest.SetStatus([gid], {'V_m': pyrngs[vp].uniform(V_init-V_range/2.,V_init+V_range/2.)})

if N_I > 0:
 node_I_info   = nest.GetStatus(nodes_I, ['global_id','vp','local'])
 local_nodes_I = [(gid,vp) for gid,vp,islocal in node_I_info if islocal]
 for gid,vp in local_nodes_I: 
   nest.SetStatus([gid], {'V_m': pyrngs[vp].uniform(V_init-V_range/2.,V_init+V_range/2.)})

# Generate connectivity matrix
def flip(p):
    return 1 if numpy.random.random() < p else 0

if N_E > 0:
 conn_EE = []
 for r in range(0,N_E):
  row = []
  for e in range(0,N_E):
   row.append(flip(p_conn_EE))
  conn_EE.append(row)

if N_E > 0 and N_I > 0:
 conn_EI = []
 for r in range(0,N_I):
  row = []
  for e in range(0,N_E):
   row.append(flip(p_conn_EI))
  conn_EI.append(row)

 conn_IE = [] # Inhib to Excit
 for r in range(0,N_E): # Each row contains synaptic connections to a neuron
  row = []
  for e in range(0,N_I):
   row.append(flip(p_conn_IE))
  conn_IE.append(row)

if N_I > 0:
 conn_II = []
 for r in range(0,N_I):
  row = []
  for e in range(0,N_I):
   row.append(flip(p_conn_II))
  conn_II.append(row)

# Count number of synapses per cell
if N_E > 0:
 n_conn_EE = numpy.sum(conn_EE,1) # excit conns received by each excit neuron
 if N_I > 0:
  n_conn_EI = numpy.sum(conn_EI,1) # excit conns received by each inhib neuron
if N_I > 0:
 if N_E > 0:
  n_conn_IE = numpy.sum(conn_IE,1)
 n_conn_II = numpy.sum(conn_II,1)

if N_E > 0 and N_I > 0: 
 print "Numbers of targets and mean/var of numbers of sources of EE, EI, IE, II connections:"
 print len(n_conn_EE), len(n_conn_EI), len(n_conn_IE), len(n_conn_II)
 print numpy.mean(n_conn_EE), numpy.mean(n_conn_EI), numpy.mean(n_conn_IE), numpy.mean(n_conn_II)
 print numpy.var(n_conn_EE), numpy.var(n_conn_EI), numpy.var(n_conn_IE), numpy.var(n_conn_II)

elif N_E > 0:
 print "Numbers of targets and mean/var of numbers of sources of EE connections:"
 print len(n_conn_EE)
 print numpy.mean(n_conn_EE)
 print numpy.var(n_conn_EE)

elif N_I > 0:
 print "Numbers of targets and mean/var of numbers of sources of II connections:"
 print len(n_conn_II)
 print numpy.mean(n_conn_II)
 print numpy.var(n_conn_II)

# Make synaptic connections
i = 0
if N_E > 0:
 nest.CopyModel("static_synapse", "excitatory") # From Excitatory
 for tgt_gid, tgt_vp in local_nodes_E: # To Excitatory
   eweights = pyrngs[tgt_vp].uniform((1-J_range/2.)*J_E, (1+J_range/2.)*J_E, n_conn_EE[i])
   edelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_EE[i])
   nest.RandomConvergentConnect(nodes_E, [tgt_gid], n_conn_EE[i],
                                weight = list(eweights), delay = list(edelay),
                                model="excitatory")
   i = i+1

if N_E >0 and N_I >0:
 i = 0
 for tgt_gid, tgt_vp in local_nodes_I: # To Inhibitory
   eweights = pyrngs[tgt_vp].uniform((1-J_range/2.)*J_E, (1+J_range/2.)*J_E, n_conn_EI[i])
   edelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_EI[i])
   nest.RandomConvergentConnect(nodes_E, [tgt_gid], n_conn_EI[i],
                                weight = list(eweights), delay = list(edelay),
                                model="excitatory")
   i = i+1

 i = 0
 nest.CopyModel("static_synapse", "inhibitory") # From Inhibitory
 for tgt_gid, tgt_vp in local_nodes_E: # To Excitatory
   iweights = pyrngs[tgt_vp].uniform((1-J_range/2.)*J_I, (1+J_range/2.)*J_I, n_conn_IE[i])
   idelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_IE[i])
   nest.RandomConvergentConnect(nodes_I, [tgt_gid], n_conn_IE[i],
                                weight = list(iweights), delay = list(idelay),
                                model="inhibitory")
   i = i+1

else:
 nest.CopyModel("static_synapse", "inhibitory") # From Inhibitory
if N_I > 0:
 i = 0
 for tgt_gid, tgt_vp in local_nodes_I: # To Inhibitory
   iweights = pyrngs[tgt_vp].uniform((1-J_range/2.)*J_I, (1+J_range/2.)*J_I, n_conn_II[i])
   idelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_II[i])
   nest.RandomConvergentConnect(nodes_I, [tgt_gid], n_conn_II[i],
                                weight = list(iweights), delay = list(idelay),
                                model="inhibitory")
   i = i+1

# Modify inputs
if Is_ext > 0:
 if N_E > 0:
  eI_elist = numpy.random.normal(I_ext, Is_ext, N_E)
  for k in range(N_E):
   nest.SetStatus([nodes_E[k]], {"I_e": eI_elist[k]})

 if N_I > 0:
  iI_elist = numpy.random.normal(I_ext, Is_ext, N_I)
  for k in range(N_I):
   nest.SetStatus([nodes_I[k]], {"I_e": iI_elist[k]})

# Make input connections
#noise=nest.Create("poisson_generator",1,{"rate": p_rate*N_neurons})

#nest.CopyModel("static_synapse_hom_wd",
#               "excitatory-input",
#               {"weight":J_E, 
#                "delay":delay})
#nest.DivergentConnect(noise,nodes,model="excitatory-input")

# Make and connect recording devices
N_rec = 50    # Number of neurons to record from
if N_E > 0:
 spikes_E=nest.Create("spike_detector",1, 
                    [{"label": "brunel-py-ex", "withtime": True,
 					"withgid": True, "to_file": True, "start":1000.}])                   
 nest.ConvergentConnect(nodes_E[:N_rec],spikes_E)
 voltmeter_E = nest.Create("voltmeter")
 nest.SetStatus(voltmeter_E,[{"to_file": True, "withtime": True, "withgid": True, "start":1000.}])
 nest.Connect(voltmeter_E, nodes_E[1:2])

if N_I > 0:
 spikes_I=nest.Create("spike_detector",1, 
                    [{"label": "brunel-py-in", "withtime": True,
	  	 			"withgid": True, "to_file": True, "start":1000.}])                   
 nest.ConvergentConnect(nodes_I[:N_rec],spikes_I)

 voltmeter_I = nest.Create("voltmeter")
 nest.SetStatus(voltmeter_I,[{"to_file": True, "withtime": True, "withgid": True, "start":1000.}])
 nest.Connect(voltmeter_I, nodes_I[1:2])

# Visualization of initial membrane potential and initial weight
# distribution only if we run on single MPI process
if nest.NumProcesses() == 1:
 if plotdistribs:
  pylab.figure()
  V_I = nest.GetStatus(nodes_I[:N_rec], 'V_m')
  pylab.hist(V_I, bins=10)
  pylab.xlabel('Membrane potential V_m [mV]')
  pylab.savefig('./figures/'+fnprefix+'Vm_dist.eps')

  pylab.figure()
  w = nest.GetStatus(nest.GetConnections(nodes_I[:N_rec],
                                          synapse_model='inhibitory'),
                     'weight')
  pylab.hist(w, bins=100)
  pylab.xlabel('Synaptic weight [pA]')
  pylab.title('Distribution of synaptic weights (%d synapses)' % len(w))
  pylab.savefig('./figures/'+fnprefix+'w_dist.eps')

  pylab.figure()
  d = nest.GetStatus(nest.GetConnections(nodes_I[:N_rec],
                                          synapse_model='inhibitory'),
                     'delay')
  pylab.hist(d, bins=100)
  pylab.xlabel('Synaptic delay [ms]')
  pylab.title('Distribution of synaptic delay (%d synapses)' % len(d))
  pylab.savefig('./figures/'+fnprefix+'d_dist.eps')

  pylab.figure()
  if N_E > 0:
   pylab.hist(eI_elist,bins=N_E/3)
   pylab.xlabel("Input current in pA")
  if N_I > 0:
   pylab.hist(iI_elist,bins=N_I/3)
   pylab.xlabel("Input current in pA")
  pylab.show()

else:
  print "Multiple MPI processes, skipping graphical output"

simtime   = 1300.  # how long shall we simulate [ms]
nest.Simulate(simtime)

# Before we compute the rates, we need to know how many of the recorded
# neurons are on the local MPI process
if N_E>0:
 events_E = nest.GetStatus(spikes_E,"n_events")
 N_rec_local_E = sum(nest.GetStatus(nodes_E[:N_rec], 'local'))
 rate_ex= events_E[0]/simtime*1000.0/N_rec_local_E
 print "Excitatory rate   : %.2f Hz" % rate_ex
 filename_E = nest.GetStatus(spikes_E,"filenames")
 os.rename(filename_E[0][0], filename_E[0][0][:-3]+"txt")
 print filename_E[0][0][:-3]
 spiketrain = numpy.loadtxt(filename_E[0][0][:-3]+"txt")  
 print spiketrain.shape, "size of spiketrain_E"

if N_I>0:
 events_I = nest.GetStatus(spikes_I,"n_events")
 print events_I, "events_I"
 N_rec_local_I = sum(nest.GetStatus(nodes_I[:N_rec], 'local'))
 rate_in= events_I[0]/simtime*1000.0/N_rec_local_I
 print "Inhibitory rate   : %.2f Hz" % rate_in
 filename_I = nest.GetStatus(spikes_I,"filenames")
 print filename_I[0][0], "before"
 os.rename(filename_I[0][0], filename_I[0][0][:-3]+"txt")
 print filename_I[0][0][:-3]
 #for filename in os.listdir("./output"):
 # if filename == filename_I[0][0][:-3]+"txt":
 #  print filename
 spiketrain = numpy.loadtxt(filename_I[0][0][:-3]+"txt")  

# Plot Results
if plotresults:
 pylab.figure()
 if N_E>0:
  nest.voltage_trace.from_device(voltmeter_E)
 if N_I>0:
  nest.voltage_trace.from_device(voltmeter_I)
  pylab.savefig('./figures/'+fnprefix+'Vm_trace.eps')
  nest.voltage_trace.show()
 
 if nest.NumProcesses() == 1:
   if N_E>0:
    nest.raster_plot.from_device(spikes_E, hist=True, title='Excitatory')
    pylab.savefig('./figures/'+fnprefix+'raster.eps')
   if N_I>0:  
    nest.raster_plot.from_device(spikes_I, hist=True, title='Inhibitory')
    pylab.savefig('./figures/'+fnprefix+'raster.eps')
 else:
   print "Multiple MPI processes, skipping graphical output"

 pylab.show()
else:
 print "Set plotresults to True to show plots here."
