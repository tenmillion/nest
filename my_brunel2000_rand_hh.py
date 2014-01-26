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

import nest
import nest.raster_plot
import nest.voltage_trace
import pylab
import numpy
from datetime import datetime
import sys

# Network parameters. Some of these are given in Brunel (2000) J.Comp.Neuro.
phi		= float(sys.argv[1])    # Default 1.
g       = float(sys.argv[2])    # Ratio of IPSP to EPSP amplitude: J_I/J_E
p_rate  = float(sys.argv[3])    # rate of external population

delay   = 1.5   # synaptic delay in ms
V_init	= -65.	# Initial membrane potential, set to same value as E_L
E_L	    = -65.
V_range = 10.	# Range of initial membrane potential
d_range = 0.	# Range of synaptic delay (0 to 1)
N_E = 400
N_I = 100
N_neurons = N_E+N_I

M_syn_EE = 1. # Average number of ex-ex syns
M_syn_EI = 1.
M_syn_IE = 1.
M_syn_II = 1.

p_conn_EE = M_syn_EE/float(N_E) # Probability of a synapse existing between ex-ex
p_conn_EI = M_syn_EI/float(N_E)
p_conn_IE = M_syn_IE/float(N_I)
p_conn_II = M_syn_II/float(N_I)

J_E  = 1. # Should be 0.1 mS/cm^2. Not sure what this corresponds to in NEST.
J_I  = -g*J_E

# Set parameters of the NEST simulation kernel
nest.SetKernelStatus({"print_time": True,
                      "local_num_threads": 1})
tdatetime = datetime.now()
tstr = tdatetime.strftime('_%m%d-%H%M-%S_')
fnprefix = "phi"+str('%.1f'%phi)+"g"+str('%.1f'%g)+"in"+str('%.0f'%p_rate)+tstr
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
node_E_info   = nest.GetStatus(nodes_E, ['global_id','vp','local'])
local_nodes_E = [(gid,vp) for gid,vp,islocal in node_E_info if islocal]
for gid,vp in local_nodes_E: 
  nest.SetStatus([gid], {'V_m': pyrngs[vp].uniform(V_init-V_range/2.,V_init+V_range/2.)})

node_I_info   = nest.GetStatus(nodes_I, ['global_id','vp','local'])
local_nodes_I = [(gid,vp) for gid,vp,islocal in node_I_info if islocal]
for gid,vp in local_nodes_I: 
  nest.SetStatus([gid], {'V_m': pyrngs[vp].uniform(V_init-V_range/2.,V_init+V_range/2.)})

# Generate connectivity matrix
def flip(p):
    return 1 if numpy.random.random() < p else 0

conn_EE = []
for r in range(0,N_E):
 row = []
 for e in range(0,N_E):
  row.append(flip(p_conn_EE))
 conn_EE.append(row)

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

conn_II = []
for r in range(0,N_I):
 row = []
 for e in range(0,N_I):
  row.append(flip(p_conn_II))
 conn_II.append(row)

# Count number of synapses per cell
n_conn_EE = numpy.sum(conn_EE,1) # excit conns received by each excit neuron
n_conn_EI = numpy.sum(conn_EI,1) # excit conns received by each inhib neuron
n_conn_IE = numpy.sum(conn_IE,1)
n_conn_II = numpy.sum(conn_II,1)

print "Numbers of targets and mean/var of numbers of sources of EE, EI, IE, II connections:"
print len(n_conn_EE), len(n_conn_EI), len(n_conn_IE), len(n_conn_II)
print numpy.mean(n_conn_EE), numpy.mean(n_conn_EI), numpy.mean(n_conn_IE), numpy.mean(n_conn_II)
print numpy.var(n_conn_EE), numpy.var(n_conn_EI), numpy.var(n_conn_IE), numpy.var(n_conn_II)

# Make connections
i = 0
nest.CopyModel("static_synapse", "excitatory") # From Excitatory
for tgt_gid, tgt_vp in local_nodes_E: # To Excitatory
  eweights = pyrngs[tgt_vp].uniform(0.5*J_E, 1.5*J_E, n_conn_EE[i])
  edelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_EE[i])
  nest.RandomConvergentConnect(nodes_E, [tgt_gid], n_conn_EE[i],
                               weight = list(eweights), delay = list(edelay),
                               model="excitatory")
  i = i+1

i = 0
for tgt_gid, tgt_vp in local_nodes_I: # To Inhibitory
  eweights = pyrngs[tgt_vp].uniform(0.5*J_E, 1.5*J_E, n_conn_EI[i])
  edelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_EI[i])
  nest.RandomConvergentConnect(nodes_E, [tgt_gid], n_conn_EI[i],
                               weight = list(eweights), delay = list(edelay),
                               model="excitatory")
  i = i+1

i = 0
nest.CopyModel("static_synapse", "inhibitory") # From Inhibitory
for tgt_gid, tgt_vp in local_nodes_E: # To Excitatory
  iweights = pyrngs[tgt_vp].uniform(0.5*J_I, 1.5*J_I, n_conn_IE[i])
  idelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_IE[i])
  nest.RandomConvergentConnect(nodes_I, [tgt_gid], n_conn_IE[i],
                               weight = list(iweights), delay = list(idelay),
                               model="inhibitory")
  i = i+1

i = 0
for tgt_gid, tgt_vp in local_nodes_I: # To Inhibitory
  iweights = pyrngs[tgt_vp].uniform(0.5*J_I, 1.5*J_I, n_conn_II[i])
  idelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, n_conn_II[i])
  nest.RandomConvergentConnect(nodes_I, [tgt_gid], n_conn_II[i],
                               weight = list(iweights), delay = list(idelay),
                               model="inhibitory")
  i = i+1

noise=nest.Create("poisson_generator",1,{"rate": p_rate*N_neurons})

nest.CopyModel("static_synapse_hom_wd",
               "excitatory-input",
               {"weight":J_E, 
                "delay":delay})
nest.DivergentConnect(noise,nodes,model="excitatory-input")

spikes=nest.Create("spike_detector",2, 
                   [{"label": "brunel-py-ex", "withtime": True,
					"withgid": True, "to_file": True, "start":1000.},
                    {"label": "brunel-py-in", "withtime": True,
					"withgid": True, "to_file": True, "start":1000.}])
                   
spikes_E=spikes[:1]
spikes_I=spikes[1:]

N_rec = 50    # Number of neurons to record from
nest.ConvergentConnect(nodes_E[:N_rec],spikes_E)
nest.ConvergentConnect(nodes_I[:N_rec],spikes_I)

voltmeter_E = nest.Create("voltmeter")
voltmeter_I = nest.Create("voltmeter")
nest.SetStatus(voltmeter_E,[{"to_file": True, "withtime": True, "withgid": True, "start":1000.}])
nest.SetStatus(voltmeter_I,[{"to_file": True, "withtime": True, "withgid": True, "start":1000.}])
nest.Connect(voltmeter_E, nodes_E[1:2])
nest.Connect(voltmeter_I, nodes_I[1:2])

# Visualization of initial membrane potential and initial weight
# distribution only if we run on single MPI process
if nest.NumProcesses() == 1:
  pylab.figure()
  V_E = nest.GetStatus(nodes_E[:N_rec], 'V_m')
  pylab.hist(V_E, bins=10)
  pylab.xlabel('Membrane potential V_m [mV]')
  pylab.savefig('./figures/'+fnprefix+'rand_Vm.eps')

  pylab.figure()
  w = nest.GetStatus(nest.GetConnections(nodes_E[:N_rec],
                                          synapse_model='excitatory'),
                     'weight')
  pylab.hist(w, bins=100)
  pylab.xlabel('Synaptic weight [pA]')
  pylab.title('Distribution of synaptic weights (%d synapses)' % len(w))
  pylab.savefig('./figures/'+fnprefix+'rand_w.eps')

  pylab.figure()
  d = nest.GetStatus(nest.GetConnections(nodes_E[:N_rec],
                                          synapse_model='excitatory'),
                     'delay')
  pylab.hist(d, bins=100)
  pylab.xlabel('Synaptic delay [ms]')
  pylab.title('Distribution of synaptic delay (%d synapses)' % len(d))
  pylab.savefig('./figures/'+fnprefix+'rand_d.eps')
else:
  print "Multiple MPI processes, skipping graphical output"

simtime   = 1300.  # how long shall we simulate [ms]
nest.Simulate(simtime)

pylab.figure()
nest.voltage_trace.from_device(voltmeter_E)
nest.voltage_trace.from_device(voltmeter_I)
nest.voltage_trace.show()

events = nest.GetStatus(spikes,"n_events")

# Before we compute the rates, we need to know how many of the recorded
# neurons are on the local MPI process
N_rec_local_E = sum(nest.GetStatus(nodes_E[:N_rec], 'local'))
rate_ex= events[0]/simtime*1000.0/N_rec_local_E
print "Excitatory rate   : %.2f Hz" % rate_ex

N_rec_local_I = sum(nest.GetStatus(nodes_I[:N_rec], 'local'))
rate_in= events[1]/simtime*1000.0/N_rec_local_I
print "Inhibitory rate   : %.2f Hz" % rate_in

if nest.NumProcesses() == 1:
  nest.raster_plot.from_device(spikes_E, hist=True, title='Excitatory')
  pylab.savefig('./figures/'+fnprefix+'rand_raster.eps')
  nest.raster_plot.from_device(spikes_I, hist=True, title='Inhibitory')
  pylab.savefig('./figures/'+fnprefix+'rand_raster.eps')
else:
  print "Multiple MPI processes, skipping graphical output"

pylab.show()
