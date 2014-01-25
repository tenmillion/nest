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
delay   = 1.5    # synaptic delay in ms
V_init	= -65.	# Initial membrane potential, set to same value as E_L
E_L	    = -65.
V_range = 10.	 # Range of initial membrane potential
d_range = 0.	# Range of synaptic delay (0 to 1)
N_E = 200
N_I = 50
N_neurons = N_E+N_I

C_E    = N_E/10 # number of excitatory synapses per neuron
C_I    = N_I/10 # number of inhibitory synapses per neuron  

J_E  = 10.
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

nest.SetDefaults("hh_psc_alpha", 
                 {"E_L": -54.402,
				  "C_m": 100.0,
				  "phi_t": phi,
                  "V_1": 40.0, # Using Hodgkin-Huxley
                  "V_2": 65.0,
                  "V_3": 65.0,
                  "V_4": 35.0,
                  "V_5": 55.0,
                  "V_6": 65.0})

#                  "E_L": 54.4, # Using Wang-Buzsaki Todo: adjust other E,g for WB too
#                  "V_1": 35.0, # Using Wang-Buzsaki
#                  "V_2": 60.0,
#                  "V_3": 58.0,
#                  "V_4": 28.0,
#                  "V_5": 34.0,
#                  "V_6": 44.0})

nodes = nest.Create("hh_psc_alpha",N_neurons)
nodes_E= nodes[:N_E]
nodes_I= nodes[N_E:]

# randomize membrane potential
node_info   = nest.GetStatus(nodes, ['global_id','vp','local'])
local_nodes = [(gid,vp) for gid,vp,islocal in node_info if islocal]
for gid,vp in local_nodes: 
  nest.SetStatus([gid], {'V_m': pyrngs[vp].uniform(V_init-V_range/2.,V_init+V_range/2.)})

nest.CopyModel("static_synapse", "excitatory")
for tgt_gid, tgt_vp in local_nodes:
  eweights = pyrngs[tgt_vp].uniform(0.5*J_E, 1.5*J_E, C_E)
  edelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, C_E)
  nest.RandomConvergentConnect(nodes_E, [tgt_gid], C_E,
                               weight = list(eweights), delay = list(edelay),
                               model="excitatory")

nest.CopyModel("static_synapse", "inhibitory")
for tgt_gid, tgt_vp in local_nodes:
  iweights = pyrngs[tgt_vp].uniform(0.5*J_I, 1.5*J_I, C_I)
  idelay = pyrngs[tgt_vp].uniform((1-d_range/2.)*delay, (1+d_range/2.)*delay, C_I)
  nest.RandomConvergentConnect(nodes_I, [tgt_gid], C_I,
                               weight = list(iweights), delay = list(idelay),
                               model="inhibitory")

#nest.CopyModel("static_synapse", "excitatory")
#for tgt_gid, tgt_vp in local_nodes:
#  eweights = numpy.absolute(pyrngs[tgt_vp].normal(J_E, numpy.sqrt(0.5*J_E), C_E))
#  nest.RandomConvergentConnect(nodes_E, [tgt_gid], C_E,
#                               weight = list(eweights), delay = delay,
#                               model="excitatory")
#
#nest.CopyModel("static_synapse", "inhibitory")
#for tgt_gid, tgt_vp in local_nodes:
#  iweights = -numpy.absolute(pyrngs[tgt_vp].normal(-J_I, numpy.sqrt(-0.5*J_I), C_I))
#  nest.RandomConvergentConnect(nodes_I, [tgt_gid], C_I,
#                               weight = list(iweights), delay = delay,
#                               model="inhibitory")

noise=nest.Create("poisson_generator",1,{"rate": p_rate})

nest.CopyModel("static_synapse_hom_wd",
               "excitatory-input",
               {"weight":J_E, 
                "delay":delay})
nest.DivergentConnect(noise,nodes,model="excitatory-input")

spikes=nest.Create("spike_detector",2, 
                   [{"label": "brunel-py-ex", "withtime": True,
					"withgid": True, "to_file": True, "start":2700.},
                    {"label": "brunel-py-in", "withtime": True,
					"withgid": True, "to_file": True, "start":2700.}])
                   
spikes_E=spikes[:1]
spikes_I=spikes[1:]

N_rec = 50    # Number of neurons to record from
nest.ConvergentConnect(nodes_E[:N_rec],spikes_E)
nest.ConvergentConnect(nodes_I[:N_rec],spikes_I)

voltmeter_E = nest.Create("voltmeter")
voltmeter_I = nest.Create("voltmeter")
nest.SetStatus(voltmeter_E,[{"to_file": True, "withtime": True, "withgid": True, "start":2700.}])
nest.SetStatus(voltmeter_I,[{"to_file": True, "withtime": True, "withgid": True, "start":2700.}])
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

simtime   = 3000.  # how long shall we simulate [ms]
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
