#! /usr/bin/env python

import nest
import nest.voltage_trace
import sys

nest.ResetKernel()
nest.SetKernelStatus({"overwrite_files": True})

V_init	= -65.
phi = float(sys.argv[1])

nest.SetDefaults("hh_psc_alpha", 
                 {"E_L": V_init,
				  "C_m": 100.0,
				  "phi_t": phi, # More: less tau
                  "V_1": 35.0, # Using Wang-Buzsaki
                  "V_2": 60.0,
                  "V_3": 58.0,
                  "V_4": 28.0,
                  "V_5": 34.0,
                  "V_6": 44.0})

neuron = nest.Create("hh_psc_alpha")
noise = nest.Create("poisson_generator", 2)

voltmeter = nest.Create("voltmeter")

nest.SetStatus(noise,[{"rate": 130000.0}, {"rate": 10000.0}])

nest.SetStatus(voltmeter,[{"to_file": True, "withtime": True}])

nest.ConvergentConnect(noise, neuron, [1.2, -1.], [1.0, 1.0])
nest.Connect(voltmeter, neuron)

nest.Simulate(500.0)

nest.voltage_trace.from_device(voltmeter)
nest.voltage_trace.show()
