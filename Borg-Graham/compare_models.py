# WB model TODO: CHANGE
# alphan = 0.01*(55+x)/(1-exp((55+x)/-10))
# betan = 0.125*exp((65+x)/-80)
# alpham = 0.1*(x+40)/(exp((x+40)/10)-1)
# betam = 4*exp((v+65)/-18)
# alphah = 0.07*exp(v+65)/-20)
# betah = 1/(exp((v+35)/-10)+1)

# Borg-Graham
#    const double alpha_m = 4.20 * std::exp(3.01*(y_[0]+38.)*96.489/8.31441/310.);
#    const double beta_m  = 4.20 * std::exp(-1.29*(y_[0]+38.)*96.489/8.31441/310.);
#    const double alpha_h = 0.20 * std::exp(-3.00*(y_[0]+42.)*96.489/8.31441/310.);
#    const double beta_h  = 0.20 * std::exp(3.00*(y_[0]+42.)*96.489/8.31441/310.);
#    const double alpha_n = 0.03 * std::exp(2.10*(y_[0]+35.)*96.489/8.31441/297.);
#    const double beta_n  = 0.03 * std::exp(-0.90*(y_[0]+35.)*96.489/8.31441/297.);


from __future__ import division
import numpy as np
import pylab as plt

## Functions
# K channel
alpha_n = lambda v: 0.03 * np.exp(2.10*(v+35.)*96.489/8.31441/297.)
beta_n  = lambda v: 0.03 * np.exp(-0.90*(v+35.)*96.489/8.31441/297.)
n_inf   = lambda v: alpha_n(v)/(alpha_n(v) + beta_n(v))

# Na channel (activating)
alpha_m = lambda v: 4.20 * np.exp(3.01*(v+38.)*96.489/8.31441/310.)
beta_m  = lambda v: 4.20 * np.exp(-1.29*(v+38.)*96.489/8.31441/310.)
m_inf   = lambda v: alpha_m(v)/(alpha_m(v) + beta_m(v))

# Na channel (inactivating)
alpha_h = lambda v: 0.20 * np.exp(-3.00*(v+42.)*96.489/8.31441/310.)
beta_h  = lambda v: 0.20 * np.exp(3.00*(v+42.)*96.489/8.31441/310.)
h_inf   = lambda v: alpha_h(v)/(alpha_h(v) + beta_h(v))

### channel activity ###
v = np.arange(-100,50) # mV
plt.figure()
plt.plot(v, m_inf(v), v, h_inf(v), v, n_inf(v))
plt.legend(('m','h','n'))
plt.title('Steady state values of ion channel gating variables')
plt.ylabel('Magnitude')
plt.xlabel('Voltage (mV)')

## setup parameters and state variables
T     = 100    # ms
dt    = 0.001 # ms
time  = np.arange(0,T+dt,dt)

## HH Parameters
V_rest  = -65      # mV
Cm      = 2.21      # uF/cm2
gbar_Na = 80.    # mS/cm2
gbar_K  = 90.    # mS/cm2
gbar_l  = 0.15   # mS/cm2
E_Na    = 45    # mV
E_K     = -90    # mV
E_l     = -61 # mV

Vm      = np.zeros(len(time)) # mV
Vm[0]   = V_rest
m      = np.zeros(len(time))
h      = np.zeros(len(time))
n      = np.zeros(len(time))
m[0]       = m_inf(V_rest)      
h[0]       = h_inf(V_rest)
n[0]       = n_inf(V_rest)

I		= 0.3 # uA/cm2

## Variable temperature
phi = np.zeros(len(time))
for i, t in enumerate(time):
  if 0 <= t: phi[i] = 1.

## Simulate Model
for i in range(1,len(time)):
  g_Na = gbar_Na*(m[i-1]**3)*h[i-1]
  g_K  = gbar_K*(n[i-1]**4)
  g_l  = gbar_l

  m[i] = m[i-1] + dt*(alpha_m(Vm[i-1])*(1 - m[i-1]) - beta_m(Vm[i-1])*m[i-1])*phi[i-1]
  h[i] = h[i-1] + dt*(alpha_h(Vm[i-1])*(1 - h[i-1]) - beta_h(Vm[i-1])*h[i-1])*phi[i-1]
  n[i] = n[i-1] + dt*(alpha_n(Vm[i-1])*(1 - n[i-1]) - beta_n(Vm[i-1])*n[i-1])*phi[i-1]

  Vm[i] = Vm[i-1] + (I - g_Na*(Vm[i-1] - E_Na) - g_K*(Vm[i-1] - E_K) - g_l*(Vm[i-1] - E_l)) / Cm * dt 

## Write to file
indices = np.arange(0,int(len(time)/10)+1,1)
matrix = np.vstack((time,Vm,m,h,n))
print "d1",matrix[0]
print "d2",matrix[:,(0,10,20)]
np.savetxt("Vmhn.txt",np.transpose(matrix[:,indices]), fmt='%3.6f', newline='\n')
  
## Plot membrane potential trace
plt.figure()
plt.plot(time, Vm, time, phi, time, m, time, h, time, n)
plt.title('Hodgkin & Huxley Example')
plt.ylabel('Membrane Potential (mV)')
plt.xlabel('Time (msec)')

plt.show()