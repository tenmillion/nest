# generate_batch_script_var_tref.py
# Script for generating bash files
# For the new simulation conditions
# For variable refractory periods
# Python 2.7
# python my_brunel2000_rand_hh.py phi Iext Iext_s JI JE MsynI MsynE NE NI direc subdir
#direcs = [excit_only, inh_only, both]
#subdir = [connectivity, weights] # for inh_only
#subdir = [connectivity, inh_weights, exc_weights] # for both
import numpy as np

phi = np.array([5*3**(k/10.) for k in np.arange(-36,1,3)]) # temperature -37 degC
print phi

Iext = 0.1
print Iext

MI = np.array([100, 60, 50, 40, 30, 20, 10, 0])    # mean connectivity in %
JI = np.array([2.1, 1.8, 1.5, 1.2, 0.9, 0.6, 0.3, 0.1]) # synaptic weight (unit nS)

print MI
print JI

MIepi = np.array([1., 0.8, 0.6, 0.4, 0.2])*25.
MEepi = np.array([1+k*3./4. for k in np.arange(0,5,1)])*10.  # multiplication factor of excitatory connections (baseline 10%, up to 40%)
JIepi = np.array([1.0, 0.8, 0.6, 0.4, 0.2])*2. # synaptic weight (unit nS)
JEepi = np.array([1.0, 1.1, 1.5, 2.0, 4.0])*5. # synaptic weight (unit nS)
#Iext = np.array([10, 1., 0.5])  # external input current

def makecommand(direc,subdir,phi_,Iext_,JI_,JE_,MsynI_,MsynE_,NI_,NE_):
 cmd='python my_brunel2000_rand_hh.py {:.3f} {:.1f} 0 {:.1f} {:.1f} {:.1f} {:.1f} {} {} '.\
     format(phi_,Iext_,JI_,JE_,MsynI_,MsynE_,NI_,NE_) + direc + ' ' + subdir
 return cmd


##########################
# Inhibitory only networks
NI=100
NE=0
direc = "inh_only"
je = 2.0 # not used
me = 10. # not used

subdir = "MI_JI_i01"
fh = open(direc+"_"+subdir+"_JI"+str(0)+".sh", 'w')
lines = 0
print >>fh, "mkdir output/"+direc
print >>fh, "mkdir figures/"+direc
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
for Jn in range(len(JI)):
 lines = 0
 for MM in MI:
  for PP in phi:
   print >>fh, makecommand(direc,subdir,PP,Iext,JI[Jn],je,MM,me,NI,NE)
   lines += 1
 fh.close()
 print "Generated shell script", direc+"_"+subdir+"_JI"+str(Jn)+".sh with", lines, "lines"
 try:
  fh = open(direc+"_"+subdir+"_JI"+str(Jn+1)+".sh", 'w')
 except:
  print "Reached end of JI:", JI[Jn], "at", Jn

