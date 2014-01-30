# generate_batch_script_new.py
# Script for generating bash files
# For the new simulation conditions
# Python 2.7
#python my_brunel2000_rand_hh.py phi Iext Iext_s JI JE MsynI MsynE NE NI direc subdir
#direcs = [excit_only, inh_only, both]
#subdir = [connectivity, weights] # for inh_only
#subdir = [connectivity, inh_weights, exc_weights] # for both
import numpy as np

phi = np.array([5*3**((i-37)/10.) for i in np.arange(10,41,3)]) # temperature
print phi

MsynI = np.array([100, 60, 40, 25, 10])    # mean connectivity in %
JI = np.array([5.0, 2.0, 1.0, 0.5, 0.1]) # synaptic weight (unit nS)

MsynIepi = np.array([1., 0.8, 0.6, 0.4, 0.2, 0.1])*25.
MsynEepi = np.array([1+i*3./5. for i in np.arange(0,6,1)])*10.  # multiplication factor of excitatory connections (baseline 10%, up to 40%)
JIepi = np.array([1.0, 0.8, 0.6, 0.4, 0.2])*2. # synaptic weight (unit nS)
JEepi = np.array([1.0, 1.1, 1.5, 2.0, 4.0])*5. # synaptic weight (unit nS)
#Iext = np.array([10, 1., 0.5])  # external input current

def makecommand(direc,subdir,phi_,Iext_,JI_,JE_,MsynI_,MsynE_,NE_,NI_):
 cmd='python my_brunel2000_rand_hh.py {:.2f} {:.1f} 0 {:.1f} {:.1f} {:.1f} {:.1f} {} {} '.\
     format(phi_,Iext_,JI_,JE_,MsynI_,MsynE_,NE_,NI_) + direc + ' ' + subdir
 return cmd


##########################
# Inhibitory only networks
NE=0
NI=100
direc = "inh_only"
lines = 0
fh = open(direc+".sh", 'w')

#MsynI and phi
subdir = "connectivity"
print >>fh, "mkdir "+subdir
j = 2.0
i = 1.0
for mm in Msyn:
 for pp in phi:
  print >>fh, makecommand(direc,subdir,pp,i,j,mm,m,NE,NI)
  lines += 1

#JI and phi
subdir = "JI_phi"
m =100
i = Iext[2] # Not decided yet
for jj in JI:
 for pp in phi:
  print >>fh, makecommand(direc,subdir,pp,i,jj,m,NE,NI)
  lines += 1

fh.close()
print "Generated shell script", direc+".sh with", lines, "lines"

##########################
# Mixed networks
NE=400
NI=100
direc = "both"
lines = 0

#Iext and phi
subdir = "Iext_phi"
j = 0.1
m = 100
fh = open(direc+".sh", 'w')
for ii in Iext[[0,2,4]]:
 for pp in phi[[0,2,4]]:
  print >>fh, makecommand(direc,subdir,pp,ii,j,m,NE,NI)
  lines += 1

#Msyn and JI
subdir = "Msyn_JI"
p = 5
i = Iext[2] # Not decided yet
for mm in Msyn[[0,2,4]]:
 for jj in JI[[0,2,4]]:
  print >>fh, makecommand(direc,subdir,p,i,jj,mm,NE,NI)
  lines += 1

#Msyn and phi
subdir = "Msyn_phi"
j = 0.1
i = Iext[2] # Not decided yet
for mm in Msyn[[0,2,4]]:
 for pp in phi[[0,2,4]]:
  print >>fh, makecommand(direc,subdir,pp,i,j,mm,NE,NI)
  lines += 1

#JI and phi
subdir = "JI_phi"
m =100
i = Iext[2] # Not decided yet
for jj in JI[[0,2,4]]:
 for pp in phi[[0,2,4]]:
  print >>fh, makecommand(direc,subdir,pp,i,jj,m,NE,NI)
  lines += 1

fh.close()
print "Generated shell script", direc+".sh with", lines, "lines"
