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

i = 1.0

MsynI = np.array([100, 60, 40, 25, 10])    # mean connectivity in %
JI = np.array([5.0, 2.0, 1.0, 0.5, 0.1]) # synaptic weight (unit nS)

MsynIepi = np.array([1., 0.8, 0.6, 0.4, 0.2])*25.
MsynEepi = np.array([1+i*3./4. for i in np.arange(0,5,1)])*10.  # multiplication factor of excitatory connections (baseline 10%, up to 40%)
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
je = 2.0 # not used
me = 10. # not used
fh = open(direc+".sh", 'w')

#MsynI and phi (JI fixed)
subdir = "MI_phi"
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
ji = 5.0
for MM in MsynI:
 for PP in phi:
  print >>fh, makecommand(direc,subdir,PP,i,ji,je,MM,me,NE,NI)
  lines += 1

#JI and phi (MsynI fixed)
subdir = "JI_phi"
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
mi = 25.
for JJ in JI:
 for PP in phi:
  print >>fh, makecommand(direc,subdir,PP,i,JJ,je,mi,me,NE,NI)
  lines += 1

fh.close()
print "Generated shell script", direc+".sh with", lines, "lines"

##########################
# Mixed networks
NE=400
NI=100
direc = "both"
lines = 0
fh = open(direc+".sh", 'w')

#MsynI and phi (MsynE, JI, JE fixed)
subdir = "MI_phi"
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
ji = 5.0
je = 2.0
me = 10.
for MM in MsynIepi:
 for PP in phi:
  print >>fh, makecommand(direc,subdir,PP,i,ji,je,MM,me,NE,NI)
  lines += 1

#MsynE and phi (MsynI, JI, JE fixed)
subdir = "ME_phi"
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
ji = 5.0
je = 2.0
mi = 25.
for MM in MsynEepi:
 for PP in phi:
  print >>fh, makecommand(direc,subdir,PP,i,ji,je,mi,MM,NE,NI)
  lines += 1

#JI and phi (MsynI, MsynE, JE fixed)
subdir = "JI_phi"
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
je = 2.0
mi = 25.
me = 10.
for JJ in JIepi:
 for PP in phi:
  print >>fh, makecommand(direc,subdir,PP,i,JJ,je,mi,me,NE,NI)
  lines += 1

#JE and phi (MsynI, MsynE, JI fixed)
subdir = "JE_phi"
print >>fh, "mkdir output/"+direc+"/"+subdir
print >>fh, "mkdir figures/"+direc+"/"+subdir
ji = 5.0
mi = 25.
me = 10.
for JJ in JEepi:
 for PP in phi:
  print >>fh, makecommand(direc,subdir,PP,i,ji,JJ,mi,me,NE,NI)
  lines += 1

fh.close()
print "Generated shell script", direc+".sh with", lines, "lines"
