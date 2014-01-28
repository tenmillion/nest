# generate_batch_script.py
# Script for generating bash files
# Python 2.7
#python my_brunel2000_rand_hh.py phi g Iext Iext_s JI Msyn NE NI direc subdir
#direcs = [excit_only, inh_only, both]
#subdir = [Iext_phi, JI_phi, Msyn_phi, Msyn_JI]
import numpy

Msyn = numpy.array([100, 80, 60, 40, 20])     # mean connectivity in %
JI = numpy.array([1.0, 0.5, 0.1, 0.05, 0.01]) # synaptic weight (unit ?)
phi = numpy.array([5*3**(numpy.float(-i)/2.) for i in numpy.arange(5)]) # temperature
Iext = numpy.array([1.6, 1.3, 1., 0.7, 0.4])  # external input current

print phi

def makecommand(direc,subdir,phi_,Iext_,JI_,Msyn_,NE_,NI_):
 cmd='python {:.2f} 4 {:.1f} 0 {:.2f} {} {} {} '.format(phi_,Iext_,JI_,Msyn_,NE_,NI_) + \
      direc + ' ' + subdir
 return cmd

##########################
# Excitatory only networks
NE=400
NI=0
direc = "excit_only"
lines = 0

#Iext and phi
subdir = "Iext_phi"
j = 0.1
m = 100
fh = open(direc+".sh", 'w')
for ii in Iext:
 for pp in phi:
  print >>fh, makecommand(direc,subdir,pp,ii,j,m,NE,NI)
  lines += 1

#Msyn and JI
subdir = "Msyn_JI"
p = 5
i = Iext[2] # Not decided yet
for mm in Msyn:
 for jj in JI:
  print >>fh, makecommand(direc,subdir,p,i,jj,mm,NE,NI)
  lines += 1

#Msyn and phi
subdir = "Msyn_phi"
j = 0.1
i = Iext[2] # Not decided yet
for mm in Msyn:
 for pp in phi:
  print >>fh, makecommand(direc,subdir,pp,i,j,mm,NE,NI)
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
# Inhibitory only networks
NE=0
NI=100
direc = "inhib_only"
lines = 0

#Iext and phi
subdir = "Iext_phi"
j = 0.1
m = 100
fh = open(direc+".sh", 'w')
for ii in Iext:
 for pp in phi:
  print >>fh, makecommand(direc,subdir,pp,ii,j,m,NE,NI)
  lines += 1

#Msyn and JI
subdir = "Msyn_JI"
p = 5
i = Iext[2] # Not decided yet
for mm in Msyn:
 for jj in JI:
  print >>fh, makecommand(direc,subdir,p,i,jj,mm,NE,NI)
  lines += 1

#Msyn and phi
subdir = "Msyn_phi"
j = 0.1
i = Iext[2] # Not decided yet
for mm in Msyn:
 for pp in phi:
  print >>fh, makecommand(direc,subdir,pp,i,j,mm,NE,NI)
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
