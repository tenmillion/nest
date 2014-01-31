import numpy as np
import pylab as plt
import glob
import re

#filename = glob.glob('output/*/*/*.txt')[0]
filename = 'output/inh_only/MI_JI_phi/T0mVphi1.340in100.0pA_JI-0.6_JE2.0+-0.0_E0I100_MsynI100.0_MsynE10.0_0131-2323-22_brunel-py-in-101-0.txt'
spikes = np.loadtxt(filename,dtype='float')

print spikes[:,0].tolist()

tstart=0
tstop=1300

left = 0.1
width = 0.8
bottom_spikes = 0.3
height_spikes = 0.6

mi = re.search('MsynI([0-9]+\.[0-9]+)',filename).group(1)	
ji = re.search('JI\-([0-9]+\.[0-9]+)',filename).group(1)
phi = re.search('phi([0-9]+\.[0-9]+)',filename).group(1)

plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
ax1 = plt.axes([left, bottom_spikes, width, height_spikes])
ax2 = plt.axes([left, 0.05, width, 0.2])

plt.suptitle(r'$\phi='+phi+', M_{syn}^I='+mi+', J^I='+ji+'$')
ax1.scatter(spikes[:,1],spikes[:,0],s=1,c='k',marker='.')
ax1.axis([tstart,tstop,0,50])			# Only recorded from 50 cells
ax1.set_yticklabels([0,'','','','',str(np.max(spikes[:,0]))],size=6)
ax1.set_yticks(np.arange(0,np.max(spikes[:,0])+1,10))
ax1.set_xticklabels([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop],size=6)
ax1.set_xticks([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])

ax2.hist(spikes[:,1], bins=np.arange(np.min(spikes[:,1]),np.max(spikes[:,1]),10), color='0.5')
ax2.axis('tight')
ax2.set_yticklabels(np.arange(0,np.max(spikes[:,0])+1,10),size=6)
ax2.set_yticks(np.arange(0,np.max(spikes[:,0])+1,10))
ax2.set_xticklabels([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop],size=6)
ax2.set_xticks([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])

plt.show()

filename = filename.replace('brunel-py-in','voltmeter')
filename = filename.replace('-101-','-102-')
filename = filename.replace('txt','dat')
vtrace = np.loadtxt(filename,dtype='float')
plt.figure(num=2, facecolor='w')
plt.title('Membrane voltage trace for '+r'$\phi='+phi+', M_{syn}^I='+mi+', J^I='+ji+'$')
plt.ylabel('Membrane voltage [mV]')
plt.xlabel('Time [msec]')
plt.plot(vtrace[:,1],vtrace[:,2],color='0.1')
plt.show()
