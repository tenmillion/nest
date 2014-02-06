import numpy as np
import pylab as plt
import glob
import re

mixed = False

#filename = glob.glob('output/*/*/*.txt')[0]
filename1 = 'output/inh_only/MI_JI_i10/T0mVphi5.000in100.0pA_JI-10.0_JE2.0+-0.0_E0I100_MsynI0.0_MsynE10.0_0204-0539-15_brunel-py-in-101-0.txt'
filename2 = '/home/yoriko-y/Dropbox/Academic/NAIST_MI/Cool/_MODELS/Nest/output/mixed/M_JE_i10/T0mVphi0.498in100.0pA_JI-10.0_JE100.0+-0.0_E400I100_MsynI40.0_MsynE40.0_0204-0716-34_brunel-py-ex-501-0.txt'

spikes1 = np.loadtxt(filename1,dtype='float')
print spikes1[:,0].tolist()
if mixed:
 spikes2 = np.loadtxt(filename2,dtype='float')
 print spikes2[:,0].tolist()

tstart=1000
tstop=1300

left = 0.1
width = 0.8
bottom_spikes = 0.3
height_spikes = 0.6

mi = re.search('MsynI([0-9]+\.[0-9]+)',filename1).group(1)	
ji = re.search('JI\-([0-9]+\.[0-9]+)',filename1).group(1)
if mixed:
 je = re.search('JE([0-9]+\.[0-9]+)',filename2).group(1)
phi = re.search('phi([0-9]+\.[0-9]+)',filename1).group(1)

plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
if mixed:
 ax1 = plt.axes([left, 0.49, width, 0.4])
 ax2 = plt.axes([left, 0.09, width, 0.4])
else:
 ax1 = plt.axes([left, bottom_spikes, width, height_spikes])
 ax2 = plt.axes([left, 0.05, width, 0.2])

if mixed:
	plt.suptitle(r'$\phi='+phi+', M_{syn}^I='+mi+', J^I='+ji+', J^E='+je+'$')
else:
#	plt.suptitle(r'$\phi='+phi+', M_{syn}^I='+mi+', J^I='+ji+'$')
	plt.suptitle(r'$\mathrm{Temp}=\mathrm{BL}-0^{\circ}C$, $M_{syn}='+mi+'$')

ax1.scatter(spikes1[:,1],spikes1[:,0],s=1,c='k',marker='.')
ax1.axis([tstart,tstop,np.min(spikes1[:,0]),np.max(spikes1[:,0])])
ax1.set_yticklabels([0,'','','',str(np.max(spikes1[:,0]))],size=6)
ax1.set_yticks(np.arange(np.min(spikes1[:,0]),np.max(spikes1[:,0])+1,10))
ax1.set_xticklabels([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop],size=6)
ax1.set_xticks([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])
ax1.set_ylabel("Neuron ID")

if mixed:
	ax2.scatter(spikes2[:,1],spikes2[:,0],s=1,c='k',marker='.')
	ax2.axis([tstart,tstop,np.min(spikes2[:,0]),np.max(spikes2[:,0])])			# Only recorded from 50 cells
	ax2.set_yticklabels([0,'','','','',str(np.max(spikes2[:,0]))],size=6)
	ax2.set_yticks(np.arange(np.min(spikes2[:,0]),np.max(spikes2[:,0])+1,10))
	ax2.set_xticklabels([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop],size=6)
	ax2.set_xticks([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])
	ax2.set_ylabel("Excitatory")
else:
	ax2.hist(spikes1[:,1], bins=np.arange(np.min(spikes1[:,1]),np.max(spikes1[:,1]),1), color='0.5')
	ax2.axis([tstart,tstop,0,20])
	ax2.set_yticklabels(np.arange(0,np.max(spikes1[:,0])+1,10),size=6)
	ax2.set_yticks(np.arange(0,np.max(spikes1[:,0])+1,10))
	ax2.set_xticklabels([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop],size=6)
	ax2.set_xticks([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])

plt.savefig('../../../Thesis/-0_'+mi+'.eps')
plt.show()

filename = filename1.replace('brunel-py-in','voltmeter')
filename = filename.replace('03-0','04-0')
filename = filename.replace('txt','dat')
vtrace = np.loadtxt(filename,dtype='float')
plt.figure(num=2, facecolor='w')
plt.title('Membrane voltage trace for '+r'$\phi='+phi+', M_{syn}^I='+mi+', J^I='+ji+'$')
plt.ylabel('Membrane voltage [mV]')
plt.xlabel('Time [msec]')
plt.plot(vtrace[:,1],vtrace[:,2],color='0.1')
plt.show()
