import numpy as np
import pylab as plt

filename='output/inh_only/MI_JI_phi/T0mVphi0.260in100.0pA_JI-0.3_JE2.0+-0.0_E0I100_MsynI20.0_MsynE10.0_0131-2116-44_brunel-py-in-101-01.txt'

spikes = np.loadtxt(filename,dtype='float')

print spikes[:,0].tolist()

fig1 = plt.figure(num=1, figsize=(5, 3), dpi=100, facecolor='w', edgecolor='k')
plt.scatter(spikes[:,1],spikes[:,0],s=1,c='k',marker='.')
plt.axis('tight')
plt.show()

fig2 = plt.figure(num=2, figsize=(5, 3), dpi=100, facecolor='w', edgecolor='k')
plt.hist(spikes[:,0])
plt.axis('tight')
plt.show()
