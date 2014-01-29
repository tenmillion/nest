#plotkappa.py
import numpy as np
import pylab as plt

kappas=np.loadtxt("phiiext.txt")
fig = plt.figure()
plt.pcolor(kappas.transpose())
plt.yticks(np.arange(0,5,1)+0.5,(160,130,100,70,40))
plt.ylabel('iext')
plt.xticks(np.arange(0,5,1)+0.5,(0.6,1.0,1.7,2.9,5.0))
plt.xlabel('phi')


plt.colorbar()
plt.show()