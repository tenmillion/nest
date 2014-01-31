from pylab import *
import numpy as np
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

dim1=sys.argv[1]
dim2=sys.argv[2] #expects phi
celltype=sys.argv[4]
directory=sys.argv[4]
d2s = np.array([5*3**((k-37)/10.) for k in np.arange(10,41,3)])

mat1=np.loadtxt('allphi_'+dim1+'_'+dim2+'_'+str(sys.argv[3])+'.txt')
mat2=np.loadtxt('lowphi_'+dim1+'_'+dim2+'_'+str(sys.argv[3])+'.txt')

npkappas=np.concatenate((mat2,mat1[:,mat2.shape[1]:]),axis=1)

fig, ax = plt.subplots()
fig = plt.gcf()
title("concatenated "+dim1+" vs "+dim2+' other dim value:'+str(sys.argv[3])+" "+celltype+"("+directory+")")
heatmap = ax.pcolor(npkappas)
ax.axis('tight')
ylabel(dim1)
xticks(np.arange(0,len(d2s),1)+0.5,np.round(d2s,1))
xlabel(dim2)
ax.set_aspect('equal')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(heatmap,cax=cax)
savefig("figures/"+"concat_heatmap_"+dim1+"_"+dim2+'_'+str(sys.argv[3])+"_"+directory+"_"+celltype+".png", facecolor='w')
show()
