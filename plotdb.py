# Script to plot from a database of results.
# Adapted for the NEST simulation.
# sample input: python plotdb.py d1 d2 celltype
# Y Yamamura Jan 29, 2014

import sqlite3 as sql
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# Read file names and param values from DB
if not os.path.isfile('output.db'):
	print 'No output.db found in ./'
	exit()
	
conn = sql.connect('output.db')
c = conn.cursor()

# DB structure:
# CREATE TABLE output (filename text PRIMARY KEY, dir text, thres int,
# phi real, iext real, ji real, je real, jis real, ni int, ne int,
# mi real, me real, type text, trial int, kappa real)
# Expected variables: phi, (iext), ji, je, mi, me

# Plot by any two dimensions of the parameter space
# Parameter range must be rectangular

print 'Creating 2D supspace...'
c.execute('DROP TABLE IF EXISTS t1')
c.execute('DROP TABLE IF EXISTS subspace')

if sys.argv[3] == 'both':
	ni = 100
	ne = 400
	celltype = 'in' # Todo: plot both in same fig
	directory = 'both'

elif sys.argv[3] == 'ex':
	ni = 0
	ne = 400
	celltype = 'ex'
	directory = 'excit_only'

else:
	ni = 100
	ne = 0
	celltype = 'in'
	directory = 'inh_only'

thres = 0 # Todo: variable thres

dim1 = sys.argv[1]
dim2 = sys.argv[2]

cmd = 'WHERE'
if (dim1 != 'phi') and (dim2 != 'phi'):
	phi = 5. #Assume this value if phi is fixed
	cmd+=' phi='+str(phi)+' AND'
if (dim1 != 'iext') and (dim2 != 'iext'):
	iext = 100.
	cmd+=' iext='+str(iext)+' AND'
if (dim1 != 'ji') and (dim2 != 'ji'):
	ji = 5.
	cmd+=' ji='+str(ji)+' AND')
if (dim1 != 'je') and (dim2 != 'je'):
	je = 2.
	cmd+=' ji='+str(je)+' AND'
if (dim1 != 'mi') and (dim2 != 'mi'):
	mi = 25.
	cmd+=' mi='+str(mi)+' AND'
if (dim1 != 'me') and (dim2 != 'me'):
	me = 10.
	cmd+=' me='+str(me)+' AND'
cmd+='thres='+str(thres)

tstart = 1000
tstop = 1300
binwidth = 1.

print cmd

c.execute('CREATE TABLE IF NOT EXISTS t1 AS SELECT * FROM output '+cmd)
c.execute("CREATE TABLE IF NOT EXISTS subspace AS SELECT * FROM t1 WHERE \
				dir=:directory AND ni=:ni AND ne=:ne AND type=:type",
				{"directory":directory, "ni": ni, "ne": ne, "type": celltype})
## If params are modified, modify down to here.

ndim1=len(c.execute('SELECT DISTINCT '+dim1+' FROM subspace').fetchall())
ndim2=len(c.execute('SELECT DISTINCT '+dim2+' FROM subspace').fetchall())
print 'Will generate', ndim1, 'by', ndim2, 'matrix of raster plots.'
print '( dim1 =', dim1, ', dim2=', dim2, ')'

print 'Reading file names...'
flist = []
tlist = []
for distinctd1 in c.execute('SELECT DISTINCT '+dim1+' FROM subspace ORDER BY '+dim1+' ASC').fetchall():
	ftemp = []
	ttemp = []
	for entry in c.execute('SELECT filename, '+dim1+', '+dim2+' FROM subspace WHERE '+dim1+'='+str(distinctd1[0])+' ORDER BY '+dim2+' DESC'):
		ftemp.append(entry[0])
		ttemp.append(dim1+'='+str(entry[1])+', '+dim2+'='+str(entry[2]))
		#print "current entry:", entry
	flist.append(ftemp)
	tlist.append(ttemp)
	#print "d1:",distinctd1
	#print "flist now:", flist
	#print "tlist now:", tlist
	
#print "flist now:", flist
print "tlist now:"
for column in tlist:
	for tuple in column:
		print tuple

print len(flist), len(flist[0])
	
# Read from files and plot
fig = plt.figure(num=1, figsize=(10, 7), dpi=100, facecolor='w', edgecolor='k')
plt.rc('xtick', labelsize=5)
plt.rc('ytick', labelsize=5)
for i in range(ndim1):
	for j in range(ndim2):
#		try
			print i, j, "Load"
			spikes = np.loadtxt(flist[i][j].encode('ascii','ignore'),dtype='float')
			print spikes[0]
			ax = fig.add_subplot(ndim2,ndim1,ndim1*j+i+1)
			ax.scatter(spikes[:,1],spikes[:,0],s=1,c='k',marker='.')
			ax.set_title(tlist[i][j],size='6')
			if celltype=='ex': # Plot exc
				ax.axis([tstart,tstop,0,50])			# Only recorded from 50 cells
				ax.set_yticklabels([0,50])
				ax.set_yticks([0,50])
			else: # Plot inh
				ax.axis([tstart,tstop,ne,ne+50])	# Only recorded from 50 cells
				ax.set_yticklabels([ne,ne+50])
				ax.set_yticks([ne,ne+50])
			ax.set_xticklabels([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])
			ax.set_xticks([tstart,tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5,tstop])
			print i, j, str(flist[i][j])
#		except:
#			print i, j, "No file yet, or something wrong with loading or plotting"
plt.suptitle(dim1+" vs "+dim2+" "+celltype+"("+directory+")")
plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=0.9)
plt.savefig(dim1+'_'+dim2+'_'+celltype+'_'+directory+'.png')
plt.show()
# conn.commit()
conn.close()
