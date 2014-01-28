# Script to plot from a database of results.
# Adapted for the NEST simulation.
# sample input: python plotdb.py d1 d2 type
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
# output (filename text PRIMARY KEY, thres int, 
#				phi real, g real, iext real, ji real, jis real,
#				ne int, ni int, msyn int, type text, trial int)
# Expected variables: phi, iext, ji, msyn

# Plot by any two dimensions of the parameter space

print 'Creating 2D supspace...'
c.execute('DROP TABLE IF EXISTS t1')
c.execute('DROP TABLE IF EXISTS t2')
c.execute('DROP TABLE IF EXISTS t3')
c.execute('DROP TABLE IF EXISTS subspace')

if sys.argv[3] == 'both':
	ni = 100
	ne = 400
	type = 'in' # Todo: plot both in same fig
	ncells = ni

elif sys.argv[3] == 'ex':
	ni = 0
	ne = 400
	type = 'ex'
	ncells = ne

else:
	ni = 100
	ne = 0
	type = 'in'
	ncells = ni

dim1 = sys.argv[1]
dim2 = sys.argv[2]
cmd = []
if (dim1 != 'phi') and (dim2 != 'phi'):
	phi = 5. #Assume this value if phi is fixed
	cmd.append('WHERE phi='+str(phi))
if (dim1 != 'iext') and (dim2 != 'iext'):
	iext = 1.
	cmd.append('WHERE iext='+str(iext))
if (dim1 != 'ji') and (dim2 != 'ji'):
	ji = 0.1
	cmd.append('WHERE ji='+str(ji))
if (dim1 != 'msyn') and (dim2 != 'msyn'):
	msyn = 100
	cmd.append('WHERE msyn='+str(msyn))

trial = 0 # Todo: get all trials
thres = 0 # Todo: variable thres
tstart = 1000
tstop = 1300

print cmd

c.execute('CREATE TABLE IF NOT EXISTS t1 AS SELECT * FROM output '+cmd[0])
c.execute('CREATE TABLE IF NOT EXISTS t2 AS SELECT * FROM t1 '+cmd[1])
c.execute("CREATE TABLE IF NOT EXISTS subspace AS SELECT * FROM t2 WHERE \
				ni=:ni AND ne=:ne AND type=:type AND trial=:trial AND thres=:thres",
				{"ni": ni, "ne": ne, "type": type, "trial": trial, "thres": thres})

# Todo: there can be a problem when different parameter values for d1 or d2
# were tested in a different context. (The data can add more distinct values
# and mess up the table layout.) Planning to devise a workaround.

ndim1=len(c.execute('SELECT DISTINCT '+dim1+' FROM subspace').fetchall())
ndim2=len(c.execute('SELECT DISTINCT '+dim2+' FROM subspace').fetchall())
print 'Will generate', ndim1, 'by', ndim2, 'matrix of raster plots.'
print '( dim1 =', dim1, ', dim2=', dim2, ')'

print 'Reading file names...'
flist = []
tlist = []
for distinctd1 in c.execute('SELECT DISTINCT '+dim1+' FROM subspace').fetchall():
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
	
# Read from files and plot
fig = plt.figure() 
for i in range(ndim1):
	for j in range(ndim2):
#		try:
			spikes = np.loadtxt(flist[i][j].encode('ascii','ignore'),dtype='float')
			print spikes[0]
			ax = fig.add_subplot(ndim2,ndim1,ndim1*j+i+1)
			ax.scatter(spikes[:,1],spikes[:,0],s=1,c='k',marker='.')
			ax.set_title(tlist[i][j],size='6')
			ax.axis([tstart,tstop,ne,ne+50]) # Only recorded from 50 cells
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_xticks([tstart+(tstop-tstart)/5,tstart+(tstop-tstart)*2/5,tstart+(tstop-tstart)*3/5,tstart+(tstop-tstart)*4/5])
			ax.set_yticks([ne+25])
			print i, j, str(flist[i][j])
#		except:
#			print i, j, "No file yet, or something wrong with loading or plotting"
plt.show()
	
# conn.commit()
# conn.close()