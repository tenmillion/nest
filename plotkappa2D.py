# -*- coding: utf-8 -*-
# plotkappa2D.py
# Script to plot from a database of results.
# Adapted for the NEST simulation.
# sample input: python plotkappa2D.py d1 d2 type
# Y Yamamura Jan 30, 2014

import sqlite3 as sql
import sys
import os
import numpy as np
import pylab as plt

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

print 'Reading kappas...'
kappas = []
d1s = c.execute('SELECT DISTINCT '+dim1+' FROM subspace ORDER BY '+dim1+' ASC').fetchall()
d2s = c.execute('SELECT DISTINCT '+dim2+' FROM subspace ORDER BY '+dim2+' ASC').fetchall()
print "Will plot", len(d1s), "by", len(d2s), "heat map of kappas"

print dim1, d1s
print dim2, d2s

for d1 in d1s:
	ktemp = []
	for d2 in d2s:
		entry = c.execute('SELECT kappa FROM subspace WHERE '+dim1+'=? AND '+dim2+'=? \
							ORDER BY trial DESC',(d1[0],d2[0])).fetchall()
		print len(entry), "trial(s) found for", dim1, d1[0], dim2, d2[0]
		ktemp.append(np.mean(entry))
	kappas.append(ktemp)
conn.close()
npkappas = np.array(kappas)

#fig = plt.figure(num=1, figsize=(10, 7), dpi=100, facecolor='w', edgecolor='k')
#plt.rc('xtick', labelsize=5)
#plt.rc('ytick', labelsize=5)

print "Kappas averaged over trials:"
print npkappas

###

fig = plt.figure(facecolor='w', edgecolor='k')
plt.title(dim1+" vs "+dim2+" "+celltype+"("+directory+")")
plt.pcolor(npkappas)
plt.yticks(np.arange(0,len(d2s),1)+0.5,np.array(d1s)[:,0])
plt.ylabel(dim1)
plt.xticks(np.arange(0,len(d1s),1)+0.5,np.array(d2s)[:,0])
plt.xlabel(dim2)

plt.colorbar()
plt.savefig("heatmap_"+dim1+"_"+dim2+"_"+directory+"_"+celltype+".png")
plt.show()
