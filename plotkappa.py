# -*- coding: utf-8 -*-
# plotkappa.py
# Script to plot overlapping line graphs from a database of results.
# Adapted for the NEST simulation.
# Almost the exact same as plotkappa2D...
# sample input: python plotkappa.py d1 d2 type
# Y Yamamura Jan 30, 2014

#plotkappa2D.py
import numpy as np
import sqlite3 as sql
import sys
import os
import pylab as plt

# Read file names and param values from DB
if not os.path.isfile('output.db'):
	print 'No output.db found in ./'
	exit()
	
conn = sql.connect('output.db')
c = conn.cursor()

testkappa = c.execute('PRAGMA index_info(kappa)')		   # Check for column kappa in output.db
nokappa = testkappa.fetchall()
print nokappa # Should be empty list if kappa doesn't exist yet
if nokappa:
	print "Please compute kappa first."
	exit()

# DB structure:
# CREATE TABLE output (filename text PRIMARY KEY, dir text, thres int,
# phi real, iext real, ji real, je real, jis real, ni int, ne int,
# mi real, me real, type text, trial int, kappa real)
# Expected variables: phi, (iext), ji, je, mi, me

# Plot by any two dimensions of the parameter space
# Parameter range must be rectangular

print 'Creating 2D supspace...'
c.execute('DROP TABLE IF EXISTS t1')
c.execute('DROP TABLE IF EXISTS t2')
c.execute('DROP TABLE IF EXISTS t3')
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

dim1 = sys.argv[1]
dim2 = sys.argv[2]
cmd = []
if (dim1 != 'phi') and (dim2 != 'phi'):
	phi = 5. #Assume this value if phi is fixed
	cmd.append('WHERE phi='+str(phi))
if (dim1 != 'iext') and (dim2 != 'iext'):
	iext = 100.
	cmd.append('WHERE iext='+str(iext))
if (dim1 != 'ji') and (dim2 != 'ji'):
	ji = 0.1
	cmd.append('WHERE ji='+str(ji))
if (dim1 != 'msyn') and (dim2 != 'msyn'):
	msyn = 100
	cmd.append('WHERE msyn='+str(msyn))

thres = 0 # Todo: variable thres
tstart = 1000
tstop = 1300
binwidth = 1.

print cmd

c.execute('CREATE TABLE IF NOT EXISTS t1 AS SELECT * FROM output '+cmd[0])
c.execute('CREATE TABLE IF NOT EXISTS t2 AS SELECT * FROM t1 '+cmd[1])
c.execute("CREATE TABLE IF NOT EXISTS subspace AS SELECT * FROM t2 WHERE \
				dir=:directory AND ni=:ni AND ne=:ne AND type=:type AND thres=:thres",
				{"directory":directory, "ni": ni, "ne": ne, "type": celltype, "thres": thres})

print 'Reading kappas...'
kappas = []
d1s = c.execute('SELECT DISTINCT '+dim1+' FROM subspace ORDER BY '+dim1+' DESC').fetchall()		# Display higher param at top
d2s = c.execute('SELECT DISTINCT '+dim2+' FROM subspace ORDER BY '+dim2+' ASC').fetchall()
print "Will plot line graphs for ", len(d1s), dim1, "values", len(d2s), "data points in", dim2

print dim1, d1s
print dim2, d2s

for d1 in d1s:
	ktemp = []
	for d2 in d2s:
		entry = c.execute('SELECT kappa FROM subspace WHERE '+dim1+'=? AND '+dim2+'=? \
							ORDER BY trial DESC',(d1[0],d2[0])).fetchall()
		print len(entry), "trial(s) found for", dim1, d1[0], dim2, d2[0]
		print entry
		ktemp.append(np.mean(entry))
	kappas.append(ktemp)
conn.close()
npkappas = np.array(kappas)

print "Kappas averaged over trials:"
print npkappas

cmds = [(np.array(d2s)[:,0].tolist(), npkappas[i].tolist()) for i in range(len(d1s))]

fig = plt.figure(facecolor='w', edgecolor='k')
ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
labels = [dim1+"="+str(np.array(d1s)[i,0]) for i in range(len(d1s))]
handles = []
i = 0
for cmd in cmds:
	i +=1
	axis1,axis2 = cmd
	handle,=ax.plot(axis1,axis2,'o-')
	handles.append(handle)
plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95)
ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title(dim1+" vs "+dim2+" "+celltype+"("+directory+")")
plt.ylabel("Kappa")
plt.xlabel(dim2)

plt.savefig("lines_"+dim1+"_"+dim2+"_"+directory+"_"+celltype+".png")
plt.show()
