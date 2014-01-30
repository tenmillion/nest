# -*- coding: utf-8 -*-
# Script to compute kappa from a database of spikes and to write the values to the database.
# Adapted for the NEST simulation.
# sample input: python kappadb.py d1 d2 type
# Y Yamamura Jan 29, 2014

import sqlite3 as sql
import sys
import os
import numpy as np
import kappamat as km #Make sure kappamat.py is in the same directory as this

# Read file names and param values from DB
if not os.path.isfile('output.db'):
	print 'No output.db found in ./'
	exit()
	
conn = sql.connect('output.db')
c = conn.cursor()

#testkappa = c.execute('PRAGMA index_info(kappa)')		  # Check for column kappa in output.db
#nokappa = testkappa.fetchall()
#print "kappa exists?", nokappa # Should be empty list if kappa doesn't exist yet
#if nokappa) == 0:
#	c.execute('ALTER TABLE output ADD COLUMN kappa real') # Add kappa column if it doesn't exist
#	print "Added new column kappa to output"
#	conn.commit()

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
	type = 'in' # Todo: plot both in same fig
	directory = 'both'

elif sys.argv[3] == 'ex':
	ni = 0
	ne = 400
	type = 'ex'
	directory = 'excit_only'

else:
	ni = 100
	ne = 0
	type = 'in'
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

trial = 0 # Todo: get all trials
thres = 0 # Todo: variable thres
tstart = 0
tstop = 1300
binwidth = 1.

print cmd

c.execute('CREATE TABLE IF NOT EXISTS t1 AS SELECT * FROM output '+cmd[0])
c.execute('CREATE TABLE IF NOT EXISTS t2 AS SELECT * FROM t1 '+cmd[1])
c.execute("CREATE TABLE IF NOT EXISTS subspace AS SELECT * FROM t2 WHERE \
				dir=:directory AND ni=:ni AND ne=:ne AND type=:type AND trial=:trial AND thres=:thres",
				{"directory":directory, "ni": ni, "ne": ne, "type": type, "trial": trial, "thres": thres})

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

print len(flist), len(flist[0])
	
# Read from files and compute kappa
kappas = []
for i in range(ndim1):
	row = []
	for j in range(ndim2):
		try:
			print i, j, "Loading"
			spikes = np.loadtxt(flist[i][j].encode('ascii','ignore'),dtype='float')
			print spikes[0], "<-first row"
		except:
			print i, j, "No file yet, or something wrong with loading"
		k = km.kappa(spikes,binwidth,tstart,tstop)
		c.execute('UPDATE output SET kappa=? WHERE filename=?', (k, flist[i][j]))
		row.append(k)
		print k
	kappas.append(row)
print kappas
filename = dim1+dim2+".txt"
np.savetxt(filename,kappas)
print "Saved", np.shape(kappas), "kappa matrix to ./"+filename
conn.commit()
conn.close()
print "Wrote kappas to DB"
