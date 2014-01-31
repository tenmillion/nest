# -*- coding: utf-8 -*-
# Script to compute kappa from a database of spikes and to write the values to the database.
# Adapted for the NEST simulation.
# sample input: python kappadb.py d1 d2 d3default type
# Y Yamamura Jan 29, 2014

import sqlite3 as sql
import sys
import os
import math
import numpy as np
import kappamat as km #Make sure kappamat.py is in the same directory as this

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

if sys.argv[4] == 'both':
	ni = 100
	ne = 400
	celltype = 'in' # Todo: plot both in same fig
	directory = 'both'

elif sys.argv[4] == 'ex':
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

dim1 = sys.argv[2]
dim2 = sys.argv[1]

cmd = 'WHERE'
if (dim1 != 'phi') and (dim2 != 'phi'):
	phi = 5. #Assume this value if phi is fixed
	cmd+=' phi='+str(phi)+' AND'
if (dim1 != 'iext') and (dim2 != 'iext'):
	iext = 100.
	cmd+=' iext='+str(iext)+' AND'
if (dim1 != 'ji') and (dim2 != 'ji'):
	ji = float(sys.argv[3])
	cmd+=' ji='+str(ji)+' AND'
if (dim1 != 'je') and (dim2 != 'je'):
	je = 2.
	cmd+=' je='+str(je)+' AND'
if (dim1 != 'mi') and (dim2 != 'mi'):
	mi = float(sys.argv[3])
	cmd+=' mi='+str(mi)+' AND'
if (dim1 != 'me') and (dim2 != 'me'):
	me = 10.
	cmd+=' me='+str(me)+' AND'
cmd+=' thres='+str(thres)

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
for distinctd1 in c.execute('SELECT DISTINCT '+dim1+' FROM subspace').fetchall():
	ftemp = []
	for distinctd2 in c.execute('SELECT DISTINCT '+dim2+' FROM subspace').fetchall():
		ftemp2=[]
		entries = c.execute('SELECT filename, trial FROM subspace WHERE '+dim1+'='+str(distinctd1[0])+\
							' AND '+dim2+'='+str(distinctd2[0])+' ORDER BY trial DESC').fetchall()
		print "found", len(entries), "trials"
		for entry in entries:
			ftemp2.append(entry[0])
			#print "current trial:", entry[1]
		ftemp.append(ftemp2)
	flist.append(ftemp)
		#print "d1:",distinctd1
		#print "flist now:", flist
	
print "the file list size:", np.array(flist).shape
	
# Read from files and compute kappa
counter = 0
for i in range(ndim1):
	row = []
	for j in range(ndim2):
		for trial in range(len(flist[i][j])):
			try:
				print i, j, trial, "Loading"
				spikes = np.loadtxt(flist[i][j][trial].encode('ascii','ignore'),dtype='float')
				#print spikes[0], "<-first row"
			except:
				print i, j, trial, "No file yet, or something wrong with loading"
				spikes = []
				print spikes
			k = km.kappa(spikes,binwidth,tstart,tstop)
			c.execute('UPDATE output SET kappa=? WHERE filename=?', (k, flist[i][j][trial]))
			#print k
			counter += 1
conn.commit()
conn.close()
print "Wrote", counter, "kappas to DB"
