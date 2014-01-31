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
import glob
import kappamat as km #Make sure kappamat.py is in the same directory as this

# Read file names and param values from DB
if not os.path.isfile('output.db'):
	print 'No output.db found in ./'
	exit()
	
# DB structure:
# CREATE TABLE output (filename text PRIMARY KEY, dir text, thres int,
# phi real, iext real, ji real, je real, jis real, ni int, ne int,
# mi real, me real, type text, trial int, kappa real)
# Expected variables: phi, (iext), ji, je, mi, me

tstart = 1000
tstop = 1300
binwidth = 1.
	
# Read from files and compute kappa
conn = sql.connect('output.db')
c = conn.cursor()

counter = 0
for f in glob.glob('output/*/*/*.txt'):	
	spikes = np.loadtxt(f,dtype='float')
	k = km.kappa(spikes,binwidth,tstart,tstop)
	if math.isnan(k):
		k = -1
		print "Dims of spike train matrix",spikes.shape
		print f[15:-15],"k corrected to -1"
	counter += 1
	print counter, k
	c.execute('UPDATE output SET kappa=? WHERE filename=?', (k, f))
conn.commit()
conn.close()
print "Wrote", counter, "kappas to DB"
