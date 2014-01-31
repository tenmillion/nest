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

# Read file names and param values from DB
if not os.path.isfile('output.db'):
	print 'No output.db found in ./'
	exit()
	
# DB structure:
# CREATE TABLE output (filename text PRIMARY KEY, dir text, thres int,
# phi real, iext real, ji real, je real, jis real, ni int, ne int,
# mi real, me real, type text, trial int, kappa real, freq real)
# Expected variables: phi, (iext), ji, je, mi, me

tstart = 1000.
tstop = 1300.
binwidth = 1.
numcells = 50.
	
# Read from files and compute kappa
conn = sql.connect('output.db')
c = conn.cursor()

#c.execute('ALTER TABLE output ADD COLUMN freq real')

counter = 0
for f in glob.glob('output/*/*/*.txt'):	
	spikes = np.loadtxt(f,dtype='float')
	e_spikes = spikes[np.where(spikes[:,1]>=tstart)]
	eligible_spikes = spikes[np.where(spikes[:,1]<=tstop)]	
	freq = (1000*np.array(eligible_spikes).shape[0]/numcells)/(tstop-tstart)
	counter += 1
	print counter, freq, "Hz"
	c.execute('UPDATE output SET freq=? WHERE filename=?', (freq, f))
conn.commit()
conn.close()
print "Wrote", counter, "freqs to DB"
