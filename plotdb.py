# Script to plot from a database of results.
# Only supports rectangular, evenly-spaced parameter subspaces for now.
# Considering adding feature to allow arbitrary spacing.
# Y Yamamura Oct 8-9, 2013

dir = 'Data/'

import sqlite3 as sql
import os
import matplotlib.pyplot as plt
import numpy as np

# Read file names and param values from DB
if not os.path.isfile(dir+'output.db'):
	print 'No output.db found in',dir
	exit()
	
conn = sql.connect(dir+'output.db')
c = conn.cursor()

# DB structure:
# output (filename text PRIMARY KEY,
# temp real, inp real, inb real, pp real, pb real, bp real )

# Plot by any two dimensions of the parameter space
# Todo: specify params from command line

dim1 = 'temp'
dim2 = 'bp'
inp = 0.015
inb = 0.015
#temp = 30
pp = 1
pb = 0.5
#bp = 1
nmda = 0.2

tstop = 500
ncells = 500

print 'Creating 2D supspace...'
c.execute('DROP TABLE IF EXISTS t1')
c.execute('DROP TABLE IF EXISTS t2')
c.execute('DROP TABLE IF EXISTS t3')
c.execute('DROP TABLE IF EXISTS t4')
c.execute('DROP TABLE IF EXISTS subspace')

#---Only need to change this part------
#c.execute('CREATE TABLE t1 AS SELECT * FROM output WHERE temp='+str(temp))
c.execute('CREATE TABLE IF NOT EXISTS t1 AS SELECT * FROM output WHERE inp='+str(inp))
c.execute('CREATE TABLE IF NOT EXISTS t2 AS SELECT * FROM t1 WHERE inb='+str(inb))
c.execute('CREATE TABLE IF NOT EXISTS t3 AS SELECT * FROM t2 WHERE pp='+str(pp))
c.execute('CREATE TABLE IF NOT EXISTS t4 AS SELECT * FROM t3 WHERE nmda='+str(nmda))
c.execute('CREATE TABLE IF NOT EXISTS subspace AS SELECT * FROM t4 WHERE pb='+str(pb))
#c.execute('CREATE TABLE subspace AS SELECT * FROM t3 WHERE bp='+str(bp))
#--------------------------------------

# Todo: there can be a problem when different parameter values for d1 or d2
# were tested in a different context. (The data can add more distinct values
# and mess up the table layout.) Planning to devise a workaround.

ndim1=len(c.execute('SELECT DISTINCT '+dim1+' FROM subspace').fetchall())
ndim2=len(c.execute('SELECT DISTINCT '+dim2+' FROM subspace').fetchall())
print 'Will generate', ndim1, 'by', ndim2, 'matrix of scatter plots.'
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
		try:
			spikes = np.transpose(np.loadtxt(dir+str(flist[i][j]),dtype='float',skiprows=1,delimiter=' ',usecols=(1,2)))
			ax = fig.add_subplot(ndim2,ndim1,ndim1*j+i+1)
			ax.scatter(spikes[0],spikes[1],s=1,c='k',marker='.')
			ax.set_title(tlist[i][j],size='6')
			ax.axis([0,tstop,0,ncells])
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			ax.set_xticks([tstop/5,tstop*2/5,tstop*3/5,tstop*4/5])
			ax.set_yticks([ncells*4/5])
			print i, j, str(flist[i][j])
		except:
			print i, j, "No file yet"
plt.show()
	
# conn.commit()
# conn.close()