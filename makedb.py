import sqlite3 as sql
import os
import glob
import re

# makedb.py modified for the NEST simulations

# Create table
if not os.path.isfile('output.db'):
	conn = sql.connect('output.db')
	c = conn.cursor()
	c.execute('''CREATE TABLE output (filename text PRIMARY KEY, thres int, phi real, g real, iext real, ji real, jis real, ne int, ni int, msyn int, type text, trial int)''')
	print "Created table"
else:
	conn = sql.connect('output.db')
	c = conn.cursor()
	
# Read filenames from directory and add table entries
# expected format of output file names:
# "T0mVphi0.6g4.0in40.0pA_JI-0.100+-0.0_E400I0_MsynII100_0128-1413-41_brunel-py-ex-401-0.txt"
#     thres    phi,   g,    iext,            ji,      jis    ne,     ni,            msyn,                                     type,  trial

nadded = 0

for x in glob.glob('output/*/*/*.txt'):
	thres_x = re.search('T([0-9]+)mV',x).group(1)
	phi_x = re.search('phi([0-9]+\.[0-9]+)',x).group(1)		#
	g_x = re.search('g([0-9]+\.[0-9]+)',x).group(1)
	iext_x = re.search('in([0-9]+\.[0-9]+)pA',x).group(1)	#
	ji_x = re.search('JI\-([0-9]+\.[0-9]+)',x).group(1)			#
	jis_x = re.search('\+\-([0-9]+\.[0-9]+)',x).group(1)
	ne_x = re.search('E([0-9]+)',x).group(1)						#
	ni_x = re.search('I([0-9]+)',x).group(1)						#
	msyn_x = re.search('MsynII([0-9]+)',x).group(1)			#
	type_x = re.search('brunel\-py\-(ex|in)',x).group(1)		#
	trial_x = re.search('([0-9]+)\.txt',x).group(1)				#
	entry = [(x.replace("\\","/"), thres_x, phi_x, g_x, iext_x, ji_x, jis_x, ne_x, ni_x, msyn_x, type_x, trial_x),]
	print entry
	
	c.execute('SELECT COUNT(filename) FROM output WHERE thres=:thres AND phi=:phi AND g=:g AND iext=:iext AND ji=:ji AND\
						jis=:jis AND ne=:ne AND ni=:ni AND msyn=:msyn AND type=:type AND trial=:trial',\
						{"thres":thres_x, "phi":phi_x, "g":g_x, "iext":iext_x, "ji":ji_x, "jis":jis_x, "ne":ne_x, "ni":ni_x, "msyn":msyn_x, "type":type_x, "trial":trial_x})
	entry=c.fetchone()[0]

#	try:
	if entry==0:
		c.executemany('INSERT INTO output VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', entry)
		nadded += 1
	else:
		print "Overlapping parameters", x
#	except:
#		print "Did not add duplicate", x

conn.commit()
print "Added", nadded, "files"
	
# Print table just created, ordered
print "Entries now in table: "
for row in c.execute('SELECT * FROM output ORDER BY type, ne, ni, iext, ji, msyn, phi, trial'):
	print row
	
conn.close()