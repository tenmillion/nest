import sqlite3 as sql
import os
import glob
import re

# makedb.py modified for the NEST simulations

# Create table
if not os.path.isfile('output.db'):
	conn = sql.connect('output.db')
	c = conn.cursor()
	c.execute('''CREATE TABLE output (filename text PRIMARY KEY, dir text, thres int, phi real, iext real, ji real, je real, jis real, ni int, ne int, mi real, me real, type text, trial text, kappa real, freq real)''')
	print "Created table"
else:
	conn = sql.connect('output.db')
	c = conn.cursor()
	
# Read filenames from directory and add table entries
# expected format of output file names:
# /home/yoriko-y/Dropbox/Academic/NAIST_MI/Cool/_MODELS/Nest/output/
# inh_only/JI_phi/T0mVphi6.950in100.0pA_JI-0.1_JE2.0+-0.0_E0I100_MsynI25.0_MsynE10.0_0130-1540-08_brunel-py-in-101-0.txt
# dir             thres    phi,   iext,  ji,   je   jis    ne, ni,  mi,     me                              type,  trial

nadded = 0
nskipped = 0
ndups = 0

for x in glob.glob('output/*/*/*.gdf'):
	os.rename(x,x[:-3]+'txt')

for x in glob.glob('output/*/*/*.txt'):
	x = x.replace("\\","/")
	dir_x = re.search('output/([a-z]+_?[a-z]*)/',x).group(1)			#
	thres_x = re.search('T([0-9]+)mV',x).group(1)
	phi_x = re.search('phi([0-9]+\.[0-9]+)',x).group(1)		#
	iext_x = re.search('in([0-9]+\.[0-9]+)pA',x).group(1)	#
	ji_x = re.search('JI\-([0-9]+\.[0-9]+)',x).group(1)			#
	je_x = re.search('JE([0-9]+\.[0-9]+)',x).group(1)			#
	jis_x = re.search('\+\-([0-9]+\.[0-9]+)',x).group(1)
	ni_x = re.search('I([0-9]+)',x).group(1)						#
	ne_x = re.search('_E([0-9]+)',x).group(1)						#
	mi_x = re.search('MsynI([0-9]+\.[0-9]+)',x).group(1)			#
	me_x = re.search('MsynE([0-9]+\.[0-9]+)',x).group(1)			#
	type_x = re.search('brunel\-py\-(ex|in)',x).group(1)		#
	trial_x = re.search('_([0-9]+\-[0-9]+\-[0-9]+)_',x).group(1)				#
	entry = [(x, dir_x, thres_x, phi_x, iext_x, ji_x, je_x, jis_x, ni_x, ne_x, mi_x, me_x, type_x, trial_x),]
	#print entry
	
	# Check for existing entries with different file name but same params
	c.execute('SELECT COUNT(*) FROM output WHERE dir=:dir AND thres=:thres AND phi=:phi AND iext=:iext AND ji=:ji AND je=:je AND\
						jis=:jis AND ni=:ni AND ne=:ne AND mi=:mi AND me=:me AND type=:type AND trial=:trial',\
						{"dir": dir_x, "thres":thres_x, "phi":phi_x, "iext":iext_x, "ji":ji_x, "je":je_x, "jis":jis_x, "ni":ni_x, \
						 "ne": ne_x, "mi":mi_x, "me":me_x, "type":type_x, "trial":trial_x})
	count=c.fetchone()[0]
	c.execute('SELECT * FROM output WHERE dir=:dir AND thres=:thres AND phi=:phi AND iext=:iext AND ji=:ji AND je=:je AND\
	jis=:jis AND ni=:ni AND ne=:ne AND mi=:mi AND me=:me AND type=:type AND trial=:trial',\
	{"dir": dir_x, "thres":thres_x, "phi":phi_x, "iext":iext_x, "ji":ji_x, "je":je_x, "jis":jis_x, "ni":ni_x, \
						 "ne": ne_x, "mi":mi_x, "me":me_x, "type":type_x, "trial":trial_x})
	overlap=c.fetchall()
	
	if count==0:
		c.executemany('INSERT INTO output VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,-1,-1)', entry)
		nadded += 1
	else:
		print "Overlapping parameters", x
		print overlap, count
		nskipped += count

conn.commit()
	
# Print table just created, ordered
#print "Entries now in table: "
#for row in c.execute('SELECT * FROM output ORDER BY type, ni, ne, iext, ji, je, mi, me, phi, trial, kappa'):
#	print row
	
conn.close()

print "Added", nadded, "files"
print "Skipped", nskipped, "files due to overlapping parameters"
