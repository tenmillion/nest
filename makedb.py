import sqlite3 as sql
import os
import glob
import re


# Create table
if not os.path.isfile('output.db'):
	conn = sql.connect('output.db')
	c = conn.cursor()
	c.execute('''CREATE TABLE output (filename text PRIMARY KEY, temp real, inp real, inb real, pp real, pb real, bp real)''')
	print "Created table"
else:
	conn = sql.connect('output.db')
	c = conn.cursor()
	
# Read filenames from directory and add table entries
# expected format of output file names:
# "s-degc%3.1f-inp%3.3f-inb%3.3f-pp%3.5f-pb%3.5f-bp%3.5f-.dat",
#     celsius, inp,     inb, p2pscale, p2bscale, b2pscale

nadded = 0

for x in glob.glob('s-*.dat'):
	temp_x = re.search('degc([0-9]+\.[0-9]+)',x).group(1)
	inp_x = re.search('inp([0-9]+\.[0-9]+)',x).group(1)
	inb_x = re.search('inb([0-9]+\.[0-9]+)',x).group(1)
	pp_x = re.search('pp([0-9]+\.[0-9]+)',x).group(1)
	pb_x = re.search('pb([0-9]+\.[0-9]+)',x).group(1)
	bp_x = re.search('bp([0-9]+\.[0-9]+)',x).group(1)

	entry = [(x, temp_x, inp_x, inb_x, pp_x, pb_x, bp_x),]
	
	try:
		c.executemany('INSERT INTO output VALUES (?,?,?,?,?,?,?)', entry)
		nadded += 1
	except:
		print "Did not add duplicate "+x

conn.commit()
print "Added", nadded, "files"
	
# Print table just created, ordered by temperature, then by inp and inb, then by pp, pb, bp
print "Entries now in table: "
for row in c.execute('SELECT * FROM output ORDER BY temp, inp, inb, pp, pb, bp'):
	print row
	
conn.close()