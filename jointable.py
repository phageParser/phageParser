#!/usr/bin/python
import sqlite3

# THIS file will parse in two files named 'DRdatabase' and 'Spacerdatabase', 
# and populate the 'crispr.sqlite' file tables Spacers and Repeats, as per 
# ISSUE #99 https://github.com/goyalsid/phageParser/issues/99

def printDict(d):
	for i in d:
		print i,d[i]



def SQL_add(elem, name, dbName, columns):
	conn = sqlite3.connect(dbName)
	c = conn.cursor()
	q = "INSERT INTO " + str(name) + " " + columns + " VALUES "+ elem+""
	#print q
	try:
		c.execute(q)
		conn.commit()
	except sqlite3.IntegrityError:
		print "Already In Database", elem
	c.close()
	conn.close()


def SQL_search(elem, name, column, dbName):
	conn = sqlite3.connect(dbName)
	c = conn.cursor()
	q = "SELECT * FROM " + name + " " + " WHERE " + column + " = '" + elem + "'"
	#print q
	c.execute(q)
	r = c.fetchone()
	c.close()
	conn.close()
	return r

def readRepeatFile(filename):
	fin = open(filename,'r')
	out_dict = {}
	flines = fin.readlines()
	it = iter(flines)
	tuple_list = zip(it,it)
	for entry in tuple_list:
		accessions =  entry[0][1:].replace('\n', '').split('|')
		for acc in accessions:
			out_dict[acc] = entry[1].replace('\n', '')
	return out_dict


def match_repeat_to_spacer(repeat_data, spacer_file_name, dbName):
	spacer_file = open(spacer_file_name,'r').readlines()
	it = iter(spacer_file)
	tuple_list = zip(it,it)
	for entry in tuple_list:
		accessions =  entry[0][1:].replace('\n', '').split('|')
		sequence = entry[1].replace('\n', '')
		#print entry[1]
		for acc in accessions:
			acc_match = "_".join(acc.split("_")[:-1])
			spacer_id =  SQL_search(sequence, 'Spacer', 'SpacerSequence', dbName)[0]
			repeat_id =  SQL_search(repeat_data[acc_match], 'Repeat', 'RepeatSequence', dbName)[0]
			try:
				match = repeat_data[acc_match]
				acc
			except KeyError:
				print "Error: Wrong Accession code"	
			SQL_add(str((spacer_id,repeat_id)), 'SpacerRepeatPair', dbName, '(SpacerID, RepeatID)')		


spacerFile = 'data/Spacerdatabase.txt'
repeatFile = 'data/DRdatabase.txt'
dbName = 'crispr.sqlite'
repeat_dict = readRepeatFile(repeatFile)
match_repeat_to_spacer(repeat_dict, spacerFile, dbName)
