#!/usr/bin/python
import sqlite3

def printDict(d):
	for i in d:
		print i,d[i]

def SQL_add(elem, name, columns, dbName):
	conn = sqlite3.connect(dbName)
	c = conn.cursor()
	if(type(elem) is list):
		# It is assumed that the length of list and of columns will always match
		q = "INSERT INTO " + str(name) + " (" + ",".join(columns) + ") VALUES ('"+ "','".join(elem) + "')"
	else:
		q = "INSERT INTO " + str(name) + columns + " VALUES " + elem
	# print q
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
	if(type(elem) is list):
		# It is assumed that the length of list and of columns will always match
		q = "SELECT * FROM " + name + " " + " WHERE "
		for i in range(0, len(elem)):
			q += str(column[i]) + " = " + str(elem[i]) + " and "
		q = q[:-5]
	else:
		q = "SELECT * FROM " + name + " " + " WHERE " + column + " = " + elem
	# print q
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

def getIdAndAdd(elem, name, column, dbName):
	item_id =  SQL_search(elem, name, column, dbName)
	if(item_id is None):
		SQL_add(elem, name, column, dbName)
		item_id =  SQL_search(elem, name, column, dbName)
	return item_id[0]


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
			spacer_id =  getIdAndAdd("('"+sequence+"')", 'Spacer', '(SpacerSequence)', dbName)
			repeat_id =  getIdAndAdd("('"+repeat_data[acc_match]+"')", 'Repeat', '(RepeatSequence)', dbName)
			try:
				match = repeat_data[acc_match]
				acc
			except KeyError:
				print "Error: Wrong Accession code"	
			getIdAndAdd([str(spacer_id),str(repeat_id)], 'SpacerRepeatPair', ['SpacerID', 'RepeatID'], dbName)


spacerFile = 'data/spacerdatabase.txt'
repeatFile = 'data/DRdatabase.txt'
dbName = 'crispr.sqlite'
repeat_dict = readRepeatFile(repeatFile)
match_repeat_to_spacer(repeat_dict, spacerFile, dbName)
