#!/usr/bin/python
import sqlite3

# THIS file will parse in two files named 'DRdatabase' and 'Spacerdatabase', 
# and populate the 'crispr.sqlite' file tables Spacers and Repeats, as per 
# ISSUE #99 https://github.com/goyalsid/phageParser/issues/99

def readInputFile(filename):
	fin = open(filename,'r')
	flines = fin.readlines()
	fileDict = dict()
	for f in flines:
		f = f.replace('\n', '')
		if(f[0] == '>'):
			accessions = f[1:].split('|')
		else:
			fileDict[f] = accessions
	return fileDict

def printDict(d):
	for i in d:
		print(i, d[i])

def getLargestID(c, name, idName):
	# NOT USED
	# SELECT * FROM Spacer ORDER BY SpacerID DESC LIMIT 1
	q = "SELECT * FROM " + name + " ORDER BY " + idName + " DESC LIMIT 1"
	c.execute(q)
	return c.fetchone()[0]

def SQL_add(elem, name, dbName):
	conn = sqlite3.connect(dbName)
	c = conn.cursor()
	colName = name+'Sequence'
	q = "INSERT INTO " + str(name) + " ("  + str(colName) + ") VALUES ('"+ elem+"')"
	try:
		c.execute(q)
		conn.commit()
	except sqlite3.IntegrityError:
		print("Already In Database", elem)
	c.close()
	conn.close()

def SQL_search(elem, name, column, dbName):
	conn = sqlite3.connect(dbName)
	c = conn.cursor()
	q = "SELECT * FROM " + name + " " + " WHERE " + column + " = '" + elem + "'"
	# print q
	c.execute(q)
	r = c.fetchone()
	c.close()
	conn.close()
	return r

def populateDB(dbName, tableName, fileName):
	inputDict = readInputFile(fileName)
	for d in inputDict:
		SQL_add(d, tableName, dbName)


repeatFile = 'data/DRdatabase.txt'
spacerFile = 'data/spacerdatabase.txt'
dbName = 'crispr.sqlite'

populateDB(dbName, 'Repeat', repeatFile)
populateDB(dbName, 'Spacer', spacerFile)
