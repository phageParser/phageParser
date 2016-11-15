#!/usr/bin/env python

""" Determine accession numbers in Organism table with missing names. 
Fetch GenBank information for these accession numbers.
Fill in missing names from GenBank information.

Adapted from bac_info.py and acc2gb.py

USAGE:
python organism_name_update.py <email> database.sqlite

where:
<email> is the email address associated with your NCBI account

DEPENDENCIES:
Biopython
"""

import sys
from Bio import Entrez
from Bio.SeqIO import parse
import sqlite3

sqlite_file = sys.argv[2]


conn = sqlite3.connect(sqlite_file)
conn.row_factory = lambda cursor, row: row[0]
c = conn.cursor()


#define email for entrez login
db           = "nuccore"
Entrez.email = sys.argv[1]



def missing_names(sqlite_file):
	'''
	Check for and return accession numbers in Organism table that don't have names
	'''
	c.execute("SELECT  Accession FROM Organism WHERE OrganismName IS NULL")
	id_list = c.fetchall()
	return id_list


def fetch_names(id_list):
    organism_names = {}

    # Doing 100 by 100 to make sure requests to NCBI are not too big
    for i in range(0, len(id_list), 100):
        j = i + 100
        if (j >= len(id_list)):
            j = len(id_list)

        sys.stderr.write("Fetching entries from %s to %s from GenBank\n" % (i, j))
        sys.stderr.flush()
        result_handle = Entrez.efetch(db=db, rettype="gb", id=id_list[i:j])

        # Populate result per organism name
        for record in parse(result_handle, 'genbank'):
            # Using NCBI name, which should match accession number passed
            organism_names[record.name] = record.annotations['organism']

    return organism_names

def insert_names(organism_names):
	for key in organism_names.keys():
		c.execute("UPDATE Organism SET OrganismName=({name}) WHERE Accession=({acc})".\
    	    format(name=organism_names[key], acc=key))	


id_list = missing_names(sqlite_file)
organism_names = fetch_names(id_list)
insert_names(organism_names)


conn.commit()
conn.close()


