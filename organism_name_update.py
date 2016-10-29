#!/usr/bin/env python

# TODO: Connect this to database!

"""Fetch GenBank names for given accessions.

Adapted from bac_info.py and acc2gb.py

USAGE:
cat <file> | python acc2gb.py <email> > <output>

where:
<file> is the name of a file containing accession numbers to download
<email> is the email address associated with your NCBI account
<output> is the name of the file you'd like to write the results to

DEPENDENCIES:
Biopython
"""

import sys
from Bio import Entrez
from Bio.SeqIO import parse


#define email for entrez login
db           = "nuccore"
Entrez.email = sys.argv[1]

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


#get accession numbers out of stdin
accession_numbers = [ l.strip() for l in sys.stdin if l.strip() ]

print(fetch_names(accession_numbers))
