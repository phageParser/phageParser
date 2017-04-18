#!/usr/bin/env python
"""Fetch GenBank entries for given accessions. 

Adapted from https://www.biostars.org/p/66921/

USAGE:
cat <file> | python acc2gb.py <email> <db> <rettype> > <output>

where:
<file> is the name of a file containing accession numbers to download
<email> is the email address associated with your NCBI account
<db> is the NCBI database ID
<rettype> is the retrieval type
<output> is the name of the file you'd like to write the results to

EXAMPLES:
Case 1: rettype = gbwithparts, db = nuccore - downloads genbank file with metadata and fasta DNA sequence (i.e. for downloading bacterial genomes with metadata)
Case 2: rettype = fasta, db = protein - downloads fasta file with protein sequence (i.e. for downloading antiCRISPR protein sequences for BLAST)

DEPENDENCIES:
Biopython

"""

import sys
from Bio import Entrez

#define email for entrez login
db           = sys.argv[2]
Entrez.email = sys.argv[1]

#get accession numbers out of stdin
accs = [ l.strip() for l in sys.stdin if l.strip() ]

#fetch
sys.stderr.write( "Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
for i,acc in enumerate(accs):
  try:
    sys.stderr.write( " %9i %s          \r" % (i+1,acc))  
    handle = Entrez.efetch(db=db, rettype=sys.argv[3], retmode="text", id=acc)
    #print output to stdout
    sys.stdout.write(handle.read())
  except:
    sys.stderr.write( "Error! Cannot fetch: %s        \n" % acc) 
