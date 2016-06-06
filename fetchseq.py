# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 14:15:58 2016

@author: Ahmed
"""

"""Fetch sequences in fasta format for given antiCRISPR protein accession numbers. 
Adapted from https://www.biostars.org/p/66921/ + acc2gb.py
USAGE:
cat <file> | python fetchseq.py <email> > <output>
where:
<file> is the name of a file containing accession numbers to download sequences for
<email> is the email address associated with your NCBI account
<output> is the name of the file you'd like to write the results to
DEPENDENCIES:
Biopython
"""

import sys
from Bio import Entrez

#define email for entrez login
db           = "protein"
Entrez.email = sys.argv[-1]

#get accession numbers out of stdin
accs = [ l.strip() for l in sys.stdin if l.strip() ]

#fetch
sys.stderr.write( "Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
for i,acc in enumerate(accs):
  try:
    sys.stderr.write( " %9i %s          \r" % (i+1,acc))  
    handle = Entrez.efetch(db=db, rettype="fasta", retmode="text", id=acc)
    #print output to stdout 
    sys.stdout.write(handle.read())
    with open('sequences.txt','a') as ourFile:
        ourFile.write(handle.read())
  except:
    sys.stderr.write( "Error! Cannot fetch: %s        \n" % acc) 
   
