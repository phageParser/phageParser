# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 15:18:42 2016

@author: Ahmed

Given a 
echo <file> | python phageblast.py <genome> > <output>
where:
<file> is the name of a fasta file containing protein sequences of the antiCRISPR genes of interest
<genome> is the name of a file containing the bacterial genome of interest 
<output> is the name of the (xml) file you'd like to write the results to  

http://www.ncbi.nlm.nih.gov/books/NBK279675/ 
^ see Table C5 for other tblastn options 

DEPENDENCIES:
Biopython
blast+
"""

from Bio.Blast.Applications import NcbitblastnCommandline 
import sys
import subprocess 

genome = sys.argv[-1]

#loads input file
seq = [l.strip() for l in sys.stdin if l.strip()]
seq = ''.join(seq)

#builds command line string 
tblastn_obj = NcbitblastnCommandline(cmd = 'tblastn', query=seq, subject=genome.read(), 
                                     evalue=10, outfmt = 5)

tblastn_cline = str(tblastn_obj) 

#executes command line string in bash 
subprocess.call(tblastn_cline, shell=True) 
