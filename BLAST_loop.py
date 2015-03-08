# This script BLASTs all the files in /spacers against a downloaded copy of the phage database.

from Bio.Blast.Applications import NcbiblastnCommandline
#help(NcbiblastnCommandline)
import os

#change this directory as needed. Note that local blast may experience issues if you're not working in the blast directory (here NCBI/blast-2.2.30+)
directory = "C:/Users/madeleine/Documents/NCBI/phageParser/data"

for fn in os.listdir("%s/spacers" %directory):
    
    query1 = "%s/spacers/%s" %(directory,fn) 
    ext = fn.index('.')
    outfile1 = fn.replace(fn[ext+1:],("txt"))
    outfile = "%s/phages/%s" %(directory,outfile1)

    # These parameters are more or less the same as the ones on PhagesDB.org
    blastn_obj = NcbiblastnCommandline(query=query1, db="phagedb", evalue=10, num_descriptions = 100, num_alignments = 100, dust = "no", task = "blastn",reward = 1, penalty = -3, out = outfile)

    stdout, stderr = blastn_obj()
