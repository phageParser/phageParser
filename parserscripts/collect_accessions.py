# -*- coding: utf-8 -*-
'''
Created on Tue Aug 15 11:26:23 2017
This script processes NCBI data exports (from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/)
to extract lists of organism accession IDs, later used for downloading complete genome sequences.

@author: mbonsma / thisisjaid
'''

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(add_help=True, description='''phageParser - collect_accessions.py - 
                                 This script processes NCBI data exports 
                                 (from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/)
                                 to generate lists of organism accession IDs, later used for downloading 
                                 complete genome sequences''')
parser.add_argument('-o', metavar='output_file', action='store', type=str, default='accessions.csv', 
                    help='Full path to results output file (default: accessions.csv)')
parser.add_argument('file', action='store', help='Full path to NCBI data export file')

args = parser.parse_args()

records = pd.read_csv(args.file,sep='\t', header=None)
records.to_csv(args.o,header=False,sep=',',columns=[1],index=False)
print('Accession file successfully written out at:',args.o)
