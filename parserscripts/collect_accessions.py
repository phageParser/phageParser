# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 11:26:23 2017

@author: madeleine
"""

import os
import numpy as np

accessions = []

with open("genomes_euks.txt", "rb") as f:
    header = f.readline()
    data = f.readlines()

for line in data:
    line = line.rsplit('\t')
    accessionlist = line[9]
    try:
        count = accessionlist.count(':')
    except:
        continue
    for i in range(count):
        i = accessionlist.index(':')
        j = accessionlist[i:].index('.')
        acc = accessionlist[i+1:i+j]
        accessions.append(acc)
        accessionlist = accessionlist[i+j+1:]
    
accessions = np.array(accessions,dtype='object')
np.savetxt("eukaryote_accessions.txt",accessions,fmt='%5s')
