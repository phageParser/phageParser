#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 12:44:58 2017

@author: madeleine
"""

'''
Issue #160

Author: Madeleine
Usage:python3 populate_casgenes.py

Populates CasProtein and OrganismCasPair tables

'''

import os
import numpy as np

casprofiles = "../data/crispr_type.csv" # table of cas proteins from Makarova et al 2015

with open(casprofiles, 'r', encoding='utf8') as f:
    profiles = np.loadtxt(f, dtype='S', delimiter='\t')

def make_profile_dict(profiles):
    """ 
    make dictionary of profile accessions to gene name and crispr types
    """
    profile_dict = {}
    for row in profiles:
        profile_dict[row[0].decode('utf8')] = [row[3].decode('utf8'),row[5].decode('utf8')]
    return profile_dict  
    
def populate_cas():
    print("Populating CasProtein table...")
    for row in profiles:
        profileID = row[0].decode('utf8')
        function = row[2].decode('utf8')
        gene = row[3].decode('utf8')
        group = row[4].decode('utf8')
        type_spec = row[5].decode('utf8')
        casprotein, created = CasProtein.objects.get_or_create(profileID = profileID, function = function, gene = gene, group = group, type_specificity = type_spec)
    print("Done.")
    return
    
def populate_organismcaspair():
    print("Populating OrganismCasPair table...")
    for fn in os.listdir("gbfiles/hmmeroutput"): # these are files split by organism accession
        data = np.loadtxt("gbfiles/hmmeroutput/%s" %fn,dtype='S') # HMMER output results
        accid = fn.rsplit('.')[0] # organism accession numberr

        organismset = Organism.objects.filter(accession=accid) # retrieve organism from database
        
        if not organismset:
            print('Organism with accid %s not found in db' % accid)
        organism = organismset[0] # this is the organism FK for the field organism
        
        querylist = []
        for row in data: # iterate over HMMER matches to cas protein profiles
            query = row[2].decode('utf8') # cds start and end
            evalue = row[4].decode('utf8')
            target_match = row[0].decode('utf8')
            casproteinset = CasProtein.objects.filter(profileID=target_match) # retrieve cas protein entry from database
            
            if not casproteinset:
                print('Cas protein with profileID %s not found in db' % target_match)
            casprotein = casproteinset[0]  # this is the cas protein FK for the field casprotein

            if query not in querylist: # we only keep 1 match if there are multiple
                querylist.append(query)
                
                try: # check if sequence is complemented - if yes, start will be > end 
                    query.index('complement(')
                    query = query[11:-1]
                except:
                    pass
                
                start, dot, end = query.rsplit('.')
                start = int(start)
                end = int(end)
                
                evalue = float(evalue)
                                  
                osrpair, created = OrganismCasPair.objects.get_or_create(organism=organism, 
                                                                         casprotein=casprotein, 
                                                                         genomic_start=start, 
                                                                         genomic_end=end, 
                                                                         evalue=evalue)                                                                    
    print("Done.")
    return

    
if __name__ == '__main__':
    import django
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "phageAPI.settings")
    django.setup()
    from restapi.models import Organism, CasProtein, OrganismCasPair
    populate_cas() 
    populate_organismcaspair()
    