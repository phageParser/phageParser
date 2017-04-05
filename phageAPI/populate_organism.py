#!/usr/bin/env python

""" Populate Organism model by using a file with accession numbers.
Fetch GenBank names for these accession numbers.
Fill in missing names from GenBank information and create models.

Adapted from organism_name_update.py

USAGE:
python populate_organism.py <email>

where:
<email> is the email address associated with your NCBI account

DEPENDENCIES:
Biopython
"""
import os
import sys
from Bio import Entrez
from Bio.SeqIO import parse
def populate_organism(org_accession_path='../data/bac_accession_list.txt'):
    def add_organism(name, accession):
        o, created= Organism.objects.get_or_create(name=name, accession=accession) #gets the object, this also checks for duplicates
        return o
    def fetch_accessions():
        lines = [line.rstrip('\n') for line in open(org_accession_path)]
        return lines
    def merge_acc_names(accession_list):
        acc_name_dict = {}
        #define email for entrez login
        db           = "nuccore"
        Entrez.email = sys.argv[1] #need to provide email for NCBI
        # Doing 100 by 100 to make sure requests to NCBI are not too big
        for i in range(0, len(accession_list), 100):
            j = i + 100
            if (j >= len(accession_list)):
                j = len(accession_list)

            sys.stderr.write("Fetching entries from %s to %s from GenBank\n" % (i, j))
            sys.stderr.flush()
            result_handle = Entrez.efetch(db=db, rettype="gb", id=accession_list[i:j])

            # Populate result per organism name
            for record in parse(result_handle, 'genbank'):
                # Using NCBI name, which should match accession number passed
                acc_name_dict[record.name] = record.annotations['organism']

        return acc_name_dict
    accession_list = fetch_accessions()
    acc_name_dict = merge_acc_names(accession_list)
    for acc in acc_name_dict.keys():
        add_organism(name=acc_name_dict[acc], accession=acc)
if __name__ == '__main__':
    print "Starting orgnanism population script"
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism
    populate_organism()
