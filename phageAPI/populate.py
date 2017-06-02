#!/usr/bin/env python

import os
import sys
import argparse
from Bio import Entrez
from Bio.SeqIO import parse

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
import django
django.setup()

from util.acc import read_accession_file
from restapi.models import Organism

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))


def populate_organism():
    def add_organism(name, accession):
        # get the object, this also checks for duplicates
        o, created = Organism.objects.get_or_create(name=name, accession=accession)
        return o

    def merge_acc_names(accession_list):
        acc_name_dict = {}
        db = "nuccore"
        # Doing 100 by 100 to make sure requests to NCBI are not too big
        for i in range(0, len(accession_list), 200):
            j = i + 200
            if (j >= len(accession_list)):
                j = len(accession_list)

            print("Fetching organism entries from %s to %s from GenBank\n" % (i, j))
            result_handle = Entrez.efetch(db=db, rettype="gb", id=accession_list[i:j])

            # Populate result per organism name
            records = parse(result_handle, 'genbank')
            count = 0
            for record in records:
                # Using NCBI name, which should match accession number passed
                acc_name_dict[record.name] = record.annotations['organism']
                count += 1
            print(count)
        return acc_name_dict

    with open(os.path.join(DATA_DIR, 'bac_accession_list.txt')) as f:
        accession_list = list(read_accession_file(f))
    acc_name_dict = merge_acc_names(accession_list)
    for acc in acc_name_dict.keys():
        add_organism(name=acc_name_dict[acc], accession=acc)


def main():
    parser = argparse.ArgumentParser(
        description='Populate the phageParser database with data from NCBI')
    parser.add_argument('email', nargs=1,
                        help='your email address (does not need to be registered, just used to identify you)')
    args = parser.parse_args()

    Entrez.email = args.email

    print("Starting organism population")
    populate_organism()


if __name__ == '__main__':
    main()
