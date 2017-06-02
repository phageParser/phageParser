#!/usr/bin/env python

import os
import sys
import argparse
import requests
import pandas
import pickle
from lxml import html, etree
from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
import django
django.setup()

from util.acc import read_accession_file
from util.prunedict import prunedict
from util import fetch
from restapi.models import (
    Organism,
    Spacer,
    Repeat,
    OrganismSpacerRepeatPair,
    AntiCRISPR
)

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))


def populate_organism():
    def add_organism(name, accession):
        # get the object, this also checks for duplicates
        o, created = Organism.objects.get_or_create(name=name, accession=accession)
        return o

    def merge_acc_names(accession_list):
        acc_name_dict = {}
        db = "nuccore"
        # Doing batches of 200 to make sure requests to NCBI are not too big
        for i in range(0, len(accession_list), 200):
            j = i + 200

            result_handle = Entrez.efetch(db=db, rettype="gb", id=accession_list[i:j])

            # Populate result per organism name
            records = SeqIO.parse(result_handle, 'genbank')
            for record in tqdm(records):
                # Using NCBI name, which should match accession number passed
                acc_name_dict[record.name] = record.annotations['organism']
        return acc_name_dict

    with open(os.path.join(DATA_DIR, 'bac_accession_list.txt')) as f:
        accession_list = list(read_accession_file(f))
    acc_name_dict = merge_acc_names(accession_list)
    for acc in acc_name_dict:
        add_organism(name=acc_name_dict[acc], accession=acc)


def get_spacerrepeatfiles():
    spath = os.path.join(DATA_DIR, "spacerdatabase.txt")
    surl = 'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/Spacer/Spacerdatabase'
    rpath = os.path.join(DATA_DIR, "repeatdatabase.txt")
    rurl = 'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/DR/DRdatabase'
    fetch.fetch(spath, surl)
    fetch.fetch(rpath, rurl)
    return spath, rpath


def repeatfiletodict(rfile):
    rdict = {}
    repeatrecords = SeqIO.parse(rfile, 'fasta')
    for record in repeatrecords:
        accessions = record.name.split('|')
        sequence = str(record.seq)
        for acc in accessions:
            rdict[acc] = {'RepeatSeq': sequence}
    return rdict


def addspacerstodict(gendict, sfile):
    spacerrecords = SeqIO.parse(sfile, 'fasta')
    for record in spacerrecords:
        accessions = record.name.split('|')
        sequence = str(record.seq)
        for acc in accessions:
            acc_elems = acc.split('_')
            order = acc_elems[-1]
            acc_id = '_'.join(acc_elems[:-1])
            try:
                if 'Spacers' in gendict[acc_id]:
                    gendict[acc_id]['Spacers'][order] = sequence
                else:
                    gendict[acc_id]['Spacers'] = {order: sequence}
            except KeyError:
                print('Error on accession id:  %s' % acc_id)
    return gendict


def addpositionstodict(gendict):
    print("Downloading position information from web...")
    for accidwithloc in tqdm(gendict):
        if 'Start' in gendict[accidwithloc]:
            continue
        accid = '_'.join(accidwithloc.split('_')[:-1])
        url = 'http://crispr.i2bc.paris-saclay.fr/crispr/crispr_db.php?checked%5B%5D={}'.format(
            accid)
        page = requests.get(url)
        htmltable = html.fromstring(page.content).xpath(
            "//table[normalize-space(@class)='primary_table']")[1]
        strtable = etree.tostring(htmltable)
        # converts to pandas df and then to numpy array then drop titles
        arrtable = pandas.read_html(strtable)[0].as_matrix()[2:]
        for row in arrtable:
            if row[0] in gendict:
                gendict[row[0]]['Start'] = row[2]
                gendict[row[0]]['Stop'] = row[3]
            else:
                if row[1] != 'questionable':
                    print("Can't find %s in local files" % row[0])
    return gendict


def populate_fromlocus(locid, locus):
    accid = '_'.join(locid.split('_')[:-1])
    organismset = Organism.objects.filter(accession=accid)
    if not organismset:
        print('Organism with accid %s not found in db' % accid)
        return
    organism = organismset[0]
    repeat, _ = Repeat.objects.get_or_create(sequence=locus['RepeatSeq'])
    posindex = int(locus['Start'])
    spacers = locus['Spacers']
    for order in sorted(spacers):
        spacer, _ = Spacer.objects.get_or_create(sequence=spacers[order])
        pairstart = posindex
        pairend = pairstart + len(spacer.sequence) + len(repeat.sequence)
        posindex = pairend
        osrpair, _ = OrganismSpacerRepeatPair.objects.get_or_create(organism=organism,
                                                                    spacer=spacer,
                                                                    repeat=repeat,
                                                                    order=int(
                                                                        order),
                                                                    genomic_start=int(
                                                                        pairstart),
                                                                    genomic_end=int(pairend))
        spacer.save()
        osrpair.save()
    repeat.save()
    organism.save()

def populate_osrpair():
    print('Downloading files and gathering online data.')
    sfile, rfile = get_spacerrepeatfiles()
    gendict = prunedict(
        addpositionstodict(
            addspacerstodict(
                repeatfiletodict(rfile), sfile)))
    with open('genedict.pickle', 'wb') as f:
        pickle.dump(gendict, f, protocol=pickle.HIGHEST_PROTOCOL)

    print('Created dictionary and dumped data to genedict.pickle')
    print("Populating Spacer, Repeat, SpacerRepeatPair, OrganismSpacerRepeatPair tables")
    for locid in tqdm(gendict):
        populate_fromlocus(locid, gendict[locid])


def populate_anticrispr():
    with open(os.path.join(DATA_DIR, 'antiCRISPR_accessions.txt')) as f:
        accession_list = list(read_accession_file(f))
    print("Fetching AntiCRISPR entries")
    result_handle = Entrez.efetch(db='protein', rettype="fasta", id=accession_list)
    for record in tqdm(SeqIO.parse(result_handle, 'fasta')):
        spacer, _ = AntiCRISPR.objects.get_or_create(
            accession=record.name,
            sequence=str(record.seq))
        spacer.save()


def main():
    parser = argparse.ArgumentParser(
        description='Populate the phageParser database with data from NCBI')
    parser.add_argument('email', nargs=1,
                        help='your email address (does not need to be registered, just used to identify you)')
    args = parser.parse_args()

    Entrez.email = args.email

    print("Starting organism population")
    populate_organism()
    print("Starting OSR population")
    populate_osrpair()
    print("Starting AntiCRISPR population")
    populate_anticrispr()


if __name__ == '__main__':
    main()
