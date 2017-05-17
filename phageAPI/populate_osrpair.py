import os
import fetch
from IPython import embed
from Bio import SeqIO
from lxml import html, etree
import requests
import pandas
import pickle
from prunedict import prunedict
from tqdm import tqdm

def get_spacerrepeatfiles(datapath):
    spath = os.path.join(datapath, "spacerdatabase.txt")
    surl = 'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/Spacer/Spacerdatabase'
    rpath = os.path.join(datapath, "repeatdatabase.txt")
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
    accidswithloc = gendict.keys()
    for accidwithloc in accidswithloc:
        if 'Start' in gendict[accidwithloc]:
            continue
        accid = '_'.join(accidwithloc.split('_')[:-1])
        url = 'http://crispr.i2bc.paris-saclay.fr/crispr/crispr_db.php?checked%5B%5D={}'.format(
            accid)
        page = requests.get(url)
        htmltable = html.fromstring(page.content).xpath(
            "//table[normalize-space(@class)='primary_table']")[1]
        strtable = etree.tostring(htmltable)
        # embed()
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
    for order in spacers:
        spacer, _ = Spacer.objects.get_or_create(sequence=spacers[order])
        spacerrepeatpair, _ = SpacerRepeatPair.objects.get_or_create(spacer=spacer, repeat=repeat)
        pairstart = posindex
        pairend = pairstart + len(spacer.sequence) + len(repeat.sequence)
        posindex = pairend
        osrpair, _ = OrganismSpacerRepeatPair.objects.get_or_create(organism=organism,
                                                                    spacerrepeatpair=spacerrepeatpair,
                                                                    order=int(order),
                                                                    genomic_start=int(pairstart),
                                                                    genomic_end=int(pairend))
        spacer.save()
        spacerrepeatpair.save()
        osrpair.save()
    repeat.save()
    organism.save()
if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism, Spacer, Repeat, SpacerRepeatPair, OrganismSpacerRepeatPair
    datapath = "../data"
    print('Downloading files and gathering online data.')
    sfile, rfile = get_spacerrepeatfiles(datapath)
    gendict = prunedict(
        addpositionstodict(
            addspacerstodict(
                repeatfiletodict(rfile), sfile)))
    with open('genedict.pickle', 'wb') as f:
        pickle.dump(gendict, f, protocol=pickle.HIGHEST_PROTOCOL)
    print('Created dictionary and dumped data to gendict.pickle')
    print("Populating Spacer, Repeat, SpacerRepeatPair, OrganismSpacerRepeatPair tables")
    for locid in tqdm(gendict):
        populate_fromlocus(locid, gendict[locid])
    embed()
