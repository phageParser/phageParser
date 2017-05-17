import os
import fetch
from IPython import embed
from Bio import SeqIO
from lxml import html, etree
import requests
import pandas
from prunedict import prunedict


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


if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism, Spacer, Repeat, SpacerRepeatPair, OrganismSpacerRepeatPair
    print "Populating spacerrepeatpair and organismspacerrepeatpair tables"
    datapath = "../data"
    sfile, rfile = get_spacerrepeatfiles(datapath)
    gendict = prunedict(
        addpositionstodict(
            addspacerstodict(
                repeatfiletodict(rfile), sfile)))
    import pickle
    with open('gendict.pickle', 'wb') as f:
        pickle.dump(gendict, f, protocol=pickle.HIGHEST_PROTOCOL)
    print('Created dictionary and dumped data to gendict.pickle')
    embed()
