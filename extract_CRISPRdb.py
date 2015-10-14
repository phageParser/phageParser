import requests
from pattern import web
import re
import csv
import os

def get_dom(url):
    html = requests.get(url).text
    dom = web.Element(html)
    return dom

def get_taxons_from_CRISPRdb():
    taxon_ids = []
    if os.path.exists('taxon_ids.csv'):
        with open('taxon_ids.csv', 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                taxon_ids.append(row[0])
            return taxon_ids
    url = "http://crispr.u-psud.fr/crispr/"
    dom_homepage = get_dom(url)
    container = dom_homepage('div[class="strainlist"]')[0]
    for link in container('a'):
        taxon_id = link.href.encode('ascii','ignore')[46:]
        taxon_ids.append(taxon_id)
        with open('taxon_ids.csv', 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([taxon_id])
    return taxon_ids

def get_sequences(taxon_id):
    url = "http://crispr.u-psud.fr/cgi-bin/crispr/SpecieProperties.cgi?Taxon_id=" + taxon_id
    dom = get_dom(url)
    table = dom('table[class="primary_table"]')[1]
    sequences = []
    for sequence in table('tr'):
        seq = []
        if sequence('th')==[]:
            seq = {}
            seq['Taxon_id'] = taxon_id
            seq['RefSeq'] = sequence('td')[1].content.encode('ascii','ignore')
            sequences.append(seq)
    return sequences

def get_loc_ids(sequence):
    url = "http://crispr.u-psud.fr/crispr/CRISPRProperties.php?RefSeq=" + sequence['RefSeq']
    dom = get_dom(url)
    div = dom('div[class="rightcontent"]')[1]
    table = div('table[class="primary_table"]')[0]
    tbody = table('tbody')[0]
    loc_ids = []
    for crispr_id in tbody('tr'):
        cell = crispr_id('td')[1]
        crispr_id_value = cell('font')[0].content.encode('ascii','ignore').replace('<br/n/>','')
        crispr_id_value = crispr_id_value.replace('\n','')
        crispr_id_value = crispr_id_value.replace('\t','')
        crispr_id_value = crispr_id_value.replace('<br />','')
        loc_ids.append(crispr_id_value.split('_')[-1])
    return loc_ids

def get_positions(sequence,loc_id):
    params = {}
    params['Taxon'] =  sequence['Taxon_id']
    crispr_id = sequence['RefSeq'] + "_" + loc_id
    params['checked[]'] = crispr_id
    r = requests.post("http://crispr.u-psud.fr/crispr/crispr.php", data=params)
    source = web.Element(r.text)
    print "Crispr id: " + crispr_id
    table = source('table[class="crispr"]')[0]
    tr = table('tr')[3]
    td = tr('td')
    begin = re.search('(?<=</span>) +\d+', td[0].content)
    begin = begin.group(0).split(' ')[-1]
    end = re.search('(?<=</span>) +\d+', td[1].content)
    end = end.group(0).split(' ')[-1]
    return [begin,end]

def get_results():
    print "Getting taxon ids..."
    taxon_ids =  get_taxons_from_CRISPRdb()
    n = len(taxon_ids)
    i = 0.
    print "Getting positions..."
    restart = False
    var = False
    if os.path.exists('results.csv'):
        f = open('last','r')
        last = f.read()
        [last_taxon_id,last_refseq,last_loc_id] = last.split(',')
        restart = True
    for taxon_id in taxon_ids:
        print "{:.2%}".format(i/n)
        print "Taxon id: " + taxon_id
        if(restart == True and last_taxon_id == taxon_id):
            var = True
            continue
        if (var == True or restart == False):
            sequences = get_sequences(taxon_id)
            for sequence in sequences:
                loc_ids = get_loc_ids(sequence)
                sequence['Loc_ids'] = loc_ids
                for loc_id in sequence['Loc_ids']:
                    [begin,end] = get_positions(sequence,loc_id)
                    result = [sequence['RefSeq'],loc_id,begin,end]
                    with open('results.csv', 'a') as csvfile:
                        writer = csv.writer(csvfile, delimiter=',')
                        writer.writerow(result)
                    f = open('last','w')
                    f.write(sequence['Taxon_id'] + "," + sequence['RefSeq'] + "," + loc_id)
        i += 1

get_results()
