import requests
from pattern import web
import re
import csv

def get_dom(url):
    html = requests.get(url).text
    dom = web.Element(html)
    return dom

def get_taxons_from_CRISPRdb():
    url = "http://crispr.u-psud.fr/crispr/"
    dom_homepage = get_dom(url)
    container = dom_homepage('div[class="strainlist"]')[0]
    crisprs = []
    for link in container('a'):
        crispr = {}
        crispr['Name'] = link.content.encode('ascii','ignore')
        crispr['Taxon_id'] = link.href.encode('ascii','ignore')[46:]
        crisprs.append(crispr)
    return crisprs

def get_genome_properties(crispr):
    url = "http://crispr.u-psud.fr/cgi-bin/crispr/SpecieProperties.cgi?Taxon_id=" + crispr['Taxon_id']
    dom = get_dom(url)
    table = dom('table[class="primary_table"]')[1]
    crispr['Ref_seqs']={}
    for sequence in table('tr'):
        if sequence('th')==[]:
            crispr['Ref_seqs'][sequence('td')[1].content.encode('ascii','ignore')]=[]
    return crispr

def get_CRISPR_properties(crispr):
    for ref_seq in crispr['Ref_seqs'].keys():
        url = "http://crispr.u-psud.fr/crispr/CRISPRProperties.php?RefSeq=" + ref_seq + "&Taxon=" + crispr['Taxon_id']
        dom = get_dom(url)
        div = dom('div[class="rightcontent"]')[1]
        table = div('table[class="primary_table"]')[0]
        tbody = table('tbody')[0]
        for crispr_id in tbody('tr'):
            cell = crispr_id('td')[1]
            crispr_id_value = cell('font')[0].content.encode('ascii','ignore').replace('<br/n/>','')
            crispr_id_value = crispr_id_value.replace('\n','')
            crispr_id_value = crispr_id_value.replace('\t','')
            crispr_id_value = crispr_id_value.replace('<br />','')
            crispr['Ref_seqs'][ref_seq].append(crispr_id_value)
    return crispr

def get_begin(source):
    begin = re.search('(?<=Crispr_begin_position: )\d+', source)
    begin = begin.group(0)
    return begin

def get_end(source):
    end = re.search('(?<=Crispr_end_position: )\d+', source)
    end = end.group(0)
    return end

def get_loc_id(source):
    loc_id = re.search('(?<=Crispr Rank in the sequence: )\d+', source)
    loc_id = loc_id.group(0)
    return loc_id

def download(source,a,b):
    directory = "files/"
    namefile =  a + "_" + b
    if not os.path.exists(directory):
        os.makedirs(directory)
    f = open(directory + namefile, 'w')
    f.write(source)
    f.close()

def get_results():
    crisprs =  get_taxons_from_CRISPRdb()
    n = len(crisprs)
    i = 0.
    print "Getting genome and CRISPR properties..."
    for crispr in crisprs:
        print str(i/n*100) + "%"
        crispr = get_genome_properties(crispr)
        crispr = get_CRISPR_properties(crispr)
        i += 1
    results = []
    i = 0.
    print "Downloading files..."
    for crispr in crisprs:
        print str(i/n*100) + "%"
        for ref_seq in crispr['Ref_seqs'].keys():
            for crispr_id in crispr['Ref_seqs'][ref_seq]:
                params = {'checked[]': crispr_id, 'Taxon': crispr['Taxon_id']}
                r = requests.post("http://crispr.u-psud.fr/crispr/crispr.php", data=params)
                file_crispr_seq = web.Element(r.text)
                table = file_crispr_seq('table[class="crisprs_table"]')[0]
                url = "http://crispr.u-psud.fr" + table('form')[-2].action
                source = requests.get(url).text
                begin = get_begin(source)
                end = get_end(source)
                main_accession_number = ref_seq
                loc_id = get_loc_id(source)
                download(source,main_accession_number,loc_id)
                results.append([main_accession_number,loc_id,begin,end])
                i += 1

    with open('results.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        for result in results:
            writer.writerow(result)
            
get_results()
