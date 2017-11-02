# parse_Genbank.py

from Bio.SeqIO import parse

genbank_file = 'data/Genbank_example.txt'

outdict = {}

result_handle = open(genbank_file)

records = parse(result_handle, 'genbank')

for x in records:
    try:
        acc = x.name
    except ValueError:
        continue
    try:
        org = x.annotations['organism']
    except:
        print("no organism name")
    try:
        tax = x.annotations['taxonomy']
    except:
        print("no taxonomy")
    accs = x.annotations['accessions']

    outdict[acc] = [org, tax]
