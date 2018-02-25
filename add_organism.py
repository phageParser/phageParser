
import sys
import pprint 
import os
import pickle
import subprocess
import json
from pprint import pprint

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
import django
django.setup()

from restapi.models import (
    Organism,
    Spacer,
    Repeat,
    LocusSpacerRepeat,
    AntiCRISPR,
    Locus
)

### Constants ###
Entrez_email = 'wenxinchen11@hotmail.com'
CRISPRFinder_output_folder = 'output_files'

def pickle_item(item, filename):

    with open(filename, 'wb') as handle:
        pickle.dump(item, handle, protocol=pickle.HIGHEST_PROTOCOL)

def get_pickle(filename):

    f = None

    if os.path.isfile(filename):
        with open(filename, 'rb') as handle:
            f = pickle.load(handle)
    return f

def add_sequence_to_db(accession, name):

    o, created = Organism.objects.get_or_create(
        name=name, accession=accession)
    
    # Avoid Duplicate Entries (Comment out for testing)
    if not created:
        return created

    # for testing only
    seq = get_pickle(accession)

    if not seq:
        seq = get_ncbi_seq(Entrez_email, 'nuccore', 'fasta', accession)

    with open('seq.fasta', 'w') as file:
        file.write(seq)

    if os.path.isdir(CRISPRFinder_output_folder):
        subprocess.Popen(['rm', '-rf', CRISPRFinder_output_folder + '/'], stdout = sys.stdout)

    crisprfinder_script = subprocess.Popen(['perl' ,'CRISPRCasFinder/CRISPRCasFinder.pl','-out', CRISPRFinder_output_folder, '-in', 'seq.fasta'], stdout=sys.stdout)
    crisprfinder_script.communicate()

    results = json.load(open(CRISPRFinder_output_folder + '/result.json'))

    for seq in results['Sequences']:
        pprint(results)

def get_ncbi_seq(email, db, rettype, accession):

    # fetch
    print("Fetching accession %s from GenBank\n" % (accession))

    Entrez.email = email

    try:
        handle = Entrez.efetch(
            db=db,
            rettype=rettype,
            retmode="text",
            id=accession
        )
        res = handle.read()

        # for testing only
        pickle_item(res, accession)
        return res

    except Exception:
        sys.stderr.write("Error! Cannot fetch: %s        \n" % accession)

if __name__ == '__main__':
    add_sequence_to_db(accession='NC_000853', name='Thermotoga maritima MSB8')