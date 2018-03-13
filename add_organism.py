
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
    
    print(o)
    # Avoid Duplicate Entries (Comment out for testing)
    if not created:
        return created

    # == for testing only ==
    seq = get_pickle(accession)
    if not seq:
    # == for testing only ==
        seq = get_ncbi_seq(Entrez_email, 'nuccore', 'fasta', accession)

    with open('seq.fasta', 'w') as file:
        file.write(seq)

    if os.path.isdir(CRISPRFinder_output_folder):
        subprocess.Popen(['rm', '-rf', CRISPRFinder_output_folder + '/'], stdout = sys.stdout)

    crisprfinder_script = subprocess.Popen(['perl' ,'CRISPRCasFinder/CRISPRCasFinder.pl','-out', CRISPRFinder_output_folder, '-in', 'seq.fasta'], stdout=sys.stdout)
    crisprfinder_script.communicate()

    results = json.load(open(CRISPRFinder_output_folder + '/result.json'))

    for seq in results['Sequences']:
        for cspr in seq['Crisprs']:
            loc_start = int(cspr['Start'])
            loc_end = int(cspr['End'])
            dr_consensus = cspr['DR_Consensus']
            orientation = cspr['Potential_Orientation'] == '+'
            locus, _ = Locus.objects.get_or_create(
                organism=o,
                genomic_start=loc_start,
                genomic_end=loc_end,
                consensus=dr_consensus,
                orientation=orientation
            )
            for i in range(1, len(cspr['Regions'])-2, 2):
                rep = cspr['Regions'][i]
                spac = cspr['Regions'][i+1]
                if rep['Type'] != 'DR':
                    print('First region in repeat/spacer pair is not a repeat!')
                    print('Organism: ', o)
                    return
                elif spac['Type'] != 'Spacer':
                    print('Second region in repeat/spacer pair is not a spacer!')
                    print('Organism: ', o)
                    return

                order = (i-1)/2 + 1 #starting at 1
                repeat, _ = Repeat.objects.get_or_create(sequence=rep['Sequence'])
                spacer, _ = Spacer.objects.get_or_create(sequence=spac['Sequence'])

                lsr, _ = LocusSpacerRepeat.objects.get_or_create(
                    locus=locus,
                    spacer=spacer,
                    repeat=repeat,
                    order=order
                )
            rep = cspr['Regions'][-2]
            if rep['Type'] != 'DR':
                print('Last region in locus is not a repeat!')
                print('Organism: ', o)
                return
            order = (len(cspr['Regions'])-3)/2 + 1
            repeat, _ = Repeat.objects.get_or_create(sequence=rep['Sequence'])
            lsr, _ = LocusSpacerRepeat.objects.get_or_create(
                locus=locus,
                spacer=None,
                repeat=repeat,
                order=order
            )

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