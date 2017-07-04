"""
usage - python blastselftargets.py [options]

Performs blastn on each organism for its spacers and writes results in xml format.

python blastselftargets.py -h will throw up arg references if need be.

Author: @cemyld

DEPENDENCIES:
Biopython
blast+
"""
import re
import os
import subprocess
import glob
from util.fetch import fetch
from tqdm import tqdm
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Blast self targeting spacers and output the files.',
                                 usage='blast.py [options]')
# optional args
parser.add_argument('-i', '--gbdir', required=False, type=str,
                    help='Directory of genbank files', default='gbdir')
parser.add_argument('-o', '--outdir', required=False, type=str,
                    help='Output directory path', default='gbdir/blastoutput')

args = parser.parse_args()



def main(gbdir, outdir):
    os.makedirs(gbdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    tempq = 'tempquery.fasta'
    tempdb = 'tempdb.fasta'
    for org in tqdm(Organism.objects.all()):
        # get genbank and convert to fasta
        fpath = os.path.join(gbdir, '{}.gb'.format(org.accession))
        if not os.path.isfile(fpath):
            print('\nFetching {} with accession {}'.format(
                org.name, org.accession))
            fetch(fpath)
        SeqIO.convert(fpath, 'genbank', tempdb, 'fasta')
        # get spacers of organism and convert to fasta
        spacers = Spacer.objects.filter(loci__organism=org)
        fastatext = ''.join(['>{}\n{}\n'.format(spacer.id, spacer.sequence) for spacer in spacers])
        with open(tempq, 'w') as f:
            f.write(fastatext)
        # run blast and save output
        outpath = os.path.join(outdir, '{}.json'.format(org.accession))
        commandargs = ['blastn', '-query', tempq,
                       '-subject', tempdb, '-out', outpath, '-outfmt', '15']
        subprocess.run(
            commandargs, stdout=subprocess.DEVNULL)

    os.remove(tempq)
    os.remove(tempdb)
if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism, Locus, Spacer
    # main(args.gbdir, args.outdir)
    main('gbfiles', 'gbfiles/blastoutput')
