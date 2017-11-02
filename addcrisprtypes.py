"""Running this script requires HMMER tool for sequence analysis
HMMER can be downloaded via www.hmmer.org
A compressed database of profiles is also needed where its path
can be specified as the hmmdbpath argument in hmmscangenbankfiles function."""
import glob
import os
import re
import subprocess

from tqdm import tqdm

from util.fetch import fetch


def fetch_genbank_files(gbdir='gbfiles'):
    os.makedirs(gbdir, exist_ok=True)
    for org in tqdm(Organism.objects.all()):
        fpath = os.path.join(gbdir, '{}.gb'.format(org.accession))
        if os.path.isfile(fpath):
            continue
        print('\nFetching {} with accession {}'.format(
            org.name,
            org.accession
        ))
        fetch(os.path.join(fpath, '{}.gb'.format(org.accession)))


def convert_genbank_to_fasta(fpath):
    with open(fpath, 'r') as f:
        file_text = f.read()

    def extract_cds():
        cds = []
        m = re.finditer('CDS\s{13}(?P<CDS>\S*)', file_text)

        for reg in m:
            element = reg.group('CDS')
            if element != "":
                cds.append(element)

        return cds

    def extract_translations():
        translations = []
        m = re.finditer(
            '/translation="(?P<translation>[\w\n ]+)"',
            file_text,
            re.MULTILINE
        )
        for reg in m:
            translation = reg.group('translation')
            translation = translation.replace("\n", "")
            translation = translation.replace(" ", "")
            translations.append(translation)

        return translations

    cds = extract_cds()
    translations = extract_translations()

    fastatext = ''.join([">{}\n{}\n".format(a, b)
                         for a, b in zip(cds, translations)])
    return fastatext


def hmmscan_genbank_files(gbdir='gbfiles', hmmdbpath='util/hmmercasdb/casdb'):
    hmm_outpath = os.path.join(gbdir, 'hmmeroutput')
    os.makedirs(hmm_outpath, exist_ok=True)
    print('Running hmmerscan on files in directory {}'.format(gbdir))
    for f in glob.glob(gbdir + '/*.gb'):
        accession = os.path.splitext(os.path.split(f)[1])[0]
        fastainput = convert_genbank_to_fasta(f).encode('utf-8')
        table_path = os.path.join(hmm_outpath, '{}.txt'.format(accession))
        commandargs = ['hmmscan', '--noali',
                       '--tblout', table_path, hmmdbpath, '-']
        subprocess.run(
            commandargs,
            input=fastainput,
            stdout=subprocess.DEVNULL
        )


if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django

    django.setup()
    from restapi.models import Organism

    fetch_genbank_files()
    hmmscan_genbank_files()
