"""Running this script requires HMMER tool for sequence analysis
HMMER can be downloaded via www.hmmer.org
A compressed database of profiles is also needed where its path
can be specified as the hmmdbpath argument in hmmscangenbankfiles function."""
import re
import os
import subprocess
import glob
from util.fetch import fetch
from tqdm import tqdm


def fetchgenbankfiles(gbdir='gbfiles'):
    os.makedirs(gbdir, exist_ok=True)
    for org in tqdm(Organism.objects.all()):
        fpath = os.path.join(gbdir, '{}.gb'.format(org.accession))
        if os.path.isfile(fpath):
            continue
        print('\nFetching {} with accession {}'.format(org.name, org.accession))
        fetch(os.path.join(fpath, '{}.gb'.format(org.accession)))


def convertgenbanktofasta(fpath):
    with open(fpath, 'r') as f:
        fileText = f.read()

    def extractCDS():
        cds = []
        m = re.finditer('CDS\s{13}(?P<CDS>\S*)', fileText)

        for reg in m:
            element = reg.group('CDS')
            if element != "":
                cds.append(element)

        return cds

    def extractTranslations():
        translations = []
        m = re.finditer(
            '/translation="(?P<translation>[\w\n ]+)"', fileText, re.MULTILINE)
        for reg in m:
            translation = reg.group('translation')
            translation = translation.replace("\n", "")
            translation = translation.replace(" ", "")
            translations.append(translation)

        return translations

    cds = extractCDS()
    translations = extractTranslations()

    fastatext = ''.join([">{}\n{}\n".format(a, b)
                         for a, b in zip(cds, translations)])
    return fastatext


def hmmscangenbankfiles(gbdir='gbfiles', hmmdbpath='util/hmmercasdb/casdb'):
    hmmoutpath = os.path.join(gbdir, 'hmmeroutput')
    os.makedirs(hmmoutpath, exist_ok=True)
    print('Running hmmerscan on files in directory {}'.format(gbdir))
    for f in glob.glob(gbdir + '/*.gb'):
        accession = os.path.splitext(os.path.split(f)[1])[0]
        fastainput = convertgenbanktofasta(f).encode('utf-8')
        table_path = os.path.join(hmmoutpath, '{}.txt'.format(accession))
        commandargs = ['hmmscan', '--noali',
                       '--tblout', table_path, hmmdbpath, '-']
        subprocess.run(
            commandargs, input=fastainput, stdout=subprocess.DEVNULL)

if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism
    fetchgenbankfiles()
    hmmscangenbankfiles()
