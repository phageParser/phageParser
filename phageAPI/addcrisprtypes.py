"""Running this script requires HMMER tool for sequence analysis
HMMER can be downloaded via www.hmmer.org
A compressed database of profiles is also needed where its path
can be specified as the hmmdbpath argument in hmmscangenbankfiles function."""
import re
import os
import subprocess
import glob
from util.fetch import fetch


def fetchgenbankfiles(fpath='gbfiles'):
    os.makedirs(fpath, exist_ok=True)
    for org in Organism.objects.all():
        print('\nFetching {} with accession {}'.format(org.name, org.accession))
        fetch(os.path.join(fpath, '{}.gb'.format(org.accession)))


def convertgenbanktofasta(fpath):
    print(fpath)
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


def hmmscangenbankfiles(fpath='gbfiles', hmmdbpath='hmmdb/testdb'):
    hmmoutpath = os.path.join(fpath, 'hmmeroutput')
    os.makedirs(hmmoutpath, exist_ok=True)
    for f in glob.glob(fpath + '/*.gb'):
        accession = os.path.splitext(os.path.split(f)[1])[0]
        fastainput = convertgenbanktofasta(f).encode('utf-8')
        table_path = os.path.join(hmmoutpath, '{}.txt'.format(accession))
        commandargs = ['hmmscan', '--noali',
                       '--tblout', table_path, hmmdbpath, '-']
        result = subprocess.run(commandargs, input=fastainput)
if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import Organism
    fetchgenbankfiles()
    hmmscangenbankfiles()
