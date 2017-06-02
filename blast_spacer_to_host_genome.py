''' Script to BLAST spacers against host genome '''
import os
import sys


def fetch(accession_id):
    '''Fetch Spacer File'''
    print("Downloading spacers for accession id -- {0}".format(accession_id))


def spacers(accession_id, months=None, blast=None, dir_path='./data/spacers/'):
    ''' gather spacers to blast '''
    if isinstance(accession_id, str):
        pass
    else:
        raise ValueError("Check Accession Number")
    blast = blast
    accession_file = accession_id + '.fasta'
    accession_file_path = os.path.join(dir_path, accession_file)
    with open(accession_file_path, 'rb') as reader:
        for line in reader:
            print(line)
            break

def main():
    """ Main script """
    spacers('NC_000853')
    return

if __name__ == '__main__':
    sys.exit(main())
