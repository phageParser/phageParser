'''
Test that decorator can verify genome is downloaded
and up-to-date
'''
from check import file_check
from datetime import datetime
import os


@file_check
def print_date(*args, **kwargs):
    '''
    print datetime after decorator
    checks if file exists and is
    up-to-date
    '''
    print "Date printed Outside of decorator: ", datetime.now()

if __name__ == '__main__':
    if os.path.exists('data/genomes/NC_000858.txt'):
        os.remove('data/genomes/NC_000858.txt')
    print("\nTest 1: File does not exists")
    print("----------------------------")
    print_date(accession_id='NC_000858', months=2)
    print("\nTest 2: File exists")
    print("-------------------")
    print_date(accession_id='NC_000858', months=2)
    print("\n")
