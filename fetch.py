""" A decorator to check if genome file exists or up-to-date.

    Example:
        The decorated function require the 'accession_id' kwarg
        and optionally the 'months' kwarg::

        >>> from check import genome_check
        >>> from datetime import datetime

        >>> @genome_check
            def print_date(*args, **kwargs):
               print "Printed by func 'print_date' ", datetime.now()

        >>> print_date(accession_id='NC_000853', months=8)

    Output:
        NC_000853 File exists and created @ Tue May 17 12:48:14 2016
        Days passed since last download 0:00:32.932481
        Printed by func 'print_date' 2016-05-17 12:48:46.932501

"""

import os
import time
import datetime
from Bio import Entrez
import functools
import requests
from tqdm import tqdm
import re
import zlib
"""
Decorator for downloading or updating files required for a function.
"""
def reqfile(func=None, path=None, url=None):
    if not func:
        return functools.partial(fetch_factory, path=path, url=url)
    @functools.wraps(test_func)
    def wrapper(*args, **kwargs):
        fetch(path, url)
        return func(*args, **kwargs)
    return wrapper
"""
Checks a given path for file and downloads if a url is given or file name
is an accession id. Also downloads the file if the remote location has a
more recent version.
"""
"""
Downloads a file to a given path. Also shows a progress bar
"""
def downfile(path, url):
    r = requests.get(url, stream=True)
    r.raise_for_status()
    total_size = int(r.headers.get('content-length', 0))
    with open(path, 'wb') as f:
        chunk_size = 64*1024
        pbar = tqdm(total=total_size, unit='B', unit_scale=True)
        for data in r.iter_content(chunk_size):
            f.write(data)
            pbar.update(chunk_size)


def fetch(path=None, url=None):

    def sync():
        if not path_exists:
            downfile(path, url)
            return
        '''
        Check last modified dates of file and url, download if url is newer.
        '''
        filemodtime = datetime.datetime.fromtimestamp(os.path.getmtime(path))
        r = requests.get(url, stream=True)
        if 'Last-Modified' not in r.headers.keys():
            print('Last modified header not found in url, downloading...')
            downfile(path, url) #no last-modified header in url, downloading file
        urlmodstr = r.headers['Last-Modified']
        urlmodtime = datetime.datetime.strptime(urlmodstr, '%a, %d %b %Y %H:%M:%S %Z')
        if filemodtime < urlmodtime:
            '''
            Url file is more recent, downloading url.'''
            print('Url is more recent than file, downloading')
            downfile(path, url)
            return
        print('File is recent, returning')
        return
    def gbsync():
        print('Trying to fetch from Entrez')
        regex = r'(\w{2}_\d{6}|\w{1}_\d{5})'
        filename = os.path.basename(path)
        matches = re.search(regex, filename)
        if not matches:
            print('Filename does not match an accession, returning from gbsync')
            return False
        else:
            acc = matches.groups()[0]
            url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=gbwithparts&tool=biopython&db=nuccore&id={}&email=example%40example.com'.format(acc)
            if not path_exists:
                print('Path given does not exist, downloading from url')
                downfile(path, url)
            else:
                '''
                Path exists, try to get date to compare with Entrez version, download if different.
                '''
                regex = r'(\d{2}-\w{3}-\d{4})'
                with open(path, 'r') as f:
                    print('Checking accession file for date')
                    fline = f.readline()
                    matches = re.search(regex, fline)
                    if not matches:
                        print('No date found in accession file, returning from gbsync')
                        return False
                    else:
                        gbdate = matches.groups()[0]
                        r = requests.get(url, stream=True)
                        if r.status_code != requests.codes.ok:
                            print('Status code is faulty from Entrez, returning from gbsync')
                            return False
                        else:
                            '''
                            Status code ok, download the first chunk and
                            decompress with gzip to get date.
                            '''
                            chunk_size = 256
                            fchunk = r.raw.read(chunk_size)
                            r.close()
                            gzip_decomp = zlib.decompressobj(16+zlib.MAX_WBITS)
                            decomp_chunk = gzip_decomp.decompress(fchunk)
                            urldate = re.search(regex, decomp_chunk).groups()[0]
                            print('Status code ok from Entrez, checking dates')
                            if gbdate != urldate:
                                print("Dates don't match, downloading")
                                downfile(path, url)
                                return True
                            else:
                                print('Dates are matching, returning from gbsync')
                                return True

    if not path:
        print('No path given, returning')
        return
    filename = os.path.basename(path)
    path_exists = os.path.isfile(path)
    if url:
        sync()
    else:
        '''
        Try to fetch the file as an accession id, if it fails print status.
        '''
        if not gbsync():
            if path_exists:
                print("File {} exists but no url is given nor an accession id found, returning".format(filename))
                return
            else:
                print("File {} does not exist and no url is given nor an accession id found, returning".format(filename))
                raise IOError('No file found')
