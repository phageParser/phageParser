import os
import time
import datetime
import functools
import zlib
import requests
import re
from tqdm import tqdm
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
Downloads a file from url to a given path.
Also shows a progress bar to keep track.
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

"""
Checks a given path for file and downloads if a url is given or file name
is an accession id. Also downloads the file if the remote location has a
more recent version. Gives an IOError if path does not exist and an accession id
can not be found.
    Usage:
        fetch('NC_000853.txt')
        fetch('DRDatabase.txt', 'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/DR/DRdatabase')
    Output:
    NC_000853.txt file would be created in the running folder, if it existed before,
    it would be updated to a newer version if available in Entrez database. DRDatabase
    would also be created in the running folder, if the url has a more recent version
    the local file would be updated.
"""
def fetch(path=None, url=None):
    """
    Checks and downloads a url if last modified date does not exist for
    the url or is more recent than local file.
    """
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
    """
    Checks and downloads an accession file from nuccore database if file does
    not exist or its dates are different from url version.
    """
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
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=gbwithparts&tool=biopython&db=nuccore&id={}&email=example%40example.com'.format(acc)
            r = requests.get(url, stream=True)
            if r.status_code != requests.codes.ok:
                print('Status code is faulty from Entrez, returning from gbsync')
                r.close()
                return False
            if not path_exists:
                print('Path given does not exist, downloading from url')
                r.close()
                downfile(path, url)
                return True
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
                    gbdate = matches.groups()[0]
                '''
                Date found in file, download the first chunk of url and
                decompress with gzip to get url date.
                '''
                chunk_size = 256
                fchunk = r.raw.read(chunk_size)
                r.close()
                gzip_decomp = zlib.decompressobj(16+zlib.MAX_WBITS)
                decomp_chunk = gzip_decomp.decompress(fchunk)
                print(decomp_chunk)
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
                print("File {} does not exist and no url is given nor an accession id found, raising IOError".format(filename))
                raise IOError('No file found')
