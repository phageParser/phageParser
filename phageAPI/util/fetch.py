import os
import time
import datetime
import functools
import zlib
import requests
import re
from tqdm import tqdm


def download(path, url, chunk_size=64 * 1024):
    """
    Downloads a file from url to a given path.

    Opens a connection to the url and downloads the stream by supplied
    chunk size while showing a progress bar in command line. Overwrites
    any file existing in the path. Can be used in other functions to download
    files without checking.
        Args:
        path (str): Path to file to be downloaded.
        url (str): Url to download the file from.
        chunk_size (int, optional): Chunk size to download by. Defaults to 64*1024.
        Notes:
        chunk_size should be a multiple of 2 for faster downloads.
        Returns:
            bool: True if the file is successfully downloaded. False otherwise.
        Examples
            >>>download('DRDatabase.txt', 'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/DR/DRdatabase')
            True
            >>>download('DRDatabase.txt', 'faulty url')
            False
    """
    try:
        r = requests.get(url, stream=True)
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        with open(path, 'wb') as f:
            chunk_size = 64 * 1024
            pbar = tqdm(total=total_size, unit='B', unit_scale=True)
            for data in r.iter_content(chunk_size):
                f.write(data)
                pbar.update(chunk_size)
        return True
    except Exception as e:
        print(e)
        return False


def fetch(path, url=None):
    """Fetches a file from either url or accession id in filename, updates file if local version is older.

    Checks a given path for file and downloads if a url is given or file name
    is an accession id. Also downloads the file if the remote location has a
    more recent version. Can be used in other functions to download
    required files from url or accession gb files.
        Args:
        path (str): Path to the file to be fetched
        url (str, optional): Url to update the local file from. Defaults to None.
        Returns:
            bool: True if the file is updated or downloaded. False otherwise.
        Examples:
            >>>import os
            >>>os.remove('NC_000853.txt')
            >>>fetch('NC_000853.txt')
            True
            >>>import os
            >>>os.remove('DRDatabase.txt')
            >>>fetch('DRDatabase.txt', 'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/DR/DRdatabase')
            True
        Result:
        NC_000853.txt file would be created in the running folder, if it existed before,
        it would be updated to a newer version if available in Entrez database. DRDatabase
        would also be created in the running folder, if the url has a more recent version
        the local file would be updated.
    """
    def sync():
        """
        Checks and downloads a url if last modified date does not exist for
        the url or is more recent than local file.
        """
        if not path_exists:
            return download(path, url)
        '''Check last modified dates of file and url, download if url is newer.'''
        filemodtime = datetime.datetime.fromtimestamp(os.path.getmtime(path))
        r = requests.get(url, stream=True)
        if 'Last-Modified' not in r.headers:
            print('Last modified header not found in url, downloading...')
            # no last-modified header in url, downloading file
            return download(path, url)
        urlmodstr = r.headers['Last-Modified']
        urlmodtime = datetime.datetime.strptime(
            urlmodstr, '%a, %d %b %Y %H:%M:%S %Z')
        if filemodtime < urlmodtime:
            '''Url file is more recent, downloading url.'''
            print('Url is more recent than file, downloading')
            return download(path, url)
        print('File is up-to-date')
        return False

    def gbsync():
        """
        Checks and downloads an accession file from nuccore database if file does
        not exist or its dates are different from url version.
        """
        print('Trying to fetch from Entrez')
        regex = r'([A-Z]{1,2}_\w+)'
        filename = os.path.basename(path)
        matches = re.search(regex, filename)
        if not matches:
            print('Filename does not match an accession')
            return False
        else:
            acc = matches.groups()[0]
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=gbwithparts&tool=biopython&db=nuccore&id={}&email=example%40example.com'.format(
                acc)
            r = requests.get(url, stream=True)
            if r.status_code != requests.codes.ok:
                print('Bad Status code returned from Entrez')
                r.close()
                return False
            if not path_exists:
                print('Path given does not exist, downloading from url')
                r.close()
                return download(path, url)
            else:
                '''Path exists, try to get date to compare with Entrez version, download if different.'''
                regex = r'(\d{2}-\w{3}-\d{4})'
                with open(path, 'r') as f:
                    print('Checking accession file for date')
                    fline = f.readline()
                    matches = re.search(regex, fline)
                    if not matches:
                        print(
                            'No date found in accession file {}, overwriting with Entrez entry'.format(filename))
                        return download(path, url)
                    gbdate = matches.groups()[0]
                '''
                Date found in file, download the first chunk of url and
                decompress with gzip to get url date.
                '''
                chunk_size = 256
                fchunk = r.raw.read(chunk_size)
                r.close()
                gzip_decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
                decomp_chunk = gzip_decomp.decompress(fchunk).decode()
                urldate = re.search(regex, decomp_chunk).groups()[0]
                if gbdate != urldate:
                    print(
                        "Dates don't match for accession file {}, downloading".format(filename))
                    return download(path, url)
                else:
                    print(
                        'Dates are matching for accession file {} returning'.format(filename))
                    return False

    filename = os.path.basename(path)
    path_exists = os.path.isfile(path)
    if url:
        return sync()
    else:
        '''Try to fetch the file as an accession id, if it fails return False.'''
        return gbsync()
