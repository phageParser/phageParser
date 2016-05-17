'''A decorator to check if file is present or up-to-date'''
import os
import time
from datetime import datetime, timedelta
from Bio import Entrez


def file_check(func):
    '''
    Decorator to check if file
    is present or up-to-date

    decorated function require these kwargs:
       'accession_id' - i.e "NC_000853"
       'months' - Options:
                    1. The number of months since file's last download
                    2. value set to NONE
       'data_dir' - directory path where the data files are stored
    '''
    database = 'nuccore'
    Entrez.email = 'example@example.com'

    def fetch(accession_id, data_dir):
        '''Download File'''
        print "Downloading {0}!".format(accession_id)
        handle = Entrez.efetch(db=database, rettype="gbwithparts",
                               id=accession_id)
        file_name = '{0}.fasta'.format(accession_id)
        file_path = os.path.join(data_dir, file_name)
        with open(file_path, 'wb') as outfile:
            for line in handle:
                outfile.write(line)
        return

    def checked(*args, **kwargs):
        '''Determine if file needs to be
           downloaded or re-downloaded
        '''
        accession_id = kwargs['accession_id']
        data_dir = kwargs['data_dir']
        months = kwargs['months']
        downloaded, date_created = is_downloaded(accession_id, data_dir)
        recent = is_recent(months, date_created, downloaded)
        if downloaded is False or recent is False:
            fetch(accession_id, data_dir)
        return func(*args, **kwargs)

    def is_downloaded(accession_id, data_dir):
        ''' Checks if file is downloaded and
        meets time framed requiredment
        Returns True or False
        '''
        accession_file_path = os.path.join(
            data_dir, '{0}.fasta'.format(accession_id))
        file_exists = os.path.exists(accession_file_path)
        if not file_exists:
            downloaded = False
            date_created = None
        else:
            downloaded = True
            date_created = time.ctime(os.path.getctime(
                accession_file_path))
            print "{1} File exists and created @ {0}".format(
                date_created, accession_id)
        return downloaded, date_created

    def is_recent(months, date_created, downloaded):
        '''
        Checks if number of days passsed since file downloaded
        is greater than time constraint specified by user
        Returns: True or False
        '''
        # Firstly, check if months parameter is defined:
        # Make sure that a non-downloaded file isn't marked as recent
        if downloaded is False:
            return False
        elif months is None:
            return True
        # if above conditions meet then proceed
        date_created = datetime.strptime(
            date_created, '%a %b %d %H:%M:%S %Y')
        if isinstance(months, int):
            # crude method that is not precise
            days = months * 30
            days_required = timedelta(days=days)
            days_passed_since_download = datetime.today() - date_created
            print "Days passed since last download", days_passed_since_download
            recent = days_passed_since_download < days_required
            if not recent:
                print "Requested file is not up-to-date!"
            return recent
        else:
            raise TypeError("Months must be an int")
        return
    return checked
