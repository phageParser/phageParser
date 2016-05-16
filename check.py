from datetime import datetime
from datetime import timedelta
import os
import time
import sys
from Bio import Entrez


class Check(object):
    def __init__(self, accession_id, months=None, dir_path='data/spacers'):
        self.months = months
        self.dir_path = dir_path
        self.accession_id = accession_id + '.fasta'
        self.accession_file_path = os.path.join(dir_path, self.accession_id)
        self.date_created = None
        self.downloaded = False
        self.recent = False

    def fetch(self):
        '''handle = Entrez.efetch(
            db='nuccore', rettype="gbwithparts", id=self.accession_id)'''
        print "fetch file"

    def check(self):
        if self.is_downloaded() is False or self.is_recent() is False:
            return self.fetch()
        else:
            return

    def is_downloaded(self):
        ''' Checks if file is downloaded and
        meets time framed requiredment
        Returns True or False
        '''
        file_exists = os.path.exists(self.accession_file_path)
        if not file_exists:
            self.downloaded = False
            return False
        else:
            self.downloaded = True
            self.date_created = time.ctime(os.path.getctime(
                self.accession_file_path))
        return True

    def is_recent(self):
        '''
        Checks if number of days passsed since file downloaded
        is greater than time constraint specified by user
        Returns: True or False
        '''
        # Firstly, check if months parameter is defined:
        # Make sure that a non-downloaded file isn't marked as recent
        if self.downloaded is False:
            return False
        elif self.months is None:
            return True
        # if above conditions meet then proceed
        date_created = datetime.strptime(
            self.date_created, '%a %b %d %H:%M:%S %Y')
        if isinstance(self.months, int):
            # crude method that is not precise
            days = self.months * 30
            days_required = timedelta(days=days)
            days_passed_since_download = datetime.today() - date_created
            self.recent = days_passed_since_download > days_required
            return self.recent
        else:
            raise TypeError("Months must be an int")
        return
