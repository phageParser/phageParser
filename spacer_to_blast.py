''' Script to BLAST spacers against host genome '''
from datetime import datetime
from datetime import timedelta
import os
import time


def fetch(accession_id):
    '''Fetch Spacer File'''
    print "Downloading spacers for accession id -- {0}".format(accession_id)


def spacers(accession_id, months=None, blast=None, dir_path='./data/spacers/'):
    ''' gather spacers to blast '''
    if isinstance(accession_id, str):
        pass
    else:
        raise ValueError("Check Accession Number")
    blast = blast
    accession_file = accession_id + '.fasta'
    accession_file_path = os.path.join(dir_path, accession_file)
    # --- Potential Bug:  Directory path needs to be cross platform compatiable
    if is_downloaded(accession_file_path, months):
        pass
    else:
        try:
            fetch(accession_id)
            return
        except:
            raise IOError("Check Accession Number")
    print "You won!, You can proceed to Blasting Spacers"
    return


def is_downloaded(accession_file_path, months):
    ''' Checks if file is downloaded and
    meets time framed requiredment
    Returns True or False
    '''
    file_exists = os.path.exists(accession_file_path)
    if not file_exists:
        return False
    else:
        date_created = time.ctime(os.path.getctime(accession_file_path))
    if months is None:
        genome_uptodate = True
    else:
        genome_uptodate = is_recent(months, date_created)
    if genome_uptodate:
        return True
    else:
        return False
    return


def is_recent(months, date_created):
    '''
    Checks if number of days passsed since file downloaded
    is greater than time constraint specified by user
    Returns: True or False
    '''
    date_created = datetime.strptime(date_created, '%a %b %d %H:%M:%S %Y')
    if isinstance(months, int) or isinstance(months, float):
        # crude method that is not precise
        days = months * 30
        days_required = timedelta(days=days)
        days_passed_since_download = datetime.today() - date_created
        time_constraint = days_passed_since_download > days_required
        return time_constraint
    else:
        raise IOError("Time variable must be an int or float")
    return
