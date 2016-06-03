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
from datetime import datetime, timedelta
from Bio import Entrez


def genome_check(func):
    """
    Decorator to check if genome file exists or`` up-to-date

    decorated function require these kwargs:
       'accession_id' - i.e "NC_000853"
       'months' - Options:
                    1. The number of months since file's last download
                    2. value set to NONE
    """

    database = 'nuccore'
    """ Database to download genome files from."""

    Entrez.email = 'example@example.com'
    """ Email address needed by NCBI to use Entrez """

    data_dir = 'data/genomes/'
    """ Directory where genome file will be saved """

    def fetch(accession_id):
        """Download and write genome file

        Args:
            accession_id (str): The genome accession ID.

        Returns:
            file:  Writes text file to genome folder.
        """

        print "Downloading {0}!".format(accession_id)
        handle = Entrez.efetch(db=database, rettype="gbwithparts",
                               id=accession_id)
        file_name = '{0}.txt'.format(accession_id)
        file_path = os.path.join(data_dir, file_name)
        with open(file_path, 'wb') as outfile:
            for line in handle:
                outfile.write(line)
        return

    def checked(*args, **kwargs):
        """Determine if file needs to be downloaded or re-downloaded.

        Args:
            **kwargs:
                accession_id (str): The genome Accession Id
                months (Optional[int]): Parameter used to determine
                    if genome is up-to-date.  Defaults to None

        Returns:
            Decorated function::

                func(*args, **kwargs)
        """
        accession_id = kwargs['accession_id']
        if 'months' not in kwargs:
            months = None
        else:
            months = kwargs['months']
        downloaded, date_created = is_downloaded(accession_id)
        recent = is_recent(months, date_created, downloaded)
        if downloaded is False or recent is False:
            fetch(accession_id)
        return func(*args, **kwargs)

    def is_downloaded(accession_id):
        """Checks if file is downloaded and specifies file creation date

        Args:
            accession_id (str): Genome accession ID

        Returns:
            bool: True if downloaded, False otherwise
            timestamp: Time is in ctime format
        """

        accession_file_path = os.path.join(
            data_dir, '{0}.txt'.format(accession_id))
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
        """ Checks if genome file meets user-specified time constraints.

        Args:
            months: (int or None): Number of months since last
                download of genome file or 'None' if not specified.
            date_created: (ctime timestamp): Date file was created.
            downloaded: (bool): True if file exists, otherwise false.

        Returns:
            bool:  True if number of days passed is less than
                time constraints, otherwise False
        """
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
