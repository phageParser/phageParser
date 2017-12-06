# -*- coding: utf-8 -*-
'''
Created on Tue Aug 15 11:26:23 2017
This script downloads from .ids files from NCBI database
(from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/)
to generate lists of organism accession IDs, later used for downloading
complete genome sequences
@author: mbonsma / thisisjaid / cemyld
'''

import pandas as pd
import argparse

ncbi_ftp_url = "ftp.ncbi.nlm.nih.gov"
ftp_path = "/genomes/GENOME_REPORTS/IDS/"


def idslist_to_dataframe(idslist):
    if len(idslist) < 1:
        print("No ids collection given")
        return ""
    df_list = []
    for ids in idslist:
        try:
            _url = 'ftp://' + ncbi_ftp_url + ftp_path + ids
            df_list.append(pd.read_csv(_url, sep='\t', header=None))
        except Exception as e:
            # might not find the filename in ftp directory
            raise e

    return pd.concat(df_list)


def ids_to_acclist(idslist):
    return list(idslist_to_dataframe(idslist)[1])


def main():
    accessions = ids_to_acclist(args.ids_file)
    for acc in accessions:
        print(acc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=True, description='''phageParser - collect_accessions.py -
                                     This script downloads from .ids files from NCBI database
                                     (from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/)
                                     to generate lists of organism accession IDs, later used for downloading
                                     complete genome sequences''')
    parser.add_argument('ids_file', action='store', type=str,
                        help='Ids files to be downloaded, a full listing can be found at ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/', nargs='+')

    args = parser.parse_args()
    main()
