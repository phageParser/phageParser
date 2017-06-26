# filterByExpect_all_v2

# for each xml output file, use NCBIXML to extract the important things
# Note: this requires blast output in xml format (set 'outfmt = 5')

'''
USAGE:
python filterByExpect.py <indir> <outdir>
'''

import os
import csv
import sys
import pandas as pd

def parse_blast(resultfile):  # takes in the BLAST result, outputs list that can be made into csv
    from Bio.Blast import NCBIXML
    result_handle = open(resultfile)
    blast_records = NCBIXML.parse(result_handle)
    csv_list = []

    header = ['Query',
              'Name', 'Length', 'Score', 'Expect',
              'QueryStart', 'QueryEnd',
              'SubjectStart', 'SubjectEnd'
              ]

    # csv_list.append(header)
    count = 0
    for blast_record in blast_records:
        '''help(blast_record.alignments[0].hsps[0])'''  # these give help info for the parts
        '''help(blast_record.alignments[0])        '''
        count += 1

        query = blast_record.query
        for alignment in blast_record.alignments:

            name = alignment.title
            length = alignment.length

            # I don't know if we will ever have more than one, so might as well
            # take the first one.
            hsp = alignment.hsps[0]
            score = hsp.score
            expect = hsp.expect
            querystart = hsp.query_start
            queryend = hsp.query_end
            subjectstart = hsp.sbjct_start
            subjectend = hsp.sbjct_end
            row = [query, name, length, score, expect,
                   querystart, queryend, subjectstart, subjectend]
            csv_list.append(row)

    result_handle.close()
    return pd.DataFrame(csv_list, columns=header)


def write_csv(dest, df):  # takes a list of lists object with each csv row as a list
    df.to_csv(dest)


def main():
    indir = str(sys.argv[1])
    outdir = str(sys.argv[2])

    for fn in os.listdir("%s/" % indir):
        ID = fn[:fn.index('.')]
        if fn == 'sorted':
            continue
        csv_list = parse_blast(fn)
        write_csv("%s/%s.csv" % (outdir, ID), csv_list)


if __name__ == '__main__':
    main()
