#!/usr/bin/env python
"""Fetch GenBank entries for given accessions. 

Adapted from https://www.biostars.org/p/66921/

DEPENDENCIES:
Biopython
"""

import sys
import textwrap
import argparse
from Bio import Entrez


def main():
    parser = argparse.ArgumentParser(
        usage='cat INPUT | python acc2gb.py EMAIL DB RETTYPE > OUTPUT',
        description=textwrap.dedent("""\
            The input file should contain accession IDs to download, one per line.

            EXAMPLE:

            cat data/antiCRISPR_accessions.txt | python acc2gb.py your@email.com protein fasta > outfile.txt

            Case 1: rettype = gbwithparts, db = nuccore - downloads genbank
                file with metadata and fasta DNA sequence (i.e. for downloading
                bacterial genomes with metadata)

            Case 2: rettype = fasta, db = protein - downloads fasta file with
                protein sequence (i.e. for downloading antiCRISPR protein sequences
                for BLAST)
            """),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('email', nargs=1,
                        help='your email address (does not need to be registered, just used to identify you)')
    parser.add_argument('db', nargs=1,
                        help='the NCBI database ID, which must be a valid Entrez database name')
    parser.add_argument('rettype', nargs=1,
                        help='the type of file to retrieve')
    args = parser.parse_args()

    Entrez.email = args.email

    # get accession numbers out of stdin
    accs = [l.strip() for l in sys.stdin if l.strip()]

    # fetch
    sys.stderr.write("Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
    for i, acc in enumerate(accs):
        try:
            sys.stderr.write(" %9i %s          \r" % (i + 1, acc))
            handle = Entrez.efetch(db=args.db, rettype=args.rettype, retmode="text", id=acc)
            # print output to stdout
            sys.stdout.write(handle.read())
        except:
            sys.stderr.write("Error! Cannot fetch: %s        \n" % acc)


if __name__ == "__main__":
    main()
