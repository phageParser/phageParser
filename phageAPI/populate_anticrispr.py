#! /usr/bin/env python

import os
from Bio import SeqIO
import textwrap


def populate(sequences, AntiCRISPR):
    for seq in sequences:
        spacer, _ = AntiCRISPR.objects.get_or_create(
            accession=seq.name,
            sequence=str(seq.seq))
        spacer.save()


def main():
    import argparse

    parser = argparse.ArgumentParser(description=textwrap.dedent("""\
        Import anticrispr sequences into the API DB.

        To use, first get the list of accession numbers from
        https://www.nature.com/articles/nmicrobiol201685. This list is
        available locally in `data/antiCRISPR_accessions.txt`, 

        The script `acc2gb.py` can then be used to download the antiCRISPR
        protein sequence in fasta format, assuming you have NICB access:

            cat data/antiCRISPR_accessions.txt | python acc2gb.py your@email.com protein fasta > anticrispr.txt

        Finally, populate the database with the accession numbers in the
        accession field and the sequences in the sequence field:

            cd phageAPI
            populate_anticrispr.py ../anticrispr.txt
        """),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sequences', metavar='FILE', nargs=1,
                        help='path to sequences file, in fasta format')
    args = parser.parse_args()

    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
    import django
    django.setup()
    from restapi.models import AntiCRISPR

    populate(SeqIO.parse(args.sequences[0], 'fasta'), AntiCRISPR)


if __name__ == '__main__':
    main()
