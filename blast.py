"""
usage - python blast.py [options]

Performs input blast and writes results in xml format.

Options for -t: blastn, blastp, tblastn, blastx, tblastx, psiblast
python blast.py -h will throw up arg references if need be.

Author: @aays

DEPENDENCIES:
Biopython
blast+
"""

import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='General purpose BLAST function.',
                                 usage='blast.py [options]')

# required args
parser.add_argument('-q', '--query', required=True, type=str,
                    help='Query sequence.')
parser.add_argument('-s', '--subject', required=True, type=str,
                    help='Ref database name or .fasta file to blast against.')
parser.add_argument('-t', '--task', required=True, type=str,
                    help='Flavor of BLAST to be used.')
parser.add_argument('-e', '--evalue', required=True, type=float,
                    help='Expect value. Default is 10')
parser.add_argument('-o', '--output', required=True, type=str,
                    help='Name of output file to write to.')

# optional args
parser.add_argument('-a', '--num_alignments', type=int,
                    help='Number of alignments. Optional.')
parser.add_argument('-r', '--reward', type=int,
                    help='Reward for nucleotide match. Optional.')
parser.add_argument('-p', '--penalty', type=int,
                    help='Penalty for nucleotide mismatch. Optional.')

args = parser.parse_args()
task = args.task
outfmt = '-outfmt 5'

# initialize command line + add task
cline = []
cline.extend([task, outfmt])

cline.extend(['-query', args.query])
cline.extend(['-evalue', str(args.evalue)])
cline.extend(['-out', str(args.output)])

# subject/db parameter
subject = os.path.splitext(args.subject)
if subject[1] == '':  # db name provided
    cline.extend(['-db', subject[0]])
elif subject[1] in ('.fasta', '.fa'):
    cline.extend(['-subject', args.subject])
else:
    sys.exit()
    print('Error in db/subject specified.')
    print('Please specify either a db name or an input .fasta file.')

# optional parameters
cline.extend(['-num_alignments', str(args.num_alignments)])
cline.extend(['-reward', str(args.reward)])
cline.extend(['-penalty', str(args.penalty)])

print('Command line BLAST+')
print('Options entered:\n')
print(' '.join(cline), '\n')
subprocess.call(cline)
