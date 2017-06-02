"""

USAGE:
python blast.py -t [blast flavour] -q [query file] -s [db/local file to blast against] -e [evalue] -o [outfile]

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

parser = argparse.ArgumentParser(description = 'General purpose BLAST function.', 
                                usage = 'blast.py [options]')

# required args
parser.add_argument('-q', '--query', required = True,
					type = str, help = 'Query sequence.')
parser.add_argument('-s', '--subject', required = True,
					type = str, help = 'Reference database, or .fasta file containing subject sequence.')
parser.add_argument('-t', '--task', required = True,
					type = str, help = 'Flavor of BLAST to be used.')
parser.add_argument('-e', '--evalue', required = True,
					type = float, help = 'Expect value. Default is 10')
parser.add_argument('-o', '--output', required = True,
					type = str, help = 'Name of output file to write to.')

# optional args
parser.add_argument('-a', '--num_alignments', required = False,
					type = int, help = 'Number of alignments. Optional.')
parser.add_argument('-r', '--reward', required = False,
					type = int, help = 'Reward for nucleotide match. Optional.')
parser.add_argument('-p', '--penalty', required = False,
					type = int, help = 'Penalty for nucleotide mismatch. Optional.')

args = parser.parse_args()
task = args.task
outfmt = '-outfmt 5'

# initialize command line + add task
cline = []
cline.append(task)

# query parameter
query = '-query ' + args.query 
evalue = '-evalue ' + str(args.evalue)
output = '-out ' + str(args.output)

# subject/db parameter
subject = os.path.splitext(args.subject)
if subject[1] == '': # db name provided
	subject = '-db ' + subject[0]
elif subject[1] == '.fasta' or subject[1] == '.fa':
	subject = '-subject ' + ''.join(subject)
else:
	sys.exit() 
	print('Error in db/subject specified.')
	print('Please specify either a db name or an input .fasta file.')

## final cline list construction
cline.extend([query, subject, outfmt, evalue, output])

### optional parameters

# num alignments
num_alignments = '-num_alignments ' + str(args.num_alignments)
cline.append(num_alignments)

# rewards and penalties
reward = '-reward ' + str(args.reward)
penalties = '-penalty ' + str(args.penalty)

# cline extension if necessary
cline.extend([reward, penalties])

cline = ' '.join(cline)

print('Command line BLAST+')
print('Options entered:\n')
print(cline, '\n')
subprocess.call(cline)


