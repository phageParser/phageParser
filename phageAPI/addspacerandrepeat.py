#!/usr/bin/python

'''
Issue #120 - testing
Author: Ahmed

Usage:
python addspacerandrepeat.py

'''


import os
import sys
from Bio import SeqIO

spacerinput = os.listdir('../data/spacers') # these are multiple files
repeatinput = '../data/DRdatabase.txt'

def populate_spacers():
	spacerfile = SeqIO.parse('../data/spacerdatabase.txt', 'fasta') 
	for entry in spacerfile: # that one file in data/
		spacer = entry.seq
		spacer, created = Spacer.objects.get_or_create(sequence = spacer)
	for file in spacerinput: # those files in data/spacers
		filename = '../data/spacers/' + file 
		handle = SeqIO.parse(filename, 'fasta')
		for entry in handle:
			spacer = entry.seq
			spacer, created = Spacer.objects.get_or_create(sequence = spacer) 

def populate_repeats():
	repeatfile = SeqIO.parse(repeatinput, 'fasta')
	for entry in repeatfile:
		repeat = entry.seq
		repeat, created = Repeat.objects.get_or_create(sequence = repeat)

if __name__ == '__main__':
	print("Initializing spacer table population.")
	import os
	import django
	os.environ.setdefault("DJANGO_SETTINGS_MODULE", "phageAPI.settings")
	django.setup()
	from restapi.models import Spacer, Repeat
	populate_spacers()
	print("Spacer table population complete.")
	print("Initializing repeat table population.")
	populate_repeats() 
	print("Repeat table population complete.")

