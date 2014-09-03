import re

phageDB = 'PhagesDB_Data.txt'

phageSource = []

# open file
ins = open( phageDB, "r" )

for line in ins:	
	# split on tab
	splittedLines = re.split(r'\t', line)

    # push to array
	phageSource.append(splittedLines)    
ins.close()

# remove the first element because this is the header
phageSource.pop(0)

def searchAccessionByPhageName(PhageName):

    for line in phageSource:

        # does the phagename exist in the array?
        if PhageName in line:

            # is In genbank true?
            if line[4] == 'True':
                return line[5]


# testing this out...
# print(searchAccessionByPhageName('Angel'))