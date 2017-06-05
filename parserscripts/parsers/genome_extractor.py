import csv

''' 
@author:         Nigesh Shakya(neaGaze)
@timestamp:      06/08/2016 @ 12:54am (UTC)
    
The Class file that examines the file downloaded via prefetch (see #2), and returns a string corresponding to the genome starting n bases before startPosition and carrying on to m bases after endPosition. If either n or m are not provided, default them to 100.
This patch pertains to issue #3 in the github wiki (https://github.com/goyalsid/phageParser/issues/3)
'''
class GenomeExtractor:

    def __init__(self, infile, **kwargs):
        blast_file = open(infile, 'r')
        self.reader = csv.reader(blast_file, dialect=csv.excel_tab)
        pass

    '''  Genome Extractor '''
    def findNeighbours(self, startPosition, endPosition, n=100, m=100):
        count = 0
        hasDataFound = False
        matchedGene = ""
        for row in self.reader:

            if row[0] == "//":
                break

            if row[0].startswith("ORIGIN"): 
                hasDataFound = True
                continue

            if hasDataFound:
                rowVals = row[0].split()
                pos = int(rowVals[0])
                if pos >= startPosition -n and pos < endPosition + m:
                    matchedGene += " ".join(rowVals)
                    matchedGene += "\n"
                     
        return matchedGene


