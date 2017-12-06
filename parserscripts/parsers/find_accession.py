import csv


class PhageFinder:
    PHAGE_NAME = 0
    CLUSTER = 1
    EXIST_IN_GENBANK = 4
    ACCESSION = 5

    def __init__(self, infile, **kwargs):

        blast_file = open(infile, 'r')
        self.reader = csv.reader(blast_file, dialect=csv.excel_tab)

    def find_by_phage(self, phage, cluster):
        for row in self.reader:
            # check if phage exists
            if phage in row[self.PHAGE_NAME] and cluster in row[self.CLUSTER]:
                # check if exists in genbank
                if 'True' in row[self.EXIST_IN_GENBANK]:
                    return row[self.ACCESSION]
        return -1
