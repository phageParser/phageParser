import csv


class GenomeExtractor:
    """
    @author:         Nigesh Shakya(neaGaze)
    @timestamp:      06/08/2016 @ 12:54am (UTC)

    The Class file that examines the file downloaded via prefetch (see #2),
    and returns a string corresponding to the genome starting n bases
    before startPosition and carrying on to m bases after endPosition. If
    either n or m are not provided, default them to 100.
    This patch pertains to issue #3 in the github wiki
    (https://github.com/goyalsid/phageParser/issues/3)
    """

    def __init__(self, infile, **kwargs):
        blast_file = open(infile, 'r')
        self.reader = csv.reader(blast_file, dialect=csv.excel_tab)

    def find_neighbours(self, start_position, end_position, n=100, m=100):
        """ Genome Extractor """
        has_data_found = False
        matched_gene = ""
        for row in self.reader:

            if row[0] == "//":
                break

            if row[0].startswith("ORIGIN"):
                has_data_found = True
                continue

            if has_data_found:
                row_vals = row[0].split()
                pos = int(row_vals[0])
                if start_position - n <= pos < end_position + m:
                    matched_gene += " ".join(row_vals)
                    matched_gene += "\n"

        return matched_gene
