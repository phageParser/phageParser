import re

from csv_file import CSVWriter


class PhageDBReader:
    EXPECT_COLUMN = 4  # zero-based

    def __init__(self, infile, **kwargs):

        blast_file = open(infile, 'r')
        self.reader = blast_file.read()

    def write(self, outfile):
        csv_file = open(outfile, 'w')
        self.writer = CSVWriter(csv_file)

        # Header
        self.writer.write_row([
            'Query',
            'Name', 'Length', 'Score', 'Expect',
            'QueryStart', 'QueryEnd',
            'SubjectStart', 'SubjectEnd'
        ])

        # Body
        for row in self.parse():
            self.writer.write_row(row)

    def filterByExpect(self, outfile, threshold=1.0):
        csv_file = open(outfile, 'w')
        self.writer = CSVWriter(csv_file)

        # Header
        self.writer.write_row([
            'Query',
            'Name', 'Length', 'Score', 'Expect',
            'QueryStart', 'QueryEnd',
            'SubjectStart', 'SubjectEnd'
        ])

        # Body
        for row in self.parse():
            if row[self.EXPECT_COLUMN] <= threshold:
                self.writer.write_row(row)

    def queries(self):
        # split file up into queries
        queries = re.split(r'Query= ', self.reader)

        # queries[0] will be the header and not needed
        queries.pop(0)

        for query in queries:
            yield query

    def sequences(self):
        for query in self.queries():
            queryNumber = float(query.splitlines()[0])

            # split query into matches
            matches = re.split(r'>', query)
            matches.pop(0)

            for match in matches:
                yield (queryNumber, match)

    def parse(self):
        for (queryNumber, sequence) in self.sequences():
            pattern = re.compile(
                r"(?P<Sequence>\w+)[^~]*(Length=)(?P<Length>\d*.\d*)(\s*)(Score = )(?P<Score>\d*.\d*)[^~]*(Expect = )(?P<Expect>\d*.\d*-*\d*)[^~]*(Query\s*)(?P<QStart>\d*)\s*[ACGT-]*\s*(?P<QEnd>\d*)[^~]*(Sbjct\s*)(?P<SubStart>\d*)\s*[ACGT-]*\s*(?P<SubEnd>\d*)[^~]*",
                re.MULTILINE)

            iterator = pattern.finditer(sequence)

            for data in iterator:
                yield [
                    float(queryNumber),
                    data.group('Sequence'),
                    int(data.group('Length')),
                    float(data.group('Score')),
                    float(data.group('Expect')),
                    int(data.group('QStart')),
                    int(data.group('QEnd')),
                    int(data.group('SubStart')),
                    int(data.group('SubEnd'))
                ]
