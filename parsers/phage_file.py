import re

from csv_file import CSVWriter


class PhageDBReader:
    def __init__(self, infile, **kwargs):

        blast_file = open(infile, 'r')
        self.reader = blast_file.read()

    def write(self, outfile):
        csv_file = open(outfile, 'w')
        self.writer = CSVWriter(csv_file)

        # Header
        self.writer.write_row([
            'Query',
            'Name', 'Score', 'Expect',
            'QueryStart', 'QueryEnd',
            'SubjectStart', 'SubjectEnd'
        ])

        # Body
        for row in self.parse():
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
                r"(?P<Sequence>\w+)[^~]*(Length = )(?P<Length>\d*.\d*)(\s*)(Score = )(?P<Score>\d*.\d*)[^~]*(Expect = )(?P<Expect>\d*.\d*)[^~]*(Query: )(?P<QStart>\d*)[ acgt]*(?P<QEnd>\d*)[^~]*(Sbjct: )(?P<SubStart>\d*)[ acgt]*(?P<SubEnd>\d*)[^~]*",
                re.MULTILINE)

            iterator = pattern.finditer(sequence)

            for data in iterator:
                yield [
                    queryNumber,
                    data.group('Sequence'),
                    data.group('Length'),
                    data.group('Score'),
                    data.group('Expect'),
                    data.group('QStart'),
                    data.group('QEnd'),
                    data.group('SubStart'),
                    data.group('SubEnd')
                ]

