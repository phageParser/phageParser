import csv
import cStringIO


class CSVReader:
    """
    Adapted from https://docs.python.org/2/library/csv.html
    A CSV reader which will iterate over lines in the CSV file "f"
    """

    def __init__(self, f, dialect=csv.excel, **kwargs):
        self.reader = csv.reader(f, dialect=dialect, **kwargs)

    def next(self):
        row = self.reader.next()
        return [s for s in row]

    def __iter__(self):
        return self


class CSVWriter:
    """
    Adapted from https://docs.python.org/2/library/csv.html
    A CSV writer which will write rows to CSV file "f"
    """

    def __init__(self, file, dialect=csv.excel, **kwargs):
        # Redirect output to a queue
        self.queue = cStringIO.StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, delimiter=",", skipinitialspace=False, **kwargs)
        self.stream = file

    def write_row(self, row):
        self.writer.writerow([s for s in row])

        # Fetch UTF-8 output from the queue ...
        data = self.queue.getvalue()

        # write to the target stream
        self.stream.write(data)

        # empty queue
        self.queue.truncate(0)
