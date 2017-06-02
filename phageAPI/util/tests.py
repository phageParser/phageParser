from unittest import TestCase
import StringIO

from util.acc import read_accession_file


def test_read_accession_file_empty():
    assert list(read_accession_file(StringIO.StringIO(""))) == []


def test_read_accession_file_check_formatting():
    f = StringIO.StringIO("""
# comment on left
123left
                        345right
        987middle              
        # comment in middle, with trailing space

        # ^^ blank line
""")
    assert list(read_accession_file(f)) == [
        '123left',
        '345right',
        '987middle',
    ]
