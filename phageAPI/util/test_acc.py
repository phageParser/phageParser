from unittest import TestCase
import io

from util.acc import read_accession_file


def test_read_accession_file_empty():
    assert list(read_accession_file(io.StringIO(""))) == []


def test_read_accession_file_check_formatting():
    f = io.StringIO("""
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
