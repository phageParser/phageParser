#!/usr/bin/env python
import time
from parsers.phage_file import PhageDBReader


if __name__ == '__main__':
    current_milli_time = lambda: int(round(time.time() * 1000))
    parser = PhageDBReader('data/blast-phagesdb.txt')
    # parser.write('output/blast-phagesdb.%s.csv' % current_milli_time())

    parser.filterByExpect('output/blast-phagesdb.%s.csv' % current_milli_time(), 0.21)