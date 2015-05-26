#!/usr/bin/env python
import os
import time
import argparse

import logging
from logging.config import fileConfig

from parsers.phage_file import PhageDBReader

fileConfig('logging.ini')
logger = logging.getLogger()

def parse_args():
    parser = argparse.ArgumentParser(description='Cleans up the results returned from phagesdb.org')
    parser.add_argument('--infile', help='the input file returned from phagesdb.org', required=True)
    parser.add_argument('--output', help='the output directory (default = output/)', required=True, default='output/')
    parser.add_argument('--threshold', help='the threshold (default = 1.0)', required=False, type=float, default=1.0)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()

    infile = args.infile
    output_dir = args.output
    outfilename = os.path.basename(os.path.splitext(infile)[0])
    threshold = args.threshold

    current_milli_time = lambda: int(round(time.time() * 1000))
    parser = PhageDBReader(infile)
    # parser.write('output/blast-phagesdb.%s.csv' % current_milli_time())

    outfile = 'output/%s.%s.csv' % (outfilename, current_milli_time())
    logger.info("Running phagesdb filter, infile=[%s], outputfile=[%s], threshold=[%.2f]" % (infile, output_dir, threshold))
    parser.filterByExpect(outfile, 0.21)