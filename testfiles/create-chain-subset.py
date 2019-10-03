#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import sys


def _main():
    parser = argparse.ArgumentParser(description="Create subset chain file")
    parser.add_argument("chain", type=argparse.FileType('r'))
    parser.add_argument(
        "-o", "--output", type=argparse.FileType('w'), required=True)
    parser.add_argument("-q", "--query", nargs="+",
                        help="Query sequence names")
    parser.add_argument("-r", "--reference", nargs="+",
                        help="Reference sequence names")
    options = parser.parse_args()

    status = False
    query_seqnames = set(options.query)
    reference_seqnames = set(options.reference)
    writer = csv.writer(options.output, delimiter=' ', quotechar=None)

    for row in csv.reader(options.chain, delimiter=' ', quotechar=None):
        if len(row) >= 1 and row[0] == 'chain':
            # print(row)
            if row[2] in reference_seqnames and row[7] in query_seqnames:
                status = True
                writer.writerow(row)
            else:
                status = False
        elif status:
            writer.writerow(row)


if __name__ == '__main__':
    _main()
