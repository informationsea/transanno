#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import json


def _main():
    parser = argparse.ArgumentParser(description="Collect mapped regions with UCSC LiftOver and create json")
    parser.add_argument('original', type=argparse.FileType('r'))
    parser.add_argument('mapped', type=argparse.FileType('r'))
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    options = parser.parse_args()

    original_reader = csv.reader(options.original, delimiter='\t')
    mapped_reader = csv.reader(options.mapped, delimiter='\t')

    original = dict()

    for row in original_reader:
        original[row[3]] = {
            'original_chrom': row[0],
            'original_start': int(row[1]),
            'original_end': int(row[2]),
            'original_strand': row[5],
            'mapped': list()
        }

    for row in mapped_reader:
        original[row[3]]['mapped'].append({
            'chrom': row[0],
            'start': int(row[1]),
            'end': int(row[2]),
            'strand': row[5]
        })

    json.dump({'lift': list(original.values())}, options.output, indent='    ')
    

if __name__ == '__main__':
    _main()
