#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pyliftover
import csv
import math
import json


def _main():
    parser = argparse.ArgumentParser(description="Create correct lift position with pylift over")
    parser.add_argument('chain')
    parser.add_argument('targetchrom')
    parser.add_argument('chromlen', type=int)
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    options = parser.parse_args()

    STEP = 3000

    lo = pyliftover.LiftOver(options.chain)

    expected_results = []
    
    for i in range(0, math.floor(options.chromlen / STEP)):
        pos = i * STEP
        mapped = lo.convert_coordinate(options.targetchrom, pos)
        expected_results.append({
            'original_chrom': options.targetchrom,
            'original_pos': pos,
            'mapped': [
                {'chrom': x[0], 'pos': x[1], 'strand': x[2]} for x in mapped
            ]
        })
    json.dump({'list': expected_results}, options.output)

    

if __name__ == '__main__':
    _main()
