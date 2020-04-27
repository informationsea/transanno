#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import math


def _main():
    parser = argparse.ArgumentParser(description="Create liftOver test regions")
    parser.add_argument('chromosome')
    parser.add_argument('chromlen', help='chromosome length', type=int)
    parser.add_argument('--step', '-s', type=int, default=3000, help='step size: (default: %(default)s)')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), required=True)
    options = parser.parse_args()

    for i in range(math.floor(options.chromlen / options.step)):
        start = i * options.step
        end = (i + 1) * options.step
        print('\t'.join([str(x) for x in [options.chromosome, start, end, 'region-{}'.format(i), '.', '+']]), file=options.output)
    

if __name__ == '__main__':
    _main()
