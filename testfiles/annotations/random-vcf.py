#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import random

def _main():
    parser = argparse.ArgumentParser(description="VCF random sampling")
    parser.add_argument('vcf', type=argparse.FileType('r'))
    parser.add_argument('--output', '-o', type=argparse.FileType('w'))
    options = parser.parse_args()

    if options.output:
        output = options.output
    else:
        output = sys.stdout

    for line in options.vcf:
        if line.startswith('#'):
            output.write(line)
        elif random.random() < 0.01:
            output.write(line)
    

if __name__ == '__main__':
    _main()
