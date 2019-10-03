#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import textwrap

def _main():
    parser = argparse.ArgumentParser(description="Concatenate FASTA")
    parser.add_argument('fasta', nargs='+', type=argparse.FileType('r'))
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), required=True)
    parser.add_argument('--name', '-n', required=True)
    options = parser.parse_args()

    seq = ''

    for one in options.fasta:
        for line in one:
            if line.startswith('>'): continue
            seq += line.strip()

    options.output.write('>' + options.name + '\n')
    for line in textwrap.wrap(seq, width=60):
        options.output.write(line + '\n')
    

if __name__ == '__main__':
    _main()
