#!/usr/bin/env python3

import argparse
import csv


def _main():
    parser = argparse.ArgumentParser('create test region from fai')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Input .fai')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'),
                        required=True, help='Output BED file')
    parser.add_argument('--chunk-size', type=int, default=1000,
                        help='Chunk size (default: %(default)s)')
    parser.add_argument('--step-size', type=int, default=500,
                        help='Chunk size (default: %(default)s)')
    options = parser.parse_args()

    reader = csv.reader(options.fai, delimiter='\t', quotechar=None)
    for row in reader:
        current_start = 0
        while current_start < int(row[1]):
            options.output.write(
                f"{row[0]}\t{current_start}\t{current_start + options.chunk_size}\tCHUNK-{row[0]}-{current_start}-{current_start + options.chunk_size}\t0\t+\n")
            current_start += options.step_size
        options.output.write(
            f"{row[0]}\t{int(row[1]) - options.chunk_size}\t{int(row[1])}\tCHUNK-{row[0]}-{int(row[1]) - options.chunk_size}-{int(row[1])}\t0\t+\n")


if __name__ == '__main__':
    _main()
