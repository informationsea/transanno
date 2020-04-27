#!/usr/bin/env python3

import argparse
import gzip


def _main():
    parser = argparse.ArgumentParser("extract chromosome")
    parser.add_argument("fasta", help="Input Gzipped FASTA",
                        type=lambda x: gzip.open(x, "rt"))
    parser.add_argument("outfasta", help="Output FASTA",
                        type=argparse.FileType('w'))
    parser.add_argument("chromosome", help="chromosome name", nargs='+')
    options = parser.parse_args()

    chromosome_set = set(options.chromosome)

    in_target_chromosome = False
    for line in options.fasta:
        if line.startswith(">"):
            name = line[1:].split()[0]
            in_target_chromosome = name in chromosome_set

        if in_target_chromosome:
            options.outfasta.write(line)


if __name__ == "__main__":
    _main()
