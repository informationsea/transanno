#!/usr/bin/env python3

import argparse
import gzip
import textwrap
import math

trans = str.maketrans("ATCG", "TAGC")


def _main():
    parser = argparse.ArgumentParser("Reverse complement of FASTA")
    parser.add_argument("fasta", help="Input FASTA",
                        type=argparse.FileType('r'))
    parser.add_argument("outfasta", help="Output FASTA",
                        type=argparse.FileType('w'))
    options = parser.parse_args()
    data = ""
    chromosome_name = ""

    for line in options.fasta:
        if line.startswith(">"):
            finish(options.outfasta, chromosome_name, data)
            chromosome_name = line[1:].strip()
            data = ""
        else:
            data += line.strip()
    finish(options.outfasta, chromosome_name, data)


def finish(output, chromosome_name, data):
    if chromosome_name:
        print(">" + chromosome_name + ":revcomp", file=output)
        rev_data = "".join(reversed(
            [x for x in data.translate(trans)]))
        for i in range(math.ceil(len(rev_data) / 60)):
            print(rev_data[60*i:60 * (i+1)], file=output)


if __name__ == "__main__":
    _main()
