#!/usr/bin/env python3

import argparse
import csv
import json
import re
import collections

REGION_NAME = re.compile(r'CHUNK-([\w_]+)-(\d+)-(\d+)')


def _main():
    parser = argparse.ArgumentParser('create result json')
    parser.add_argument('mapped', type=argparse.FileType(
        'r'), help='Input mapped BED')
    parser.add_argument('unmapped', type=argparse.FileType(
        'r'), help='Input unmapped BED')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'),
                        required=True, help='Output JSON file')
    options = parser.parse_args()

    result = collections.defaultdict(list)

    reader = csv.reader(options.mapped, delimiter='\t', quotechar=None)
    for row in reader:
        matches = REGION_NAME.match(row[3])
        if not matches:
            continue
        result[(matches.group(1), int(matches.group(2)), int(matches.group(3)))].append(
            (row[0], int(row[1]), int(row[2]), row[5]))

    result_keys = list(result.keys())
    result_keys.sort()

    output_result = []
    for one in result_keys:
        output_result.append({
            "chrom": one[0],
            "start": one[1],
            "end": one[2],
            "mapped": [{"chrom": x[0], "start": x[1], "end": x[2], "strand": x[3]} for x in result[one]]
        })

    json.dump({"list": output_result}, options.output, indent="  ")


if __name__ == '__main__':
    _main()
