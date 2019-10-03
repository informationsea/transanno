#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
import collections
import itertools


def adaptive_open(mode: str):
    def _open(filename: str):
        if filename.endswith('.gz'):
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)
    return _open

def _main():
    parser = argparse.ArgumentParser(description="Check success rate")
    parser.add_argument('original', type=adaptive_open('rt'), help='(input) original VCF')
    parser.add_argument('--names', '-n', nargs='+', required=True, help='Names')
    parser.add_argument('--success', '-s', nargs='+', required=True, type=adaptive_open('rt'), help='(Input) Succeeded VCF')
    parser.add_argument('--failed', '-f', nargs='+', required=True, type=adaptive_open('rt'), help='(Input) Failed VCF')
    parser.add_argument('--output', '-o', required=True, type=adaptive_open('wt'))
    options = parser.parse_args()

    original_ids = set()
    for line in options.original:
        if line.startswith('#'): continue
        original_ids.add(line.split('\t')[2])

    success_ids = collections.defaultdict(set)
    failed_ids = collections.defaultdict(set)
    partial_success_ids = collections.defaultdict(set)

    writer = csv.writer(options.output)
    writer.writerow(['name', 'lifted', 'partial-lifted', 'unmapped'])

    for one_name, one_success, one_fail in zip(options.names, options.success, options.failed):
        print('loading', one_name)
        for line in one_success:
            if line.startswith('#'): continue
            success_ids[one_name].add(line.split('\t')[2])
        for line in one_fail:
            if line.startswith('#'): continue
            failed_ids[one_name].add(line.split('\t')[2])
        
        partial_success_ids[one_name] = success_ids[one_name] & failed_ids[one_name]
        failed_ids[one_name] -= partial_success_ids[one_name]
        success_ids[one_name] -= partial_success_ids[one_name]

        writer.writerow((one_name, len(success_ids[one_name]), len(partial_success_ids[one_name]), len(failed_ids[one_name])))

    writer.writerow([])
    writer.writerow([])
    writer.writerow(['names', 'lifted'])

    name_set = set(options.names)

    for num in range(1, len(options.names) + 1):
        for one in itertools.combinations(options.names, num):
            print(one)

            selected_names = set(one)
            
            common_set = set(success_ids[one[0]])
            for x in one[1:]:
                common_set &= success_ids[x]

            for x in name_set - selected_names:
                common_set -= success_ids[x]

            writer.writerow([','.join(one), len(common_set)])




if __name__ == '__main__':
    _main()
