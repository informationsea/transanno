#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import gzip
import collections
import os.path
import csv


def adaptive_open(path, mode = 'rt', **param):
    if path.endswith('.gz'):
        return gzip.open(path, mode, **param)
    return open(path, mode, **param)


def _main():
    parser = argparse.ArgumentParser(description="Check lift over result concordance (variant ID is required)")
    parser.add_argument('original', help='(Input) original VCF')
    parser.add_argument('rightdata', help='(Input) Correct lifted VCF')
    parser.add_argument('inputs', nargs='+', help='(Input) lifted VCF')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), required=True, help='(Output) result CSV')
    parser.add_argument('--failed-vcf', '-f', type=argparse.FileType('w'), required=True, help='(Output) failed VCF')
    parser.add_argument('--combination', '-c', type=argparse.FileType('w'), required=True, help='(Output) result CSV')
    options = parser.parse_args()

    original = dict()
    with adaptive_open(options.original, 'rt', errors='ignore') as f:
        for line in f:
            if line.startswith('#'): continue
            if line.startswith('chrMT'): continue
            if line.startswith('MT'): continue
            if line.startswith('NW_'): continue
            elements = line.split('\t')
            original[elements[2]] = (norm_chrom(elements[0]), elements[1], elements[3], set(elements[4].split(',')))

    print('Total:', len(original))

    rightdata = dict()
    with adaptive_open(options.rightdata, 'rt', errors='ignore') as f:
        for line in f:
            if line.startswith('#'): continue
            if line.startswith('chrMT'): continue
            if line.startswith('NM_'): continue
            elements = line.split('\t')
            rightdata[elements[2]] = (norm_chrom(elements[0]), elements[1], elements[3], set(elements[4].split(',')))
            

    writer = csv.DictWriter(options.output, fieldnames=['name', 'lifted', 'lifted-correct', 'lifted-incorrect', 'lifted-unknown', 'failed', 'failed-incorrect', 'failed-unknown'])
    writer.writeheader()

    incorrect_variant_dict = collections.defaultdict(set)
    failed_variants = collections.defaultdict(set)

    for one_input in options.inputs:
        with adaptive_open(one_input, 'rt', errors='ignore') as f:
            name = os.path.basename(one_input)
            print('start loading', name)

            processed_ids = set()
            success_variants = 0
            correct_variants = 0
            incorrect_variants = 0
            unknown_variants = 0

            for line in f:
                if line.startswith('#'): continue
                if line.startswith('chrMT'): continue
                if line.startswith('NM_'): continue
                elements = line.split('\t')
                variant_id_tmp = elements[2]
                one_variant = (norm_chrom(elements[0]), elements[1], elements[3], set(elements[4].split(',')))

                for variant_id in variant_id_tmp.split(';'):
                    if variant_id in processed_ids:
                        print('Duplicated ID:', variant_id)
                        continue
                    processed_ids.add(variant_id)
                    
                    success_variants += 1
    
                    if variant_id not in rightdata:
                        unknown_variants += 1
                    elif rightdata[variant_id] == one_variant:
                        correct_variants += 1
                    else:
                        incorrect_variants += 1
                        incorrect_variant_dict[name].add(variant_id)

            failed_variants[name] = set(original.keys()) - processed_ids

            failed_incorrect = failed_variants[name] & rightdata.keys()

            writer.writerow({
                'name': name,
                'lifted': success_variants,
                'lifted-correct': correct_variants,
                'lifted-incorrect': incorrect_variants,
                'lifted-unknown': unknown_variants,
                'failed': len(failed_variants[name]),
                'failed-incorrect': len(failed_incorrect),
                'failed-unknown': len(failed_variants[name]) - len(failed_incorrect)
            })

            print('Lifted variants', name, success_variants)
            print('          Correct', name, correct_variants)
            print('        Incorrect', name, incorrect_variants)
            print('          Unknown', name, unknown_variants)
            print('Failed variants', name, len(original) - success_variants)


    all_failed_incorrect_variants = set()
    for k, v in failed_variants.items():
        all_failed_incorrect_variants |= v
    for v in incorrect_variant_dict.values():
        all_failed_incorrect_variants |= v
        
    all_failed_incorrect_variant_list = [(x, original[x]) for x in all_failed_incorrect_variants]
    all_failed_incorrect_variant_list.sort(key=lambda x: (chrom_priority(x[1][0]), int(x[1][1])))

    all_chrom = set([x[1][0] for x in all_failed_incorrect_variant_list])
    all_chrom_list = list(all_chrom)
    all_chrom_list.sort(key=lambda x: chrom_priority(x))

    print('##fileformat=VCFv4.1', file=options.failed_vcf)
    print('##INFO=<ID=FAILED,Number=.,Type=String,Description="Failed">', file=options.failed_vcf)
    print('##INFO=<ID=INCORRECT,Number=.,Type=String,Description="Incorrect">', file=options.failed_vcf)
    for one in all_chrom_list:
        print('##contig=<ID={}>'.format(one), file=options.failed_vcf)
    print('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']), file=options.failed_vcf)

    fail_pattern = collections.defaultdict(int)
    
    for one in all_failed_incorrect_variant_list:
        incorrect = [x[0] for x in incorrect_variant_dict.items() if one[0] in x[1]]
        failed = [x[0] for x in failed_variants.items() if one[0] in x[1]]

        info = ''
        if failed:
            info += 'FAILED=' + ','.join(failed)
        if incorrect:
            if info:
                info += ';'
            info += 'INCORRECT=' + ','.join(incorrect)

        fail_pattern[(frozenset(failed), frozenset(incorrect))] += 1
        
        print('\t'.join([one[1][0], one[1][1], one[0], one[1][2], ','.join(one[1][3]), '.', '.', info]), file=options.failed_vcf)

    combination = csv.DictWriter(options.combination, fieldnames=['failed', 'incorrect', 'count'])
    combination.writeheader()
        
    for k, v in fail_pattern.items():
        combination.writerow({'failed': ','.join(k[0]), 'incorrect': ','.join(k[1]), 'count': v})
                    

def norm_chrom(name):
    if name.startswith('chr'):
        return name
    return 'chr' + name


def chrom_priority(name):
    if name == 'chrX':
        return 23
    if name == 'chrY':
        return 24
    return int(name[3:])
    
    

if __name__ == '__main__':
    _main()

