#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import argparse
import csv
import re

CIGER = re.compile(r'(\d+)([A-Z])')
 
def _main():
    parser = argparse.ArgumentParser(description="minimap to liftover chain")
    parser.add_argument('paf', type=argparse.FileType('r'), default='./data/01-draft.hg19.cs.paf', nargs='?', help='default:%(default)s')
    #parser.add_argument('--output-chain-to-reference', type=argparse.FileType('w'), default="output/minimap2-chain-01draftToHg19.chain", help='default:%(default)s')
    parser.add_argument('--output-chain-to-target', type=argparse.FileType('w'), default="output/minimap2-chain-hg19To01draft.chain", help='default:%(default)s')
    options = parser.parse_args()

    csv.field_size_limit(1000000000)

    reader = csv.reader(options.paf, delimiter='\t')
    #writer1 = csv.writer(options.output_chain_to_reference, delimiter='\t')
    writer2 = csv.writer(options.output_chain_to_target, delimiter='\t')

    count = 0
    
    for row in reader:
        count += 1
        #if row[4] != '+': continue
        if row[4] == '+':
            writer2.writerow(['chain', '4900', row[0], row[1], '+', row[2], row[3], row[5], row[6], row[4], row[7], row[8], count])
            #writer1.writerow(['chain', '4900', row[5], row[6], row[4], row[7], row[8], row[0], row[1], '+', row[2], row[3], count])
        else:
            seqlen = int(row[6])
            writer2.writerow(['chain', '4900', row[0], row[1], '+', row[2], row[3], row[5], row[6], row[4], seqlen - int(row[8]), seqlen - int(row[7]), count])
            ##writer1.writerow(['chain', '4900', row[5], row[6], row[4], seqlen - int(row[8]), seqlen - int(row[7]), row[0], row[1], '+', row[2], row[3], count])
        
        ciger = [x for x in row if x.startswith('cg:Z:')][0][5:]
        pos = 0

        #print(row[:9])

        ciger_list = []

        while True:
            m = CIGER.search(ciger, pos)
            if m == None:
                break
            #print(m.groups())
            ciger_list.append(m.groups())
            if m.pos != pos:
                print('ERROR', m.pos)
                break
            pos = m.end()

        if row[4] == '-':
            ciger_list.reverse()

        #print(int(row[2]) + int(ciger_list[0][0]), int(row[7]) + int(ciger_list[0][0]))

        chain_group = []
        current_group = []
        for i, one in enumerate(ciger_list):
            if one[1] == 'M' and current_group:
                chain_group.append(current_group)
                current_group = list()
            current_group.append(one)
        #chain_group.append(one)
        last = one

        for one in chain_group:
            matched = int(one[0][0])
            deletion = sum([int(x[0]) for x in one if x[1] == 'D'])
            insertion = sum([int(x[0]) for x in one if x[1] == 'I'])
            writer2.writerow([matched, insertion, deletion])
        writer2.writerow([last[0]])
        writer2.writerow([])

#        if row[4] == '+':
#            for one in chain_group:
#                matched = int(one[0][0])
#                deletion = sum([int(x[0]) for x in one if x[1] == 'D'])
#                insertion = sum([int(x[0]) for x in one if x[1] == 'I'])
#                writer1.writerow([matched, deletion, insertion])
#            writer1.writerow([last[0]])
#            writer1.writerow([])
#        else:
#            # TODO fill here
#            pass
##            for i in range(len(chain_group)):
##                j = len(chain_group) - i - 1
##                matched = last[0] if i == 0 else chain_group[j - 1]
##                deletion = sum([int(x[0]) for x in chain_group[j] if x[1] == 'D'])
##                insertion = sum([int(x[0]) for x in chain_group[j] if x[1] == 'I'])
##                writer1.writerow([matched, deletion, insertion])
##            writer1.writerow([last[0]])

        #print(chain_group)
        #count += 1
        #if count > 1: break

 
if __name__ == '__main__':
    _main() 
