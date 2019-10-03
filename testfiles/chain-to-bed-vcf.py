#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import argparse
import csv
import re
import sys
import pyfaidx

WHITE_SCAPE = re.compile(r'\s')
 
def _main():
    parser = argparse.ArgumentParser(description="Create BED and VCF file from chain file")
    parser.add_argument('chain', type=argparse.FileType('r'))
    parser.add_argument('--track-name', '-n', required=True)
    parser.add_argument('--bed-output', '-b', type=argparse.FileType('w'), required=True)
    parser.add_argument('--diff-bed-output', '-d', type=argparse.FileType('w'))
    parser.add_argument('--vcf-output', '-v', type=argparse.FileType('w'), required=True)
    parser.add_argument('--reference-fasta', '-r', required=True)
    parser.add_argument('--target-fasta', '-t', required=True)
    options = parser.parse_args()

    reference_fasta = pyfaidx.Fasta(options.reference_fasta)
    target_fasta = pyfaidx.Fasta(options.target_fasta)

    print('track name="{0} alignment" description="{0} alignment" visibility=2 itemRgb="On"'.format(options.track_name), file=options.bed_output)
    if options.diff_bed_output:
        print('track name="{0} difference" description="{0} difference" visibility=2 itemRgb="On"'.format(options.track_name), file=options.diff_bed_output)

    print('##fileformat=VCFv4.2', file=options.vcf_output)
    print('##INFO=<ID=TARGET_CHROM,Number=1,Type=String,Description="Target sequence chromosome">', file=options.vcf_output)
    print('##INFO=<ID=TARGET_POS,Number=1,Type=Integer,Description="Target sequence position">', file=options.vcf_output)
    print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">', file=options.vcf_output)
    print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">', file=options.vcf_output)
    print('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">', file=options.vcf_output)
    print('##ALT=<ID=DEL,Description="Deletion">', file=options.vcf_output)
    print('##ALT=<ID=INS,Description="Insertion of novel sequence">', file=options.vcf_output)
    print('##ALT=<ID=INDEL,Description="Deletion and Insertion of novel sequence">', file=options.vcf_output)

    for one in reference_fasta.keys():
        print('##contig=<ID={},length={}>'.format(one, len(reference_fasta[one])), file=options.vcf_output)
    
    print('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO', file=options.vcf_output)

    bed_writer = csv.writer(options.bed_output, delimiter='\t')
    if options.diff_bed_output:
        diff_bed_writer = csv.writer(options.diff_bed_output, delimiter='\t')
    vcf_output = csv.writer(options.vcf_output, delimiter='\t', lineterminator='\n')

    for line_num, line in enumerate(options.chain):
        if not line.strip():
            continue
        
        row = WHITE_SCAPE.split(line.strip())
        
        if row[0] == 'chain':
            current_chain = row

            if row[2] not in reference_fasta:
                reference_chrom = None
                continue

            reference_chrom = row[2]
            reference_start = int(row[5])
            reference_stop = int(row[6])
            
            target_chrom = row[7]
            target_length = int(row[8])
            target_start = int(row[10])
            target_stop = int(row[11])
            target_strand = row[9]
            
            if row[4] != '+':
                sys.exit('Invalid chain file at line: ' + line_num)

            if row[9] == '-':
                name = 'chain:-{}:{}-{}'.format(row[7], int(row[8]) - int(row[11]), int(row[8]) - int(row[10]))
            else:
                name = 'chain:+{}:{}-{}'.format(row[7], row[10], row[11])

            bed_writer.writerow([row[2], row[5], row[6], name, 0, row[9], row[5], row[6], '0,104,183'])
        elif reference_chrom:
            #print(row[0])
            reference_current = reference_start + int(row[0])
            target_current = target_start + int(row[0])

            if options.diff_bed_output:
                reference_seq = str(reference_fasta[reference_chrom][reference_start:reference_current]).upper()
                target_seq = str(get_target_seq(target_fasta, target_chrom, target_length, target_strand, target_start, target_current)).upper()
    
                if reference_seq != target_seq:
                    for i, (rb, tb) in enumerate(zip(reference_seq, target_seq)):
                        if rb != tb:
                            target_this = target_start + i
                            target_pos = get_target_pos(target_length, target_strand, target_this)
                            
                            name = 'SNV:{}{}:{}:{}>{}'.format(target_strand, target_chrom, target_pos, rb, tb)
                            diff_bed_writer.writerow([reference_chrom, reference_start + i, reference_start + i + 1, name, 0, target_strand, reference_start + i, reference_start + i + 1, '238,120,0'])
                            vcf_output.writerow([reference_chrom, reference_start + i + 1, '.', rb, tb, '.', '.', 'TARGET_CHROM={};TARGET_POS={}'.format(target_chrom, target_pos+1)])
                        
                    #print('diff seq at line: {}'.format(line_num))
                    #print(reference_seq[:10], target_seq[:10])
                    #print(reference_seq[-10:], target_seq[-10:])
                    #sys.exit('Invalid chain file at line: {}: different sequence'.format(line_num))

            reference_start = reference_current
            target_start = target_current

            if len(row) >= 3:
                reference_current = reference_start + int(row[1])
                target_current = target_start + int(row[2])

                start_offset = 0
                if row[1] == '0' or row[2] == '0':
                    start_offset = -1
                
                if row[1] == '0':
                    reference_seq = str(reference_fasta[reference_chrom][(reference_start + start_offset):reference_start]).upper()
                else:
                    reference_seq = str(reference_fasta[reference_chrom][(reference_start + start_offset):reference_current]).upper()

                if row[2] == '0':
                    target_seq = str(get_target_seq(target_fasta, target_chrom, target_length, target_strand, target_start + start_offset, target_start)).upper()
                else:
                    target_seq = str(get_target_seq(target_fasta, target_chrom, target_length, target_strand, target_start + start_offset, target_current)).upper()

                sv_type = None
                sv_len = None
                sv_ref = None
                if len(reference_seq) > 100 or len(target_seq) > 100:
                    sv_len = len(target_seq) - len(reference_seq)
                    sv_ref = reference_seq[0]
                    if len(reference_seq) == 1:
                        sv_type = 'INS'
                    elif len(target_seq) == 1:
                        sv_type = 'DEL'
                    else:
                        sv_type = 'INDEL'

                if len(reference_seq) > 100:
                    reference_seq = '[LEN:{}]'.format(len(reference_seq))
                if len(target_seq) > 100:
                    target_seq = '[LEN:{}]'.format(len(target_seq))

                target_pos = get_target_region(target_length, target_strand, target_start, target_current)

                
                if row[1] == '1' and row[2] == '1':
                    color = '238,120,0'
                    name = 'SNV:{}{}:{}:{}>{}'.format(target_strand, target_chrom, target_pos[0] + start_offset, reference_seq, target_seq)
                elif row[1] == '0':
                    color = '199,0,103'
                    name = 'INS:{}{}:{}-{}:{}>{}'.format(target_strand, target_chrom, target_pos[0] + start_offset, target_pos[1], reference_seq, target_seq)
                elif row[2] == '0':
                    color = '215,0,53'
                    name = 'DEL:{}{}:{}-{}:{}>{}'.format(target_strand, target_chrom, target_pos[0] + start_offset, target_pos[1], reference_seq, target_seq)
                else:
                    color = '234,85,80'
                    name = 'INDEL:{}{}:{}-{}:{}>{}'.format(target_strand, target_chrom, target_pos[0] + start_offset, target_pos[1], reference_seq, target_seq)
                    
                if options.diff_bed_output:
                    diff_bed_writer.writerow([reference_chrom, reference_start + start_offset, max(reference_current, reference_start + 1), name, 0, target_strand, reference_start, max(reference_current, reference_start + 1), color])

                if sv_type:
                    vcf_output.writerow([reference_chrom, reference_start + start_offset + 1, '.', sv_ref, '<' + sv_type + '>', '.', '.', 'TARGET_CHROM={};TARGET_POS={};END={};SVTYPE={};SVLEN={}'.format(target_chrom, target_start + start_offset, reference_current, sv_type, sv_len)])
                else:
                    if reference_seq != target_seq:
                        vcf_output.writerow([reference_chrom, reference_start + start_offset + 1, '.', reference_seq, target_seq, '.', '.', 'TARGET_CHROM={};TARGET_POS={}'.format(target_chrom, target_start + start_offset)])
                
                reference_start = reference_current
                target_start = target_current
            

def get_target_pos(target_length, target_strand, pos):
    if target_strand == '-':
        return target_length - pos
    return pos

def get_target_region(target_length, target_strand, pos_start, pos_end):
    if target_strand == '-':
        return (target_length - pos_end, target_length - pos_start)
    else:
        return (pos_start, pos_end)

def get_target_seq(target_fasta, target_chrom, target_length, target_strand, pos_start, pos_end):
    pos_start, pos_end = get_target_region(target_length, target_strand, pos_start, pos_end)
    if target_strand == '-':
        return target_fasta[target_chrom][pos_start:pos_end].reverse.complement
    return target_fasta[target_chrom][pos_start:pos_end]
 
if __name__ == '__main__':
    _main() 
