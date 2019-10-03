#!/bin/bash -xe

./liftOver -multiple -minMatch=0.3 before-normalized-regions.bed ../lifttest/seq-a--to--seq-b.chain before-normalized-regions-lifted.bed before-normalized-regions-rejected.bed
python3 ../../chain-to-bed-vcf.py --track-name "seq-a lift test" --bed-output before-normalized.chain.bed --diff-bed-output before-normalized.chain.diff.bed --vcf-output before-normalized.chain.vcf --reference-fasta ../sequence/seq-a.fa --target-fasta ../sequence/seq-b.fa before-normalized.chain
bcftools sort -o before-normalized.chain.sorted.vcf.gz -O z before-normalized.chain.vcf
bcftools index -t before-normalized.chain.sorted.vcf.gz
bcftools norm -f ../sequence/seq-a.fa before-normalized.chain.vcf|bcftools sort -o before-normalized.chain.normed.sorted.vcf.gz -O z
bcftools index -t before-normalized.chain.normed.sorted.vcf.gz
bcftools norm -f ../sequence/seq-a.fa -O v -o after-normalized.vcf before-normalize.vcf