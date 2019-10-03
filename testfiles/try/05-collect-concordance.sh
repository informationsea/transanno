#!/usr/bin/env bash

set -xe

if [ "$JOB_ID" != "" ]; then
    cd $(dirname $0)
fi


function check_command() {
    if !(type "$1" > /dev/null 2>&1); then
        echo "$1 is required to run this script"
        exit 1
    fi
}

check_command python3

GENOMES=(AK1v2 HX1.1 CHM1_1.1 hg38)
VCFS=(common_all_20180423.indels common_all_20180423.snps common_all_20180423)

mkdir -p output/crossmap

for one_genome in ${GENOMES[@]}; do
    for one_vcf in ${VCFS[@]}; do
        python3 scripts/check-common.py inputs/${one_vcf}.vcf.gz --names transanno picard crossmap --success output/vcf/liftvcf-${one_genome}-${one_vcf}-success.vcf.gz output/picard-vcf/picard-${one_genome}-${one_vcf}.success.vcf.gz output/crossmap/crossmap-${one_genome}-${one_vcf}.vcf.gz --fail output/vcf/liftvcf-${one_genome}-${one_vcf}-fail.vcf.gz output/picard-vcf/picard-${one_genome}-${one_vcf}.fail.vcf.gz output/crossmap/crossmap-${one_genome}-${one_vcf}.vcf.unmap -o output/compare/common-lifted-${one_genome}-${one_vcf}.txt &
    done
    wait
done

python3 scripts/check-concordance.py inputs/common_all_20180423.indels.vcf.gz inputs/common_all_20180418_GRCh38.indels.vcf.gz output/vcf/liftvcf-hg38-common_all_20180423.indels-success.vcf.gz output/picard-vcf/picard-hg38-common_all_20180423.indels.success.vcf.gz -o output/compare/compare-hg38-common_all_20180423.indels.txt -f output/compare/compare-hg38-common_all_20180423.indels.failed.vcf -c output/compare/compare-hg38-common_all_20180423.indels.combination.csv

python3 scripts/check-concordance.py inputs/common_all_20180423.snps.vcf.gz inputs/common_all_20180418_GRCh38.snps.vcf.gz output/vcf/liftvcf-hg38-common_all_20180423.snps-success.vcf.gz output/picard-vcf/picard-hg38-common_all_20180423.snps.success.vcf.gz -o output/compare/compare-hg38-common_all_20180423.snps.txt -f output/compare/compare-hg38-common_all_20180423.snps.failed.vcf -c output/compare/compare-hg38-common_all_20180423.snps.combination.csv

python3 scripts/check-concordance.py inputs/common_all_20180423.vcf.gz inputs/common_all_20180418_GRCh38.vcf.gz output/vcf/liftvcf-hg38-common_all_20180423-success.vcf.gz output/picard-vcf/picard-hg38-common_all_20180423.success.vcf.gz -o output/compare/compare-hg38-common_all_20180423.txt -f output/compare/compare-hg38-common_all_20180423.failed.vcf -c output/compare/compare-hg38-common_all_20180423.combination.csv
