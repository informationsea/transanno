#!/usr/bin/env bash
#$ -S /bin/bash -cwd -N run-picard -pe def_slot 20

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

check_command singularity
check_command bcftools
check_command bgzip

GENOMES=(AK1v2 HX1.1 JG1 CHM1_1.1 hg38)
VCFS=(common_all_20180423.indels common_all_20180423.snps common_all_20180423)

if [ ! -e inputs/crossmap.simg ];then
    singularity build inputs/crossmap.simg docker://informationsea/crossmap
fi

for vcf in ${VCFS[@]}; do
    mkdir -p inputs/vcf_mono_allelic
    if [ ! -e inputs/vcf_mono_allelic/${vcf}.vcf.gz ];then
        bcftools norm -m -any -N inputs/${vcf}.vcf.gz|sed -e 's/^\([1-9XYM]\)/chr\1/'|sed -e 's/##contig=<ID=/##contig=<ID=chr/'|gzip -c > inputs/vcf_mono_allelic/${vcf}.vcf.gz &
    fi
done

wait

for one in ${GENOMES[@]}; do
    mkdir -p output/crossmap
    for vcf in ${VCFS[@]}; do
        singularity exec --no-home inputs/crossmap.simg CrossMap.py vcf output/hg19-to-${one}.chain.gz inputs/vcf_mono_allelic/${vcf}.vcf.gz inputs/${one}.fa output/crossmap/crossmap-${one}-${vcf}.vcf > output/crossmap/crossmap-stdout-${one}-${vcf}.txt 2> output/crossmap/crossmap-stderr-${one}-${vcf}.txt && bgzip output/crossmap/crossmap-${one}-${vcf}.vcf &
    done
done

wait
