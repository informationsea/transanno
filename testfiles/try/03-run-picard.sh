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

check_command java
check_command curl

GENOMES=(AK1v2 HX1.1 JG1 CHM1_1.1 hg38)
VCFS=(common_all_20180423.indels common_all_20180423.snps common_all_20180423 clinvar_20190923 clinvar_20190923.indels clinvar_20190923.snps)

if [ ! -e inputs/picard.jar ]; then
    curl -o inputs/picard.jar -L https://github.com/broadinstitute/picard/releases/download/2.20.8/picard.jar
fi

for one in ${GENOMES[@]} hg19; do
    if [ ! -e inputs/${one}.dict ]; then
        java -Dpicard.useLegacyParser=false -jar inputs/picard.jar CreateSequenceDictionary -R inputs/${one}.fa -O inputs/${one}.dict
    fi
done

for vcf in ${VCFS[@]}; do
    mkdir -p inputs/vcf_with_chr
    if [ ! -e inputs/vcf_with_chr/${vcf}.vcf.gz ];then
        gzip -dc inputs/${vcf}.vcf.gz|sed -e 's/^\([1-9XYM]\)/chr\1/'|sed -e 's/##contig=<ID=/##contig=<ID=chr/'|gzip -c > inputs/vcf_with_chr/${vcf}.vcf.gz &
    fi
done

wait

mkdir -p output/picard-vcf

for one in ${GENOMES[@]}; do
    for vcf in ${VCFS[@]}; do
        java -Xmx5G -Xms5G -Dpicard.useLegacyParser=false -jar inputs/picard.jar LiftoverVcf -I inputs/vcf_with_chr/${vcf}.vcf.gz -O output/picard-vcf/picard-${one}-${vcf}.success.vcf.gz -CHAIN output/hg19-to-${one}.chain.gz -REJECT output/picard-vcf/picard-${one}-${vcf}.fail.vcf.gz -R inputs/${one}.fa --WRITE_ORIGINAL_POSITION true --WRITE_ORIGINAL_ALLELES true --RECOVER_SWAPPED_REF_ALT true --ALLOW_MISSING_FIELDS_IN_HEADER true --DISABLE_SORT true > output/picard-vcf/stdout-${one}-${vcf}.txt 2> output/picard-vcf/stderr-${one}-${vcf}.txt &
    done
    wait
done
