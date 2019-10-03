#!/usr/bin/env bash
#$ -cwd -N run-transanno -pe def_slot 20

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

GENOMES=(AK1v2 HX1.1 JG1 CHM1_1.1 hg38)
VCFS=(common_all_20180423.indels common_all_20180423.snps common_all_20180423 clinvar_20190923 clinvar_20190923.indels clinvar_20190923.snps)

check_command cargo
check_command bcftools

#pushd ../..
#cargo build --release
#popd

cp ../../target/release/transanno .

mkdir -p output
./transanno --version > output/transanno-version.txt

for one in ${GENOMES[@]}; do
    ./transanno minimap2-to-chain ./inputs/hg19-to-${one}.paf -o output/hg19-to-${one}.chain.gz &
done
wait

mkdir -p output/vcf

for one in ${GENOMES[@]}; do
    for vcf in ${VCFS[@]}; do
        ./transanno liftvcf --chain output/hg19-to-${one}.chain.gz --query inputs/${one}.fa --reference inputs/hg19.fa --vcf inputs/${vcf}.vcf.gz --output output/vcf/liftvcf-${one}-${vcf}-success.vcf.gz --fail output/vcf/liftvcf-${one}-${vcf}-fail.vcf.gz 2> output/vcf/stderr-${one}-${vcf}.txt &
    done
done

wait

mkdir -p output/sorted-vcf

for one in ${GENOMES[@]}; do
    for vcf in ${VCFS[@]}; do
        bcftools sort -O z -o output/sorted-vcf/liftvcf-${one}-${vcf}-success.vcf.gz output/vcf/liftvcf-${one}-${vcf}-success.vcf.gz &
    done
done

wait

for one in ${GENOMES[@]}; do
    for vcf in ${VCFS[@]}; do
        bcftools view -O z -o output/sorted-vcf/liftvcf-${one}-${vcf}-fail.vcf.gz output/vcf/liftvcf-${one}-${vcf}-fail.vcf.gz &
    done
done

wait

for one in ${GENOMES[@]}; do
    for vcf in ${VCFS[@]}; do
        bcftools index -t output/sorted-vcf/liftvcf-${one}-${vcf}-success.vcf.gz &
    done
done

wait

for one in ${GENOMES[@]}; do
    for vcf in ${VCFS[@]}; do
        bcftools index -t output/sorted-vcf/liftvcf-${one}-${vcf}-fail.vcf.gz &
    done
done

wait
