#!/usr/bin/env bash

set -xe

cd $(dirname $0)

mkdir -p inputs

function check_command() {
    if !(type "$1" > /dev/null 2>&1); then
        echo "$1 is required to run this script"
        exit 1
    fi
}

check_command minimap2
check_command samtools
check_command bcftools
check_command gzip
check_command curl
check_command cargo
check_command python3

function download_gzip() {
    URL=$2
    FILE=$1
    if [ ! -e $FILE ]; then
        curl -o ${FILE}.gz -L $URL
        gzip -d ${FILE}.gz
    fi
}

function download() {
    URL=$2
    FILE=$1
    if [ ! -e $FILE ]; then
        curl -o ${FILE} -L $URL
    fi
}

function create_faidx() {
    if [ ! -e ${1}.fai ];then
        samtools faidx $1
    fi
}

download_gzip inputs/HX1.1.fa ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/708/065/GCA_001708065.2_HX1_1.1/GCA_001708065.2_HX1_1.1_genomic.fna.gz
download_gzip inputs/AK1v2.fa ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/750/385/GCA_001750385.2_AK1_v2/GCA_001750385.2_AK1_v2_genomic.fna.gz
download_gzip inputs/JG1.fa https://jmorp.megabank.tohoku.ac.jp/dj1/datasets/tommo-jg1.0.0.beta-20190424/files/JG1.0.0beta.fa.gz
download_gzip inputs/CHM1_1.1.fa ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/306/695/GCA_000306695.2_CHM1_1.1/GCA_000306695.2_CHM1_1.1_genomic.fna.gz
download_gzip inputs/hg19.fa http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
download_gzip inputs/hg38.fa http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

create_faidx inputs/HX1.1.fa
create_faidx inputs/AK1v2.fa
create_faidx inputs/JG1.fa
create_faidx inputs/CHM1_1.1.fa
create_faidx inputs/hg19.fa
create_faidx inputs/hg38.fa

download inputs/common_all_20180423.vcf.gz ftp://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
download inputs/common_all_20180423.vcf.gz.tbi ftp://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi

download inputs/common_all_20180418_GRCh38.vcf.gz ftp://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
download inputs/common_all_20180418_GRCh38.vcf.gz.tbi ftp://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi


download inputs/clinvar_20190923.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190923.vcf.gz
download inputs/clinvar_20190923.vcf.gz.tbi ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190923.vcf.gz.tbi

download inputs/clinvar_20190923_GRCh38.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20190923.vcf.gz
download inputs/clinvar_20190923_GRCh38.vcf.gz.tbi ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20190923.vcf.gz.tbi

for one in common_all_20180423 common_all_20180418_GRCh38 clinvar_20190923 clinvar_20190923_GRCh38; do
    if [ ! -e inputs/${one}.snps.vcf.gz ]; then
        bcftools view -o inputs/${one}.snps.vcf.gz -O z inputs/${one}.vcf.gz -v snps
        bcftools index -t inputs/${one}.snps.vcf.gz
    fi

    if [ ! -e inputs/${one}.indels.vcf.gz ]; then
        bcftools view -o inputs/${one}.indels.vcf.gz -O z inputs/${one}.vcf.gz -v indels
        bcftools index -t inputs/${one}.indels.vcf.gz
    fi
done

for one in HX1.1 AK1v2 JG1 hg38 CHM1_1.1; do
    if [ ! -e inputs/hg19-to-${one}.paf ];then
        minimap2 -cx asm5 --cs inputs/${one}.fa inputs/hg19.fa  > inputs/hg19-to-${one}.paf
    fi
done
