#!/bin/bash

set -xeu -o pipefail

script_dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
cd "${script_dir}"

DOWNLOAD_DIR=testfiles/ucsc

function download() {
    FILENAME="$(basename $1)"
    if [ ! -f "${DOWNLOAD_DIR}/${FILENAME}" ]; then
        curl -o "${DOWNLOAD_DIR}/${FILENAME}" "$1"
    fi
}

mkdir -p "${DOWNLOAD_DIR}"

download https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
download https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

if [ "$(uname)" == 'Darwin' ]; then
    download http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver
else
    download http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
fi
chmod +x "${DOWNLOAD_DIR}/liftOver"

download https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
download https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

if [ ! -f "${DOWNLOAD_DIR}/hg19.fa" ]; then
    gzip -dc "${DOWNLOAD_DIR}/hg19.fa.gz" > "${DOWNLOAD_DIR}/hg19.fa"
fi

if [ ! -f "${DOWNLOAD_DIR}/hg38.fa" ]; then
    gzip -dc "${DOWNLOAD_DIR}/hg38.fa.gz" > "${DOWNLOAD_DIR}/hg38.fa"
fi

test ! -f "${DOWNLOAD_DIR}/hg19.fa.fai" && samtools faidx "${DOWNLOAD_DIR}/hg19.fa"
test ! -f "${DOWNLOAD_DIR}/hg38.fa.fai" && samtools faidx "${DOWNLOAD_DIR}/hg38.fa"
test ! -f "${DOWNLOAD_DIR}/hg19-regions.bed" && python3 testfiles/ucsc-scripts/create-test-regions.py testfiles/ucsc/hg19.fa.fai --output testfiles/ucsc/hg19-regions.bed
test ! -f "${DOWNLOAD_DIR}/hg38-regions.bed" && python3 testfiles/ucsc-scripts/create-test-regions.py testfiles/ucsc/hg38.fa.fai --output testfiles/ucsc/hg38-regions.bed

./testfiles/ucsc/liftOver -multiple -minMatch=0.1 -bedPlus=6 -tab "${DOWNLOAD_DIR}/hg19-regions.bed" "${DOWNLOAD_DIR}/hg19ToHg38.over.chain.gz" "${DOWNLOAD_DIR}/hg19-regions-mapped-to-hg38.bed" "${DOWNLOAD_DIR}/hg19-regions-mapped-to-hg38-unmapped.bed"
./testfiles/ucsc/liftOver -multiple -minMatch=0.1 -bedPlus=6 -tab "${DOWNLOAD_DIR}/hg38-regions.bed" "${DOWNLOAD_DIR}/hg38ToHg19.over.chain.gz" "${DOWNLOAD_DIR}/hg38-regions-mapped-to-hg19.bed" "${DOWNLOAD_DIR}/hg38-regions-mapped-to-hg19-unmapped.bed"

python3 testfiles/ucsc-scripts/create-result-json.py --output >(gzip -c > testfiles/ucsc/hg38-regions-mapped-to-hg19.jsonl.gz) testfiles/ucsc/hg38-regions-mapped-to-hg19.bed testfiles/ucsc/hg38-regions-mapped-to-hg19-unmapped.bed
python3 testfiles/ucsc-scripts/create-result-json.py --output >(gzip -c > testfiles/ucsc/hg19-regions-mapped-to-hg38.jsonl.gz) testfiles/ucsc/hg19-regions-mapped-to-hg38.bed testfiles/ucsc/hg19-regions-mapped-to-hg38-unmapped.bed
