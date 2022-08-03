#!/bin/bash

set -xue

script_dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
cd "${script_dir}"

pushd testfiles/genomes

if [ ! -f GRCh37/GRCh37.chr22.genome.fa ]; then
    xz -dk GRCh37/GRCh37.chr22.genome.fa.xz
fi

if [ ! -f GRCh37/GRCh37.chr22.revcomp.genome.fa ]; then
    xz -dk GRCh37/GRCh37.chr22.revcomp.genome.fa.xz
fi

if [ ! -f GRCh38/GRCh38.chr22.genome.fa ]; then
    xz -dk GRCh38/GRCh38.chr22.genome.fa.xz
fi

popd
