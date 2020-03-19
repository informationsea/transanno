#!/usr/bin/env bash

set -xe

if [ ! -e testfiles/genome/hg19-chr22.fa ];then
    xz -dk testfiles/genome/hg19-chr22.fa.xz
fi

if [ ! -e testfiles/genome/hg38-chr22.fa ];then
    xz -dk testfiles/genome/hg38-chr22.fa.xz
fi
