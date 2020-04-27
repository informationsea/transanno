#!/bin/bash

set -eu

CHROM=$1
START=$2
END=$3

if [ "$3" = "" ]; then
    echo "./lift-helper.sh GRCh38-chromosome GRCh38-start GRCh38-end"
fi

TEMPDIR=$(mktemp -d liftover.XXXXXX)
echo $TEMPDIR
echo "$CHROM $START $END . . +" > $TEMPDIR/original-pos.bed
./liftOver $TEMPDIR/original-pos.bed chain/GRCh38-to-GRCh37.chr22.chain $TEMPDIR/mapped.bed $TEMPDIR/rejected.bed
cat $TEMPDIR/mapped.bed $TEMPDIR/rejected.bed
rm -rf $TEMPDIR