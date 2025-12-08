# transanno

[![Build](https://github.com/informationsea/transanno/actions/workflows/rust.yml/badge.svg)](https://github.com/informationsea/transanno/actions/workflows/rust.yml)
[![GitHub](https://img.shields.io/github/license/informationsea/transanno)](https://www.gnu.org/licenses/gpl-3.0.html)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/informationsea/transanno)](https://github.com/informationsea/transanno)

accurate VCF/GFF3/GTF LiftOver tool for new genome assemblies

## Download

* [GitHub Release](https://github.com/informationsea/transanno/releases)
* [Docker Hub](https://hub.docker.com/r/informationsea/transanno)

## Author

* OKAMURA, Yasunobu

## License

* GPL 3 or later
* Parts of dependent libraries are licensed under Apache 2.0/MIT

## Features

* Create chain file from [minimap2](https://github.com/lh3/minimap2) result
* LiftOver
    * VCF file (reorder INFO and FORMAT tags, and rewrite GT data)
    * GENCODE/Ensembl GFF3/GTF
* Create VCF and BED from chain file to compare two genome assemblies

## Results

### case 1: map [dbSNP 151 common](ftp://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz) from GRCh37 to GRCh38

| Tool    | Mapped variants| Correctly mapped variants| Incorrectly mapped variants| Unmapped variants| Incorrectly unmapped variants|
|---------|---------------:|-------------------------:|---------------------------:|-----------------:|-----------------------------:|
|transanno|      37,892,118|                37,256,006|                      29,730|            14,713|                         8,669|
|Picard   |      37,807,777|                37,182,636|                      28,765|            99,054|                        83,004|

* A chain file was created with minimap2 and transanno.
* A definition of correctly mapped variant is [dbSNP 151 common GRCh38 version](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20190923.vcf.gz).
* Since not all variants are included in GRCh38 version of dbSNP, a summation of a number of correctly mapped variant and a number of incorrectly mapped variant is not equal to a number of mapped variant.

### case 2: map dbSNP and ClinVar to AK1v2 from GRCh37

| Tool    | Dataset        |Mapped Variants|Unmapped Variants|
|---------|----------------|--------------:|----------------:|
|transanno|dbSNP 151 common|     37,518,510|          388,321|
|Picard   |dbSNP 151 common|     36,478,197|        1,428,634|
|CrossMap |dbSNP 151 common|     34,726,126|        3,180,705|

## How to use

### Create chain file from minimap2 result

1. Prepare two genome assemblies.
   In this document, we call reference as a source assembly that original VCF/GFF3/GTF coordinates, and query as a target assembly that new VCF/GFF3/GTF coordinates.
2. Run minimap2
    * `minimap2 -cx asm5 --cs QUERY_FASTA.fa REFERENCE_FASTA.fa > PAF_FILE.paf`
3. Run transanno to create chain file
    * `transanno minimap2chain PAF_FILE.paf --output CHAINFILE.chain`

Thanks to minimap2, you can create new chain file in 30 minute.

### Convert VCF File

1. Prepare a VCF file, a query FASTA, a reference FASTA, a chain file.
    * You do not need to add `chr` prefix to a contig name
    * Index files for FASTA are required (create it with `samtools faidx`)
2. Run transanno to convert coordinates
    * `liftvcf -m --chain CHAINFILE.chain  -o SUCCEEDED.vcf.gz --query QUERY_FASTA.fa --reference REFERENCE_FASTA.fa --vcf INPUT_VCF.vcf.gz --fail FAILED.vcf.gz`
    * Input files can be compressed with gzip if a file name is ends with `.gz` or `.bgz`
    * transanno will compress output files if a file name is ends with `.gz`

#### VCF Notes

* transanno can handle correctly multiallelic sites, asterisk, dot in an ALT column.
* Output files are not compressed with bgzip
    * Please recompress with `bgzip` to add index to `FAILED.vcf.gz`
* Output files are not sorted.
    * If you want to load converted VCF files into genome browser, please sort and create index first.
       * e.g. `bcftools sort -O z -o SUCCEEDED.sorted.vcf.gz SUCCEEDED.vcf.gz && bcftools index -t SUCCEEDED.sorted.vcf.gz`
* transanno swaps REF and ALT if reference allele was changed.
    * Use `--noswap` to disable swapping REF and ALT.
    * If you want to convert ClinVar or COSMIC, `--noswap` option is recommended.
* transanno rewrite allele frequency, INFO tags, FORMAT tags and a GT tag if REF and ALT columns are swapped.
    * If you want to disable these rewriting these fields, please add options.
* transanno prefer same contig name if a variant was multi-mapped.
    * When a variant in hg19 `chr11` can be converted to hg38 `chr11` and `chr11_KI270721v1_random`, transanno will convert to hg38 `chr11`.
    * If you want to disable this behavior, add `--do-not-prefer-same-contig-when-multimap` options.
* Multi-mapped variants will be failed to convert.
    * If you want to allow multi-map, add `--allow-multi-map`
* transanno perform left align to chain before converting for accurate conversion
    * If you want to disable this behavior, add `--no-left-align-chain`
    * This step takes few seconds if you have a fast storage.
    * You can create lift aligned chain file with `transanno chain-left-align`, but output of this command can be invalid chain.

### Convert GENCODE/Ensembl GFF3/GTF

1. Prepare GENCODE or Ensembl GFF3/GTF file, a query FASTA, a reference FASTA, a chain file.
2. Run transanno
   `transanno liftgene GENCODE.gtf.gz --chain CHAINFILE.chain --failed FAILED.gtf.gz --output SUCCEEDED.gtf.gz`

#### GFF3/GTF Notes

* transanno is not general converter for GFF3/GTF file.
    * GFF3/GTF file must be fetched from GENCODE or Ensembl.
    * If you want to convert GFF3/GFT entries independently, please use official UCSC's liftOver tool.
* transanno do not modify frame column when CDS contains INDEL
    * This is one of known issue, and will be fixed

### Create VCF and BED from chain file

1. Prepare a query FASTA, a reference FASTA, a chain file.
    * FASTA files must be indexed.
2. Run transanno
   `transanno chain-to-bed-vcf CHAINFILE.chain.gz --output-query-bed QUERY_REGION.bed --query QUERY_FASTA.fa --output-query-vcf QUERY_DIFF.vcf.gz --output-reference-bed REFERENCE_REGION.bed --reference REFERENCE_FASTA.fa --output-reference-vcf REFERENCE_DIFF.vcf.gz`

### Notes

* VCF and BED outputs are not sorted.
* VCF and BED files are not compressed with bgzip.

## How to build

1. Install [Rust Toolchain](https://www.rust-lang.org/)
2. Run `./prepare-test-files.sh` to create test files
3. Run `cargo test && cargo build --release`

## Algorithms

### Variant LiftOver

TODO

### Gene LiftOver

1. Lift each exon, CDS, UTR features
    * Lift a region with chain file
    * If length of lifted feature was 50% smaller than original length, mark as failed.
    * If length of lifted feature was 150% larger than original length, mark as failed.
    * If lift over was failed, mark as failed.
    * Multi-map is allowed in this step
2. Lift transcript
    1. Collect lifted exons, CDS, UTR and other features for each transcript.
    2. If one or more feature were marked failed, mark the transcript as failed.
    3. Find a chromosome which is included in all lifted features.
        * If two or more chromosomes were found, mark this transcript as failed.
        * If no chromosome was found, mark this transcript as failed.
        * If a part of features were lifted to multiple chromosome, remove lifted region which is not lifted to selected chromosome from candidates.
    4. Check multi-mapped features
        * If some of features were lifted to multiple region in one chromosome, mark this transcript as failed.
    5. Check the order of exons
        * Do not care order of CDS, UTR, start codon, stop codon in this step.
        * If order of exons was changed from original order, mark this transcript as failed.
    6. Check strand of all features.
        * If any feature have a different strand from other features, mark this transcript as failed.
    7. Make new region from features.
        * Find smallest start position from lifted features and use the position as a transcript start
        * Find largest end position from lifted features and use the position as a transcript end
    8. If length of transcript was 50% smaller than original length, mark as failed.
    9. If length of transcript was 150% larger than original length, mark as failed.
3. Lift gene
    1. Collect lifted transcripts
    2. Collect all chromosome name from lifted transcripts
    3. If two or more chromosomes were found in lifted transcript, pick up most common chromosome.
        * Any transcript lifted into minor chromosome will be mark as failed.
    4. Check strand of all transcripts.
        * If any transcript have a different strand from other transcripts, mark transcripts with minor strand as failed.
    5. Make new region from transcripts.
        * Find smallest start position from lifted transcripts and use the position as a gene start
        * Find largest end position from lifted transcripts and use the position as a gene end
