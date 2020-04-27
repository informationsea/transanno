all: dbsnp-151-common-GRCh38.vcf.gz dbsnp-151-common-GRCh37.vcf.gz \
     dbsnp-151-common-GRCh38.chr22.vcf.gz dbsnp-151-common-GRCh37.chr22.vcf.gz


%.chr22.vcf.gz: %.vcf.gz
	zegrep '^(#|22)' $< | bgzip -c > $@
	bcftools index -t $@

dbsnp-151-common-GRCh38.vcf.gz:
	curl -o dbsnp-151-common-GRCh38.vcf.gz -L https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz

dbsnp-151-common-GRCh37.vcf.gz:
	curl -o dbsnp-151-common-GRCh37.vcf.gz -L https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz

.PHONY: all