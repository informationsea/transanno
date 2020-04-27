all: gnomad.genomes.r3.0.sites.chr22.vcf.gz gnomad.genomes.r3.0.sites.chr22.vcf.gz.tbi \
     gnomad.genomes.r3.0.sites.chr22.aconly.vcf.gz gnomad.genomes.r3.0.sites.chr22.aconly.vcf.gz.tbi \
	 gnomad.genomes.r3.0.sites.chr22.aconly.subset.vcf.gz gnomad.genomes.r3.0.sites.chr22.aconly.subset.vcf.gz.tbi \
	 gnomad.genomes.r3.0.sites.chr22.aconly.subset.merged.vcf.gz gnomad.genomes.r3.0.sites.chr22.aconly.subset.merged.vcf.gz.tbi \
	 gnomad.genomes.r3.0.sites.header.vcf

clean:
	-rm	*.merged.vcf.gz
	-rm	*.subset.vcf.gz
	-rm	*.aconly.vcf.gz
	-rm *.tbi

%.merged.vcf.gz: %.vcf.gz %.vcf.gz.tbi
	bcftools norm -m +any -O z -o $@ $<
	bcftools index -t $@

%.subset.vcf.gz: %.vcf.gz %.vcf.gz.tbi
	bcftools view -f PASS -O z -o $@ -r chr22:15687860-15725139,chr22:24063387-24088324 $<
	bcftools index -t $@

%.aconly.vcf.gz: %.vcf.gz
	./vcfutils rewrite-info -i AC AN AF -- $< | bgzip -c > $@
	bcftools index -t $@

gnomad.genomes.r3.0.sites.header.vcf: gnomad.genomes.r3.0.sites.chr22.vcf.gz
	zegrep '^#' $< > $@

gnomad.genomes.r3.0.sites.chr22.vcf.gz:
	curl -o gnomad.genomes.r3.0.sites.chr22.vcf.gz -L https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.chr22.vcf.bgz

gnomad.genomes.r3.0.sites.chr22.vcf.gz.tbi:
	curl -o gnomad.genomes.r3.0.sites.chr22.vcf.gz.tbi -L https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.chr22.vcf.bgz.tbi

.PHONY: all clean