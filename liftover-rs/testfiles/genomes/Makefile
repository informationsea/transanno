.PHONY: all clean

all: GRCh38/GRCh38.chr22.genome.fa.xz \
     GRCh37/GRCh37.chr22.genome.fa.xz \
     GRCh37/GRCh37.chr22.revcomp.genome.fa.xz \
     chain/GRCh38-to-GRCh37.chr22.chain chain/GRCh38-to-GRCh37.chr22.revcomp.chain chain/GRCh37-to-GRCh38.chr22.chain \
     chain/GRCh38-to-GRCh37.over.bed chain/GRCh38-to-GRCh37.over.vcf chain/GRCh38-to-GRCh37.over.sorted.vcf.gz \
     chain/GRCh37-to-GRCh38.over.bed chain/GRCh37-to-GRCh38.over.vcf chain/GRCh37-to-GRCh38.over.sorted.vcf.gz \
     chain/GRCh38-to-GRCh37.revcomp.over.bed chain/GRCh38-to-GRCh37.revcomp.over.vcf chain/GRCh38-to-GRCh37.revcomp.over.sorted.vcf.gz \
     chain/GRCh38-to-GRCh37.test-position.json \
     chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.bed \
     chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.bed \
     chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.json \
     chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.json

clean:
	-rm */*.fa
	-rm */*.vcf.gz
	-rm */*.chain
	-rm */*.paf
	-rm */*.vcf
	-rm */*.bed
	-rm */*.tbi
	-rm */*.fai
	-rm */*.xz
	-rm */*.json
	-rm chain/mapped-region/*.bed

%.sorted.vcf.gz: %.vcf
	bcftools sort --output-file $@ --output-type z $<
	bcftools index -t $@

chain/GRCh38-to-GRCh37.over.bed chain/GRCh38-to-GRCh37.over.vcf: scripts/chain-to-bed-vcf.py GRCh37/GRCh37.chr22.genome.fa GRCh38/GRCh38.chr22.genome.fa chain/GRCh38-to-GRCh37.chr22.chain
	python3 scripts/chain-to-bed-vcf.py --track-name 'GRCh38 to GRCh37' --bed-output chain/GRCh38-to-GRCh37.over.bed --vcf-output chain/GRCh38-to-GRCh37.over.vcf --reference-fasta GRCh38/GRCh38.chr22.genome.fa --target-fasta GRCh37/GRCh37.chr22.genome.fa chain/GRCh38-to-GRCh37.chr22.chain

chain/GRCh38-to-GRCh37.revcomp.over.bed chain/GRCh38-to-GRCh37.revcomp.over.vcf: scripts/chain-to-bed-vcf.py GRCh37/GRCh37.chr22.revcomp.genome.fa GRCh38/GRCh38.chr22.genome.fa chain/GRCh38-to-GRCh37.chr22.revcomp.chain
	python3 scripts/chain-to-bed-vcf.py --track-name 'GRCh38 to GRCh37 revcomp' --bed-output chain/GRCh38-to-GRCh37.revcomp.over.bed --vcf-output chain/GRCh38-to-GRCh37.revcomp.over.vcf --reference-fasta GRCh38/GRCh38.chr22.genome.fa --target-fasta GRCh37/GRCh37.chr22.revcomp.genome.fa chain/GRCh38-to-GRCh37.chr22.revcomp.chain


chain/GRCh37-to-GRCh38.over.bed chain/GRCh37-to-GRCh38.over.vcf: scripts/chain-to-bed-vcf.py GRCh37/GRCh37.chr22.genome.fa GRCh38/GRCh38.chr22.genome.fa chain/GRCh37-to-GRCh38.chr22.chain
	python3 scripts/chain-to-bed-vcf.py --track-name 'GRCh37 to GRCh38' --bed-output chain/GRCh37-to-GRCh38.over.bed --vcf-output chain/GRCh37-to-GRCh38.over.vcf --reference-fasta GRCh37/GRCh37.chr22.genome.fa --target-fasta GRCh38/GRCh38.chr22.genome.fa chain/GRCh37-to-GRCh38.chr22.chain

chain/GRCh38-to-GRCh37.test-position.json: chain/GRCh38-to-GRCh37.chr22.chain scripts/create-currect-lift-position.py
	python3 scripts/create-currect-lift-position.py -o $@ $< chr22 50818468

chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.json: chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.bed chain/GRCh38-to-GRCh37.test-region.GRCh38.bed scripts/create-test-regions-json.py
	python3 scripts/create-test-regions-json.py -o $@ chain/GRCh38-to-GRCh37.test-region.GRCh38.bed chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.bed

chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.json: chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.bed chain/GRCh38-to-GRCh37.test-region.GRCh38.bed scripts/create-test-regions-json.py
	python3 scripts/create-test-regions-json.py -o $@ chain/GRCh38-to-GRCh37.test-region.GRCh38.bed chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.bed


chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.bed: chain/GRCh38-to-GRCh37.test-region.GRCh38.bed chain/GRCh38-to-GRCh37.chr22.chain liftOver
	./liftOver -minMatch=0.9 -multiple -noSerial $< chain/GRCh38-to-GRCh37.chr22.chain chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.bed chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.unmapped.bed

chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.bed: chain/GRCh38-to-GRCh37.test-region.GRCh38.bed chain/GRCh38-to-GRCh37.chr22.revcomp.chain liftOver
	./liftOver -minMatch=0.9 -multiple -noSerial $< chain/GRCh38-to-GRCh37.chr22.revcomp.chain chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37-revcomp.bed chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.GRCh37-revcomp.unmapped.bed

chain/GRCh38-to-GRCh37.test-region.GRCh38.bed: scripts/create-test-regions.py
	python3 scripts/create-test-regions.py -o $@ chr22 50818468

%.chain: %.paf scripts/minimap2-to-liftover-chain.py
	python3 scripts/minimap2-to-liftover-chain.py --output $@ $<

chain/GRCh38-to-GRCh37.chr22.paf: GRCh37/GRCh37.chr22.genome.fa GRCh38/GRCh38.chr22.genome.fa minimap2
	./minimap2 -cx asm5 --cs GRCh37/GRCh37.chr22.genome.fa GRCh38/GRCh38.chr22.genome.fa > $@

chain/GRCh38-to-GRCh37.chr22.revcomp.paf: GRCh37/GRCh37.chr22.revcomp.genome.fa GRCh38/GRCh38.chr22.genome.fa minimap2
	./minimap2 -cx asm5 --cs GRCh37/GRCh37.chr22.revcomp.genome.fa GRCh38/GRCh38.chr22.genome.fa > $@

chain/GRCh37-to-GRCh38.chr22.paf: GRCh37/GRCh37.chr22.genome.fa GRCh38/GRCh38.chr22.genome.fa minimap2
	./minimap2 -cx asm5 --cs GRCh38/GRCh38.chr22.genome.fa GRCh37/GRCh37.chr22.genome.fa > $@

%.xz: %
	xz -k $<

%.chr22.genome.fa: %.primary_assembly.genome.fa.gz scripts/extract-chromosome.py
	python3 scripts/extract-chromosome.py $< $@ chr22 GL000217.1 KI270879.1

GRCh37/GRCh37.chr22.revcomp.genome.fa: GRCh37/GRCh37.chr22.genome.fa scripts/reverse-complement.py
	python3 scripts/reverse-complement.py $< $@

GRCh38/GRCh38.primary_assembly.genome.fa.gz:
	curl -o $@ -L "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz"

GRCh37/GRCh37.primary_assembly.genome.fa.gz:
	curl -o $@  -L "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"

minimap2:
	curl -OL https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
	tar xjf minimap2-2.17_x64-linux.tar.bz2
	cp minimap2-2.17_x64-linux/minimap2 .
	rm -rf minimap2-2.17_x64-linux

liftOver:
	curl -OL "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"
	chmod +x liftOver
