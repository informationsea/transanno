.PHONY: all clean clean-test

all: gencode.v33lift37.basic.annotation.chr22.test.regions.gff3 gencode.v33.basic.annotation.chr22.test.regions.gff3 \
     gencode.v33.basic.annotation.chr22.gff3.zst gencode.v33lift37.basic.annotation.chr22.gff3.zst \
	 gencode.v33.basic.annotation.chr22.gtf.zst gencode.v33lift37.basic.annotation.chr22.gtf.zst \
	 Homo_sapiens.GRCh38.99.ensembl.chr22.gtf.zst Homo_sapiens.GRCh38.99.ensembl.chr22.gff3.zst
clean:
	rm *.gff3

clean-test:
	rm *.test.regions.gff3

%.test.regions.gff3: %.gff3.zst
	zstd -dc $< | egrep '(ENST00000343518|ENST00000610390|ENST00000400585)'  > $@

%.chr22.gff3.zst: %.gff3.gz
	zegrep '^chr22' $< | zstd -c > $@

%.chr22.gtf.zst: %.gtf.gz
	zegrep '^chr22' $< | zstd -c > $@

%.ensembl.chr22.gff3.zst %.chr.gff3.gz
	zegrep '^22' $< | zstd -c > $@

%.ensembl.chr22.gtf.zst: %.chr.gtf.gz
	zegrep '^22' $< | zstd -c > $@


%: %.zst
	zstd -dk $<	

gencode.v33lift37.basic.annotation.gff3.gz:
	curl -OL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.basic.annotation.gff3.gz

gencode.v33lift37.basic.annotation.gtf.gz:
	curl -OL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.basic.annotation.gtf.gz

gencode.v33.basic.annotation.gff3.gz:
	curl -OL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.basic.annotation.gff3.gz

gencode.v33.basic.annotation.gtf.gz:
	curl -OL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.basic.annotation.gtf.gz
	
Homo_sapiens.GRCh38.99.chr.gtf.gz:
	curl -OL ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.chr.gtf.gz

Homo_sapiens.GRCh38.99.chr.gff3.gz:
	curl -OL ftp://ftp.ensembl.org/pub/release-99/gff3/homo_sapiens/Homo_sapiens.GRCh38.99.chr.gff3.gz