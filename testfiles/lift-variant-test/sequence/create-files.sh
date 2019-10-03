python3 concat-fasta.py -o seq-a.fa seq-a/gencode.v30.*.fa -n seq-a
python3 concat-fasta.py -o seq-b.fa seq-b/gencode.v30.*.fa -n seq-b
./minimap2 -cx asm5 --cs seq-b.fa seq-a.fa > seq-a--to--seq-b.paf
rm seq-a--to--seq-b.chain.gz
python3 minimap2-to-liftover-chain.py --output-chain-to-target seq-a--to--seq-b.chain seq-a--to--seq-b.paf
python3 ../../chain-to-bed-vcf.py --track-name 'test chain' --bed-output seq-a--to--seq-b.bed --vcf-output seq-a--to--seq-b.vcf -r seq-a.fa -t seq-b.fa seq-a--to--seq-b.chain --diff-bed-output seq-a--to--seq-b.diff.bed
bcftools sort -O z -o seq-a--to--seq-b.sorted.vcf.gz seq-a--to--seq-b.vcf
bcftools index -t seq-a--to--seq-b.sorted.vcf.gz
bcftools norm -O z -o seq-a--to--seq-b.sorted.norm.vcf.gz -f ./seq-a.fa seq-a--to--seq-b.sorted.vcf.gz
bcftools index -t seq-a--to--seq-b.sorted.norm.vcf.gz
rm seq-a--to--seq-b.vcf
gzip seq-a--to--seq-b.chain
