ml samtools/1.9
wget http://jungle.unige.ch/QTLtools_examples/HG00381.chr22.bam
mv HG00381.chr22.bam HG00381.bam
samtools index HG00381.bam
rm HG00381.chr22.bam

ml tabix
wget http://jungle.unige.ch/QTLtools_examples/genotypes.chr22.vcf.gz
tabix genotypes.chr22.vcf.gz
