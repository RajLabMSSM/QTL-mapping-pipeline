ml bedtools

ref=/sc/arion/projects/ad-omics/data/references/
bedtools intersect -b regions.bed -a $ref/hg38_reference/GENCODE/gencode.v30.annotation.gtf.gz > test.gtf

