ref=/sc/orga/projects/ad-omics/data/references/
zless $ref/hg38_reference/GENCODE/gencode.v30.annotation.gtf.gz | awk ' $1 == "chr1" && $4 >= 108364568 && $5 <= 112364535'  > test_chr1.gtf 
