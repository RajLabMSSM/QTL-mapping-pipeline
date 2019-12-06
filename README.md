# QTL mapping pipeline

authors: Sanan Venkatesh and Jack Humphrey
Towfique Raj lab, Mount Sinai


## dependencies:

Python:

`ml python/3.7.3`

R

`ml R/3.6.0`


## running on test data

```
# just run sh test_pipeline.sh
cd test

sh create_test_gtf.sh  
Rscript create_test_gct_files.R  
sh create_test_vcf.sh

cd ..
snakemake
```
