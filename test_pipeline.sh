ml python/3.7.3
ml R/3.6.0

cd test/
Rscript create_test_gct_files.R  
sh create_test_gtf.sh  
sh create_test_vcf.sh
cd ..
snakemake
