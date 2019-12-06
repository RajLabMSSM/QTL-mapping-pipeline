# QTL mapping pipeline

authors: Sanan Venkatesh and Jack Humphrey
Towfique Raj lab, Mount Sinai


## conda recipe

```
conda create -c bioconda -c conda-forge -n QTL-pipeline python=3.7 snakemake=5.8 tabix=0.2.5 bcftools=1.9
pip install feather-format
pip install bx-python
pip install scipy
```

## running on test data

```
# just run sh test_pipeline.sh
conda activate QTL-pipeline
cd test

sh create_test_gtf.sh  
sh create_test_vcf.sh

cd ..
snakemake --configfile test/config.yaml -pr
```


# Match BAMS to VCF

## running test data

```
cd test/match_bam_to_genotypes/
sh download_test_files.sh
cd ../../

conda activate QTL-pipeline
snakemake -s match_bams_to_genotypes.smk --configfile test/match_bam_to_genotypes/config.yaml -pr
# results files will be in test/match_bam_to_genotypes/results
```
