QTL mapping pipeline
=====================

authors: Sanan Venkatesh and Jack Humphrey
Towfique Raj lab, Mount Sinai

This pipeline is based on the [GTEx QTL mapping pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) largely authored by Francois Auguet. A number of scripts from that pipeline have been copied over here. All python scripts have been ported to python3 to work with Snakemake.

I cannot guarantee that this pipeline will work with your data. Please raise any bugs or errors as issues on this repo.

Jack


## conda recipe

YAML files are now procided in the repo.

tensorQTL is now in it's own environment 

```
conda env create -f QTL-pipeline.yml
conda env create -f tensorqtl.yml
```

## TensorQTL dev version

You must install the development version of TensorQTL as the pip release has several bugs that break nominal QTL mapping.
Therefore you must manually install this from where we keep it in the communal software folder once you've got your conda environment:

(This assumes you have access to ad-omics), otherwise clone the repo from https://github.com/broadinstitute/tensorqtl

```
conda create -n tensorqtl python=3.7 pip snakemake
conda activate tensorqtl
ml R/3.6.2 # or the version you are using
pip3 install --upgrade tensorqtl

#cd /sc/arion/projects/H_ad-omics/data/software/tensorqtl
#pip install -r install/requirements.txt .
#pip install --upgrade pandas
#pip install --upgrade cython
#pip install --upgrade rpy2
#pip install --upgrade snakemake
```

## running on test data

```
# just run sh test_pipeline.sh
conda activate QTL-pipeline
cd test

sh create_test_gtf.sh  
sh create_test_vcf.sh

cd ..
snakemake --configfile test/config.yaml -pr --config mode=eQTL --cores 4
snakemake --configfile test/config.yaml -pr --config mode=sQTL --cores 4
# Runs with all modes (interaction + trans + conditional)
snakemake --configfile test/config.yaml -pr --config mode=eQTL interaction=True trans=True conditional_qtls=True --cores 4
```

# Input Files

Examples of each input file are in `example/`. The master file is the `config.yaml`, which sets the paths to all other input files. 

#### `dataCode:`  
  The name of the experiment or dataset. This will be used as prefix name for the results folder created inside `results/`. 
  
#### `sampleKey:` 
  The sample key file contains two columns, sample_id and participant_id. sample_id refers to the RNA-seq samples and must match the column names in your count matrix and the names of your junction files. participant_id refers to the donor ID in your VCF and must match. 

  See: example/Cerebellum_sample_key.txt
  
#### `genotypePCs:` 
  A matrix of genotype PCs where the columns are first `ID` and then the sample IDs (matching the `sample_id` column of the sample key), and the rows are principal components of the genotypes. The first column gives the name of the PCs. 
 
 See example/Cerebellum_genotype_PCs.txt 
 
#### `covariateFile:`
  A matrix of other covariates you want to include such as sex and sequencing platform. Again the columns should be ID and the sample_id values and the rows should be the covariates.
  
  See example/Cerebellum_wgs_covariates.txt 

#### `junctionFileList:`
   This file is a list of paths to junction files for splicing QTLs. These junctions are created using regtools, which comes as part of the [leafcutter-pipeline](https://github.com/RajLabMSSM/leafcutter-pipeline).
  
   See: example/Cerebellum_junction_paths.txt 
  
#### `PEER_values: `
  Number of PEER values to use. Can be a single number of a list of numbers. If a list of numbers provided then QTL mapping will be performed with each value of PEER factors. Lists are encoded in YAML like so:

  ```
  PEER_values:
    - 15
    - 20
    - 25
    - 30
  ```
  
#### `VCF:` 
  Path to filtered joint VCF.

#### `countMatrixRData:` 
  Path to count matrix, an RData file created by collating gene expression values from RSEM. I use the output of [collate_samples](https://github.com/RajLabMSSM/RNA-pipelines/tree/master/collate_samples).
  This matrix can contain more samples than you plan to analyse as the matrix will be filtered.

#### `interaction:` 
  (optional) Flag for interaction-eQTL analysis. If set `True` (default is `False`), other the parameters `interaction_name` and `interaction_file` must be specified. See more below.

#### `interaction_name:`
  A name identifing the interaction. This string will be used to name the results files, allowing multiple interaction runs using the same data (just need to update this parameter). You can also, specify a list with multiple interaction files to run in parallel, just make sure to match the `interaction_name` with `interaction_file`. See `test/config.yaml` for an example.

#### `interaction_file:`
  A text file with 2 columns. The first column has the participant_id (matching sampleKey file) and the second has the interaction values (must be numeric!).

#### `trans:`
  (optional) Flag for trans-eQTL analysis. If set `True`, in addition to cis results, a file with nominal associations between all phenotypes and genotypes will be generated (default is `False`). The output is in txt.gz format, with four columns: phenotype_id, variant_id, pval, maf. As default only associations with p-value < 1e-4 will be reported (change in the Snakefile if necessary).

#### `conditional_qtls:`
  (optional) Flag for conditional-eQTL analysis. This mode maps conditionally independent cis-QTLs using the stepwise regression procedure. It will run on top of the cis-eQTL mapping results.

# Selecting modes

The QTL-mapping-pipeline currently supports 3 modes of execution:
 * eQTL - map expression quantitative trait loci
 * sQTL - map splicing quantitative trait loci
 * mbv - genotype your RNA samples and match to your VCF using `QTLtools mbv`.

The mode of execution is set by the `mode` variable. This can either be hardcoded in your config.yaml or can be set when you run snakemake using:

```
snakemake -s Snakefile --configfile config.yaml --config mode=<mode>
```

The _mbv_ mode requires a different config.yaml, an example of which is in `example/`. It requires the following:

* bamFolder:
  A folder containing all BAM files you want to match
* bamSuffix: ".bam"
  The suffix added to sample_id to make the name of the file. Assumed to be '.bam'.
* dataCode: 
  The name of the dataset or experiment
* VCF: 
  The path to the VCF file.
  
 # Parallel execution and chunking
 
 The pipeline has been designed to be run in parallel on the MSSM HPC cluster. The script that connects the snakemake code to the LSF job scheduler is found in `snakejob`.  If you do not have access to the `als-omics` account then you will have to modify this. 
 
 `snakejob` has the following options:
  * -s the relative path to the Snakefile.
  * -c the relative path to the config.yaml.
  * -n an optional flag to set if you want to do a [dry run](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#useful-command-line-arguments) - advisable before executing the full pipeline.
  * -m which mode you want to use.

Example:

```
./snakejob -s Snakefile -c config.yaml -m sQTL -n
```

The resources each step is allocated is set in `cluster.yaml`. 

# UNLIMITED POWER - run multiple QTL analyses simultaneously

Building upon `snakejob`, the script `run_all_QTLs_parallel.sh` takes as input list of paths to config files. For each file it then submits a separate snakejob job to the cluster for both eQTLs and sQTLs. So all your jobs can be run at the same time!

usage:

```
sh run_all_QTLs_parallel.sh -c config_file_list.txt -s Snakefile
``` 
 
# To be documented:

 # Features coming soon
  
