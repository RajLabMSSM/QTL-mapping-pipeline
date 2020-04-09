QTL mapping pipeline
=====================

authors: Sanan Venkatesh and Jack Humphrey
Towfique Raj lab, Mount Sinai

This pipeline is based on the [GTEx QTL mapping pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) largely authored by Francois Auguet. A number of scripts from that pipeline have been copied over here. All python scripts have been ported to python3 to work with Snakemake.

I cannot guarantee that this pipeline will work with your data. Please raise any bugs or errors as issues on this repo.

Jack


## conda recipe

YAML file is now in the repo 

```
conda env create -f QTL-pipeline.yml
```

## TensorQTL dev version

You must install the development version of TensorQTL as the pip release has several bugs that break nominal QTL mapping.
Therefore you must manually install this from where we keep it in the communal software folder once you've got your conda environment:

(This assumes you have access to ad-omics), otherwise clone the repo from https://github.com/broadinstitute/tensorqtl

```
conda activate QTL-pipeline
cd /sc/hydra/projects/ad-omics/data/software/tensorqtl
pip install -r install/requirements.txt .

```

## running on test data

```
# just run sh test_pipeline.sh
conda activate QTL-pipeline
cd test

sh create_test_gtf.sh  
sh create_test_vcf.sh

cd ..
snakemake --configfile test/config.yaml -pr --config mode=eQTL
```


# Input Files

Examples of each input file are in `example/`. The master file is the `config.yaml`, which sets the paths to all other input files. 


#### `dataCode:`  
  The name of the experiment or dataset
  
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
 

# To be documented:

* Interaction QTLs

 # Features coming soon
 
 * trans-QTLs using tensorQTL
  
