# QTL mapping pipeline
# Jack Humphrey 
import glob
import pandas as pd
import os

leafcutter_dir = "/sc/arion/projects/ad-omics/data/software/leafcutter/"
GTF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.gtf" # cannot be gzipped
GTFexons = GTF + ".exons.txt.gz" 
QTLtools = "/hpc/packages/minerva-centos7/qtltools/1.2/bin/QTLtools"

nPerm = 10000 # number of permutations of the permutation pass

#shell.prefix('export PS1="";source activate QTL-pipeline; ml qtltools/1.2; ml R/3.6.0;')
R_VERSION = "R/3.6.0"
shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge; ml qtltools/1.2; ml {R_VERSION}; conda activate QTL-pipeline;')

# interaction mode - set default to False
if "interaction" not in config.keys():
    config["interaction"] = False
interaction = bool(config["interaction"])

print(interaction)

# common config variables - all modes require these
mode = config["mode"]
dataCode = config["dataCode"]
VCF = config["VCF"]
VCFstem = VCF.split(".vcf.gz")[0]
sample_key = ""
genotypePCs = ""
covariateFile = ""
PEER_values = ""
countMatrixRData = "" 
junctionFileList = ""
phenotype_matrix = ""
bamFolder = ""
bamSuffix = ""
BAM_SAMPLES = []
interaction_string = ""
group_string = ""

if(mode == "mbv"):
    bamFolder = config["bamFolder"]
    bamSuffix = config["bamSuffix"]
    BAM_SAMPLES = [os.path.basename(x).strip(bamSuffix) for x in glob.glob(bamFolder + "/*" + bamSuffix)]
    dataCode = dataCode + "_mbv"
    final_output = prefix + "_mbv_summary.txt"
else:
    sample_key = config["sampleKey"]
    genotypePCs = config["genotypePCs"]
    covariateFile = config["covariateFile"]
    PEER_values = config["PEER_values"]



if( interaction is True ):
    print(" * interaction mode selected")
    if "interaction_name" not in config.keys():
        sys.exit("config.yaml does not contain interaction_name value")

# Interaction QTLs
# requires an interaction_file and an interaction_name in the config.yaml
if(interaction is True):
    interaction_name = config["interaction_name"]
    interaction_file = config["interaction_file"]
    interaction_string = " --interaction {interaction_file}  --maf_threshold_interaction 0.05 "

# expression QTLs
if(mode == "eQTL"):
    #PEER_values = [15,30]
    group_by_values = ["gene"]
    #PEER_values = config["PEER_values"]
    dataCode = dataCode + "_expression" 
    #if(interaction is True):
    #    dataCode = dataCode  + "_interaction_" + interaction_name
    outFolder = "results/" + dataCode + "/"
    prefix = outFolder + dataCode
    phenotype_matrix = prefix + ".expression.bed.gz"
    phenotype_tensorQTL_matrix = prefix + ".phenotype.tensorQTL.bed.gz"
    final_output = [ expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_{group_by}.cis_qtl.txt.gz", PEER_N = PEER_values, group_by = group_by_values), \
                     expand( outFolder + "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}.cis_qtl_nominal_tabixed.tsv.gz", PEER_N = PEER_values, group_by = "gene" ) ]
    countMatrixRData = config["countMatrixRData"]

# splicing QTLs
if(mode == "sQTL"):
    group_by_values = ["cluster"] # should be either 'gene' or 'cluster'
    # PEER values for sQTLs hard-coded at 20
    #PEER_values = [10,20]
    #PEER_values = config["PEER_values"]
    dataCode = dataCode + "_splicing"
    #if(interaction is True):
    #    dataCode = dataCode  + "_interaction_" + interaction_name
    outFolder = "results/" + dataCode + "/"
    prefix = outFolder + dataCode
    junctionFileList = config["junctionFileList"]
    phenotype_matrix = prefix + ".leafcutter.bed.gz"
    phenotype_tensorQTL_matrix = prefix + ".phenotype.tensorQTL.bed.gz"
    final_output = [ expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_{group_by}.cis_qtl.txt.gz", PEER_N = PEER_values, group_by = group_by_values), \
                     expand( outFolder + "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}.cis_qtl_nominal_tabixed.tsv.gz", PEER_N = PEER_values, group_by = "gene" ) ]    
    group_file = prefix + ".group.tensorQTL.{group_by}.bed.gz"
    group_string = " --phenotype_groups " + group_file

print(" * QTL-mapping pipeline *")
print(" Jack Humphrey 2019-2020 ")
print(" * Data code is : %s " % dataCode)
print(" * Mode selected is: %s" % mode)


rule all:
    input:
        final_output

rule collapseGTF:
    input:
        GTF
    output:
        "results/collapsed.gtf"
    params:
        script = "scripts/collapse_annotation.py"
    shell:
        "python {params.script} {input} {output} "

rule getExonsFromGTF:
    input:
        GTF
    params:
        out_folder = os.path.dirname(GTF),
        script = "scripts/get_exon_table.R"
    output:
        GTF + ".exons.txt.gz"
    shell:
        "ml R/3.6.0; "
        "Rscript {params.script} "
        " --gtf {input} "
        " --outFolder {params.out_folder} "

rule createGCTFiles:
    input:
        counts = countMatrixRData,
        key = sample_key,
        gtf = "results/collapsed.gtf"
    output:
        counts_gct_file = prefix + "_counts.gct",
        tpm_gct_file = prefix + "_tpm.gct"
    shell:
        "ml R/3.6.0; "
        "Rscript scripts/create_GCT_files.R "
        " --counts {input.counts} "
        " --key {input.key} "
        " --gtf {input.gtf} "
        " --outFileCounts {output.counts_gct_file} "
        " --outFileTPM {output.tpm_gct_file} "

rule VCF_chr_list:
    input:
        VCF = VCFstem + ".vcf.gz",
        index = VCFstem + ".vcf.gz.tbi"
    output:
        "results/vcf_chr_list.txt"
    shell:
        "tabix -l {input.VCF} > {output}"

# script from GTEX pipeline
# performs leafcutter clustering
# filtering of junctions by missingness and variability
# runs leafcutter_prepare_phenotype.py 
rule prepareSplicing:
    input:
        # expecting gzipped junction files with extension {sample}.junc.gz
        sample_key = sample_key,
        junction_file_list = junctionFileList, # from config - a file listing full paths to each junction file 
        exon_list = GTFexons, # hg38 exons from gencode v30 with gene_id and gene_name
        genes_gtf = GTF # full GTF or just gene starts and ends?
    output:
        counts = prefix + "_perind.counts.gz",
        counts_numers = prefix + "_perind_numers.counts.gz",
        clusters_pooled = prefix + "_pooled.gz",
        clusters_refined = prefix + "_refined.gz",
#       clusters_to_genes = prefix + ".leafcutter.clusters_to_genes.txt",
        phenotype_groups = prefix + ".leafcutter.phenotype_groups.txt",
        leafcutter_bed = prefix + ".leafcutter.bed.gz",
        leafcutter_bed_index = prefix + ".leafcutter.bed.gz.tbi",
        leafcutter_pcs = prefix + ".leafcutter.PCs.txt"
    params: 
        leafcutter_dir = "scripts/leafcutter/", # all leafcutter scripts hosted in a folder - some had to be converted py2 -> py3
        script = "scripts/sqtl_prepare_splicing.py",
        min_clu_reads = 30,
        min_clu_ratio = 0.001,
        max_intron_len = 500000,
        num_pcs = 10 # must be at least the number of samples!
    shell:  
        "ml R/3.6.0;"
        "ml tabix;"
        "python {params.script} "
                " {input.junction_file_list} "
                " {input.exon_list} "
                " {input.genes_gtf} "
                " {prefix} "
                " {input.sample_key} "
        " --min_clu_reads {params.min_clu_reads} "
                " --min_clu_ratio {params.min_clu_ratio} "
                " --max_intron_len {params.max_intron_len} "
                " --num_pcs {params.num_pcs} " 
        " --leafcutter_dir {params.leafcutter_dir}; "
        "rm *sorted.gz; " # clean up directory
        "rm {prefix}_perind.counts.filtered.gz.phen_* "

# from Francois Auguet at GTEX
# eQTLs - first prepare expression from RSEM count matrix
# then run PEER
# then add PEER factors to known covariates
rule prepareExpression:
    input:
        vcf_chr_list = "results/vcf_chr_list.txt",
        tpm_gct = prefix + "_tpm.gct",
        counts_gct =  prefix + "_counts.gct",
        gtf = "results/collapsed.gtf",
        sample_key = sample_key
    output:
        prefix + ".expression.bed.gz"
    params:
        script = "scripts/eqtl_prepare_expression.py"
    shell:
        "module load python/3.7.3;"
        " python {params.script} {input.tpm_gct} {input.counts_gct} {input.gtf} "
        " {input.sample_key} {input.vcf_chr_list} {prefix} "
        " --tpm_threshold 0.1 "
        " --count_threshold 6 "
        " --sample_frac_threshold 0.2 "
        " --normalization_method tmm "

       
rule runPEER:
    input:
        phenotype_matrix # either expression or splicing counts
        #prefix + ".expression.bed.gz"
    params:
        script = "scripts/run_PEER.R",
        num_peer = "{PEER_N}"
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt"
    run:
        if int(wildcards.PEER_N) > 0:
            shell("ml R/3.6.0; ")
            shell("Rscript {params.script} {input} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} {params.num_peer}")
        else:
            shell("touch {output}")

rule combineCovariates:
    input:
        pheno = prefix + ".phenotype.tensorQTL.{group_by}.bed.gz",
        geno =  genotypePCs,
        peer =  outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt",
        covariates = covariateFile
    output:
        cov_df = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.{group_by}.combined_covariates.txt"
    params:
        num_peer = "{PEER_N}",
        script = "scripts/combine_covariates.py",
        logNomFolder = outFolder + "peer{PEER_N}/logNomFolder",
        logPerFolder = outFolder + "peer{PEER_N}/logPerFolder"
    run:
        if int(wildcards.PEER_N) > 0:
            peerFile = "--add_covariates " + input.peer
        else:
            peerFile = ""
        shell("python {params.script} {peerFile} \
             --genotype_pcs {input.geno} {input.covariates} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}.{wildcards.group_by} ")
        # make sure combined covariate file has column names in same order as phenotype file
        phenotype_df = pd.read_csv(input.pheno, sep='\t', dtype={'#chr':str, '#Chr':str})
        covariate_df = pd.read_csv(output.cov_df, sep = "\t" )
        covariate_df = covariate_df[ ["ID"] + list(phenotype_df.columns[4:]) ]
        # write out
        covariate_df.to_csv(output.cov_df, header = True, index = False, sep = "\t")


# deprecated - may be able to reuse with tensorQTL but permutation step already performs qvalue testing
#rule summariseQTLtoolsResults:
#    input:
        #nominal_files = expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk{CHUNK}.nominals.txt", CHUNK = chunk_range, allow_missing=True),
        #permutation_files = expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk{CHUNK}.permutations.txt", CHUNK = chunk_range, allow_missing=True)
#    output:
#        full_nom = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.nominal.full.txt.gz",
#        full_perm = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.permutation.full.txt.gz",
#        sig = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt"
#    params:
#        nominal_wildcard = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk*nominals.txt",
#        permutation_wildcard  = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk*permutations.txt", 
#        script = "scripts/runFDR_cis.R",
#        file_prefix = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes",
#        # cat together all nominal and permutation files
#    shell:
#        "ml R/3.6.0;"
#        "cat {params.nominal_wildcard} | gzip -c > {output.full_nom};"
#        "cat {params.permutation_wildcard} | gzip -c > {output.full_perm};"
#        "Rscript {params.script} {output.full_perm} 0.05 {params.file_prefix};"


## TENSORQTL -----------------------------------------------------------------------

# tensorQTL expects fastQTL formatted phenotype files
# these don't include group info
# therefore take the phenotype file from QTLtools and strip out group info into a separate file
# this group info table should be two columns - phenotype and group
rule preparePhenotypeForTensorQTL:
    input:
        pheno = phenotype_matrix
    output:
        pheno = prefix + ".phenotype.tensorQTL.{group_by}.bed.gz",
        group = prefix + ".group.tensorQTL.{group_by}.bed.gz"
    run:
        # pandas - read in expression bed 
        # drop group column
        # write out table
        # create separate table with just phenotype and group ID - look into
        phenotype_df = pd.read_csv(input.pheno, sep='\t', dtype={'#chr':str, '#Chr':str})
        gene_id = phenotype_df.columns[3]
        # this should be the "geneid" for eQTls
        # for splicing QTLs either gene or cluster
        group_id = phenotype_df.columns[4]
        # sort by group ID
        if( wildcards.group_by == "cluster" ):
            # split phenotype_id by ":"
            split_df = phenotype_df[gene_id].str.split(":", n = 4, expand = True)
            # create new group ID - gene:clusterID
            phenotype_df[group_id] = split_df[4] + ":" + split_df[3]
        
        phenotype_df = phenotype_df.sort_values(by = group_id)
        # use column numbers not names - leafcutter table doesn't use same column names but positions are the same
        
        # create group and ID table - phenotype group (used for splicing)
        phenotype_group_df = phenotype_df[[gene_id, group_id]]

        # drop group and strand from phenotype table
        phenotype_df.drop(['strand', group_id], axis=1, inplace=True)
        phenotype_df.to_csv(output.pheno, index = False, sep = "\t")
        # no header!
        phenotype_group_df.to_csv(output.group, header = False, index = False, sep = "\t")

# extract the needed samples from the VCF 
rule getParticipants:
    input:
        txt = sample_key
    output:
        txt = prefix + "_participants.txt"
    run:
        sk = pd.read_csv(input.txt, sep = "\t")
        participants = sk[["participant_id"]]
        participants.to_csv(output.txt, index = False, header = False, sep = "\t")

# tensorQTL requires genotypes in PLINK format
# convert using plink2
# test VCF has multiallelic SNPs, hence the forcing with max-alleles
# in real data we've previously excluded multiallelics
# should go back to at some point
rule VCFtoPLINK:
    input:
        vcf = VCFstem + ".vcf.gz",
        participants = prefix + "_participants.txt"
    output:
        prefix + "_genotypes.fam"
    params:
        stem = prefix + "_genotypes"
    shell:
        "ml plink2; "
        "plink2 --make-bed "
        "--output-chr chrM "
        "--max-alleles 2 "
        "--keep {input.participants} "
        "--maf 0.01 "
        "--allow-extra-chr "
        "--max-maf 0.9975 "
        "--vcf {input.vcf} "
        "--out {params.stem} "

rule tensorQTL_cis:
    input:
        genotypes = prefix + "_genotypes.fam",
        phenotypes = prefix + ".phenotype.tensorQTL.{group_by}.bed.gz",
        covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.{group_by}.combined_covariates.txt",
        groups = prefix + ".group.tensorQTL.{group_by}.bed.gz"
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_{group_by}.cis_qtl.txt.gz"
    params:
        stem = prefix + "_genotypes",
        num_peer = "{PEER_N}",
        group = "{group_by}",
        group_string = group_string
    run:
        if interaction is False:
            shell( "python3 -m tensorqtl {params.stem} {input.phenotypes} \
             {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}_{params.group} \
             --covariates {input.covariates} \
             --mode cis ")
        if interaction is True:
            # the same 
            shell( "python3 -m tensorqtl {params.stem} {input.phenotypes} \
             {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}_{params.group} \
             --covariates {input.covariates} \
             --mode cis \
             --interaction {interaction_file}  --maf_threshold_interaction 0.05 ; \
             Rscript {params.script} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}_{params.group}.cis_qtl_top_assoc.txt.gz  {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}_{params.group}.cis_qtl.txt.gz"
        )

rule tensorQTL_cis_nominal:
    input:
        #perm_res = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.cis_qtl.txt.gz",
        genotypes = prefix + "_genotypes.fam",
        phenotypes = prefix + ".phenotype.tensorQTL.{group_by}.bed.gz",
        covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.{group_by}.combined_covariates.txt"
    output:
        expand( outFolder + "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}.cis_qtl_pairs.chr{CHROM}.parquet", CHROM = list(range(1,23)),  allow_missing=True )
    params:
        group = "{group_by}",
        stem = prefix + "_genotypes",
        num_peer = "{PEER_N}",
    run:
        if interaction is False:
            shell( "python3 -m tensorqtl {params.stem} {input.phenotypes} \
             {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}_{params.group} \
             --covariates {input.covariates} \
             --mode cis_nominal ")
        if interaction is True:
            # the same 
            shell( "python3 -m tensorqtl {params.stem} {input.phenotypes} \
             {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer}_{params.group} \
             --covariates {input.covariates} \
             --mode cis_nominal \
             --interaction {interaction_file}  --maf_threshold_interaction 0.05")
             
rule mergeNominalResult:
    input:
        expand( outFolder + "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}.cis_qtl_pairs.chr{CHROM}.parquet", CHROM = list(range(1,23)),  allow_missing=True )
    output:
        outFolder + "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}.cis_qtl_nominal_tabixed.tsv.gz",
        outFolder + "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}.cis_qtl_nominal_tabixed.tsv.gz.tbi", 
    params:
        prefix = "peer{PEER_N}/" + dataCode +"_peer{PEER_N}_{group_by}",
        script = "scripts/merge_nominal_results.R"
    
    shell:
        " ml R/3.6.0; ml arrow;"
        " Rscript {params.script} --vcf {VCF} --out_folder {outFolder} --prefix {params.prefix} "    

## MBV - MATCH BAM TO VARIANTS ---------------------------------------------------------------
# THIS SHOULD BE RUN BEFORE QTL MAPPING TO CHECK FOR SAMPLE SWAPS
# do on Chromosome 1 - big enough to properly identify sample swaps, hopefully

# index each bam file if not already
rule indexBam:
    input:
        bamFolder + "{sample}.bam"
    output:
        bamFolder + "{sample}.bam.bai"
    shell:
        "ml samtools/1.9;"
        "samtools index {input}"

# run match BAM to variants
rule matchBAM2VCF:
    input:
        bam = bamFolder + "{sample}.bam",
        vcf = VCF
    output:
        outFolder + "output/{sample}.bamstat.txt"
    shell:
        "ml qtltools/1.2;"
        "QTLtools mbv --bam {input.bam} --vcf {input.vcf} --filter-mapping-quality 150 --out {output}"

rule summariseResults:
    input:
        files = expand(outFolder + "output/{sample}.bamstat.txt", sample = BAM_SAMPLES)
    output:
        prefix + "_mbv_summary.txt"
    shell:
        "set +o pipefail;"
        "for i in {input.files};"
        "do cat $i | sort -k9nr,10nr | head -1 | awk -v i=$i \'{{print i, $0}}\'  ;"
        "done > {output};"


