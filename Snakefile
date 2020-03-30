# QTL mapping pipeline
# Jack Humphrey 
import glob
import os

leafcutter_dir = "/sc/orga/projects/ad-omics/data/software/leafcutter/"
GTF = "/sc/orga/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.annotation.gtf" # cannot be gzipped
GTFexons = GTF + ".exons.txt.gz" 
QTLtools = "/hpc/packages/minerva-centos7/qtltools/1.2/bin/QTLtools"

nPerm = 10000 # number of permutations of the permutation pass
chunk_number = 22 * 40 # at least as many chunks as there are chromosomes
chunk_range = range(1,chunk_number + 1)

shell.prefix('export PS1="";source activate snakemake; ml qtltools/1.2; ml R/3.6.0;')

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

if(mode == "mbv"):
    bamFolder = config["bamFolder"]
    bamSuffix = config["bamSuffix"]
    BAM_SAMPLES = [os.path.basename(x).strip(bamSuffix) for x in glob.glob(bamFolder + "/*" + bamSuffix)]
    dataCode = dataCode + "_mbv"
    outFolder = "results/" + dataCode + "/"
    prefix = outFolder + dataCode
    final_output = prefix + "_mbv_summary.txt"
else:
    sample_key = config["sampleKey"]
    genotypePCs = config["genotypePCs"]
    covariateFile = config["covariateFile"]
    PEER_values = config["PEER_values"]

if(mode == "eQTL"):
    PEER_values = [10]
    dataCode = dataCode + "_expression"
    outFolder = "results/" + dataCode + "/"
    prefix = outFolder + dataCode
    countMatrixRData = config["countMatrixRData"]
    phenotype_matrix = prefix + ".expression.bed.gz"
    internal_covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt"
    final_output = expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt", PEER_N = PEER_values)
    grouping_param = ""
if(mode == "sQTL"):
    PEER_values = [15]
    dataCode = dataCode + "_splicing"
    junctionFileList = config["junctionFileList"]
    outFolder = "results/" + dataCode + "/"
    prefix = outFolder + dataCode
    phenotype_matrix = prefix + ".leafcutter.bed.gz" 
    internal_covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt" #outFolder + ".leafcutter.PCs.txt"
    grouping_param =" --grp-best "
    final_output = expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt", PEER_N = PEER_values) #prefix + "_results.genes.significant.txt" 

# if interaction requested then use TensorQTL and include script that matches interaction values to covariates and samples

print(" * QTL-mapping pipeline *")
print(" Jack Humphrey 2019-2020 ")
print(" * Data code is : %s " % dataCode)
print(" * Mode selected is: %s" % mode)


rule all:
    input:
        final_output
        #prefix + ".leafcutter.PCs.txt"
        #expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + ".cis_qtl.txt.gz", PEER_N = PEER_values),
        #expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + ".cis_nominal_qtl.txt.gz", PEER_N = PEER_values)
        #expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt", PEER_N = PEER_values),
        #expand(outFolder + "peer{PEER_N}/" + dataCode +'_peer{PEER_N}_chunk{CHUNK}.nominals.txt', PEER_N = PEER_values, CHUNK = chunk_range),
        #expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.nominal.full.txt.gz", PEER_N = PEER_values)

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
    shell:
        "ml R/3.6.0; "
        "Rscript {params.script} {input} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} {params.num_peer}"

if covariateFile != "":
    covariate_string = " --add_covariates " + covariateFile 
else:
    covariate_string = ""

rule combineCovariates:
    input:
        geno =  genotypePCs,
        peer =  outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt"
        #covariates = covariateFile
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt"
    params:
        covariates = covariate_string,
        num_peer = "{PEER_N}",
        script = "scripts/combine_covariates.py",
        logNomFolder = outFolder + "peer{PEER_N}/logNomFolder",
        logPerFolder = outFolder + "peer{PEER_N}/logPerFolder"
    shell:
        "python {params.script} {input.peer} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} "
            " --genotype_pcs {input.geno} "
            "{params.covariates} ;"
        "mkdir -p {params.logNomFolder};"
        "mkdir -p {params.logPerFolder};"

#Changes to QTLtools_nominal and QTLtools_permutation do the following
#implements a do while loop to keep repeating the initial qtltools command to remove phenotypes without enough variants
#If no output is created (happens when there is no genotype within the region), creates an empty file (last if statement in the shell script)
#Additionally, it seems that a single log file will account for the logs for both nominal and permutation. They should be the same, and they should remove the same phenotypes
rule QTLtools_nominal:
    input:
        phenotypes = phenotype_matrix,
        #expression = prefix + ".expression.bed.gz",
        covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt",
        vcf = VCFstem + ".vcf.gz"
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk{CHUNK}.nominals.txt"
    params:
        pval_threshold = 1,
        logNomFolder = outFolder + "peer{PEER_N}/logNomFolder",
        chunk_num = "{CHUNK}",
        chunk_max = chunk_number
    shell:
        "touch {params.logNomFolder}/Chunk{params.chunk_num}_exclude_phenotypes.txt;" #create phenotype exclusion file beforehand
        "success=false;" #condition for do while loop
        "until [ \"$success\" = true ]; "
        "do {{"
        " {QTLtools} cis --vcf {input.vcf} --bed {input.phenotypes} --cov {input.covariates} "
        " --nominal {params.pval_threshold} "
        "--out {output} --chunk {params.chunk_num} {params.chunk_max} "
        " --normal --exclude-phenotypes {params.logNomFolder}/Chunk{params.chunk_num}_exclude_phenotypes.txt " # > {params.logNomFolder}/Chunk{params.chunk_num}_log.txt"
        " && success=true; }}" #&& means the second command will only execute if the first executes without error, || responds to the error
        " || {{ grep 'Processing phenotype' {params.logNomFolder}/Chunk{params.chunk_num}_log.txt > {params.logNomFolder}/Chunk{params.chunk_num}_log2.txt;" #greps the last phenotype processed. This one brought up the error
        "sed -n \"$(wc -l < {params.logNomFolder}/Chunk{params.chunk_num}_log2.txt)p\" {params.logNomFolder}/Chunk{params.chunk_num}_log2.txt | cut -d'[' -f 2 | cut -d']' -f 1 " #performs string parsing on the statement around the phenotype
        ">> {params.logNomFolder}/Chunk{params.chunk_num}_exclude_phenotypes.txt; }}; done;" #adds this phenotype to a file containing all phenotypes to exclude for this chunk
        "if [[ ! -f {output} ]]; then touch {output}; fi;" #handles errors that simply exit the code without creating an output


rule QTLtools_permutation:
    input:
        phenotypes = phenotype_matrix,
        #expression = prefix + ".expression.bed.gz",
        covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt",
        vcf = VCFstem + ".vcf.gz"
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk{CHUNK}.permutations.txt"
    params:
        permutations = nPerm,
        logPerFolder = outFolder + "peer{PEER_N}/logPerFolder",
        chunk_num = "{CHUNK}",
        chunk_max = chunk_number
    shell:
        "touch {params.logPerFolder}/Chunk{params.chunk_num}_exclude_phenotypes.txt;" #create phenotype exclusion file beforehand
        "success=false;" #condition for do while loop
        "until [ \"$success\" = true ]; do {{" 
        " {QTLtools} cis --vcf {input.vcf} --bed {input.phenotypes} --cov {input.covariates} "
        " {grouping_param} " # for sQTLs - permute by gene
        " --permute {params.permutations} "
        "--out {output} --chunk {params.chunk_num} {params.chunk_max} "
        " --normal "  #--exclude-phenotypes {params.logPerFolder}/Chunk{params.chunk_num}_exclude_phenotypes.txt "
        " > {params.logPerFolder}/Chunk{params.chunk_num}_log.txt"
        " && success=true; }}" #&& means the second command will only execute if the first executes without error, || responds to the error
        " || {{ grep 'Processing phenotype' {params.logPerFolder}/Chunk{params.chunk_num}_log.txt > {params.logPerFolder}/Chunk{params.chunk_num}_log2.txt;" #greps the last phenotype processed. This one brought up the error
        "sed -n \"$(wc -l < {params.logPerFolder}/Chunk{params.chunk_num}_log2.txt)p\" {params.logPerFolder}/Chunk{params.chunk_num}_log2.txt | cut -d'[' -f 2 | cut -d']' -f 1 " #performs string parsing on the statement around the phenotype
        ">> {params.logPerFolder}/Chunk{params.chunk_num}_exclude_phenotypes.txt; }}; done;" #adds this phenotype to a file containing all phenotypes to exclude for this chunk
        "if [[ ! -f {output} ]]; then touch {output}; fi;" #handles errors that simply exit the code without creating an output

rule summariseQTLtoolsResults:
    input:
        nominal_files = expand(outFolder + "peer{PEER_N}/" + dataCode + '_peer{PEER_N}_chunk{CHUNK}.nominals.txt', CHUNK = chunk_range, allow_missing=True),
        permutation_files = expand(outFolder + "peer{PEER_N}/" + dataCode + '_peer{PEER_N}_chunk{CHUNK}.permutations.txt', CHUNK = chunk_range, allow_missing=True)
    output:
        full_nom = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.nominal.full.txt.gz",
        full_perm = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.permutation.full.txt.gz",
        sig = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt"
    params:
        nominal_wildcard = outFolder + "peer{PEER_N}/" + dataCode + '_peer{PEER_N}_chunk*nominals.txt',
        permutation_wildcard  = outFolder + "peer{PEER_N}/" + dataCode + '_peer{PEER_N}_chunk*permutations.txt', 
        script = "scripts/runFDR_cis.R",
        file_prefix = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes",
        # cat together all nominal and permutation files
    shell:
        "ml R/3.6.0;"
        "cat {params.nominal_wildcard} | gzip -c > {output.full_nom};"
        "cat {params.permutation_wildcard} | gzip -c > {output.full_perm};"
        "Rscript {params.script} {output.full_perm} 0.05 {params.file_prefix};"

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


        #summarizes removed chunks/genotypes
#        "touch {params.logNomFolder}/removed_genotypes_chunks.txt;"
#        "touch {params.logNomFolder}/removed_genotypes_loci.txt;"
#        "grep \"EXITED:\" {params.logNomFolder}/* | rev | cut -d\"/\" -f1 | rev | cut -d\"_\" -f1 > {params.logNomFolder}/removed_genotypes_chunks.txt;"
#        "grep -A1 \"VCF\" $(grep \"EXITED:\" {params.logNomFolder}/* | cut -d\":\" -f1) | awk 'NR % 3 == 2' | cut -d'[' -f2 | cut -d']' -f1 > {params.logNomFolder}/removed_genotypes_loci.txt;"        
        #"ls {params.logNomFolder}/ | grep \"EXITED:\" | cut -d\"_\" -f1 > {params.logNomFolder}/removed_genotypes_chunks.txt;"
        #"ls {params.logNomFolder}/ | grep \"EXITED:\" | cut -d\":\" -f1 | grep -A1 \"VCF\" | awk 'NR % 3 == 2' | cut -d'[' -f2 | cut -d']' -f1 > {params.logNomFolder}/removed_genotypes_loci.txt;"
#        "paste {params.logNomFolder}/removed_genotypes_chunks.txt {params.logNomFolder}/removed_genotypes_loci.txt > {params.logNomFolder}/Removed_Genotypes.txt;"
 #       "rm {params.logNomFolder}/removed_genotypes_chunks.txt;"
  #      "rm {params.logNomFolder}/removed_genotypes_loci.txt;"

        #summarizes removed _exclude_phenotypes
   #     "touch {params.logNomFolder}/removed_phenotypes_chunks.txt;"
    #    "touch {params.logNomFolder}/removed_phenotypes_loci.txt;"
     #   "touch {params.logNomFolder}/removed_phenotypes_list.txt;"
      #  "find {params.logNomFolder}/*exclude_phenotypes.txt -not -empty -ls | rev | cut -d'/' -f1 | rev | cut -d'_' -f1 > {params.logNomFolder}/removed_phenotypes_chunks.txt;"
       # "grep -A1 \"Reading phenotype data in\" {params.logPerFolder}/\"$(find {params.logNomFolder}/*exclude_phenotypes.txt -not -empty -ls | rev |cut -d'/' -f1 | rev | cut -d'_' -f1)\"_log.txt | awk 'NR % 3 == 2' | cut -d'[' -f2 | cut -d']' -f1 > {params.logNomFolder}/removed_phenotypes_loci.txt;"
        #"tr '\n' ',' < {params.logPerFolder}/$(find {params.logNomFolder}/*exclude_phenotypes.txt -not -empty -ls | rev | cut -d'/' -f1 | rev) > {params.logNomFolder}/removed_phenotypes_list.txt;"
        #"paste {params.logNomFolder}/removed_phenotypes_chunks.txt {params.logNomFolder}/removed_phenotypes_loci.txt {params.logNomFolder}/removed_phenotypes_list.txt > {params.logNomFolder}/Removed_Phenotypes.txt;"
        #"rm {params.logNomFolder}/removed_phenotypes_chunks.txt;"
        #"rm {params.logNomFolder}/removed_phenotypes_loci.txt;"
        #"rm {params.logNomFolder}/removed_phenotypes_list.txt;"

        #summarizes removed chunks/genotypes
        #"touch {params.logPerFolder}/removed_genotypes_chunks.txt;"
        #"touch {params.logPerFolder}/removed_genotypes_loci.txt;"
        #"grep \"EXITED:\" {params.logPerFolder}/* | rev | cut -d\"/\" -f1 | rev | cut -d\"_\" -f1 > {params.logPerFolder}/removed_genotypes_chunks.txt;"
       # "grep -A1 \"VCF\" $(grep \"EXITED:\" {params.logPerFolder}/* | cut -d\":\" -f1) | awk 'NR % 3 == 2' | cut -d'[' -f2 | cut -d']' -f1 > {params.logPerFolder}/removed_genotypes_loci.txt;"
        #"ls {params.logPerFolder}/ | grep \"EXITED:\" | cut -d\"_\" -f1 > {params.logPerFolder}/removed_genotypes_chunks.txt;"
        #"ls {params.logPerFolder}/ | grep \"EXITED:\" | cut -d\":\" -f1 | grep -A1 \"VCF\" | awk 'NR % 3 == 2' | cut -d'[' -f2 | cut -d']' -f1 > {params.logPerFolder}/removed_genotypes_loci.txt;"
       # "paste {params.logPerFolder}/removed_genotypes_chunks.txt {params.logPerFolder}/removed_genotypes_loci.txt > {params.logPerFolder}/Removed_Genotypes.txt;"
        #        "rm {params.logPerFolder}/removed_genotypes_chunks.txt;"
        #"rm {params.logPerFolder}/removed_genotypes_loci.txt;"

        #summarizes removed _exclude_phenotypes
        #"touch {params.logPerFolder}/removed_phenotypes_chunks.txt;"
        #"touch {params.logPerFolder}/removed_phenotypes_loci.txt;"
        #"touch {params.logPerFolder}/removed_phenotypes_list.txt;"
       # #"find {params.logPerFolder}/*exclude_phenotypes.txt -not -empty -ls | rev | cut -d'/' -f1 | rev | cut -d'_' -f1 > {params.logPerFolder}/removed_phenotypes_chunks.txt;"
        #"grep -A1 \"Reading phenotype data in\" {params.logPerFolder}/\"$(find {params.logPerFolder}/*exclude_phenotypes.txt -not -empty -ls | rev |cut -d'/' -f1 | rev | cut -d'_' -f1)\"_log.txt | awk 'NR % 3 == 2' | cut -d'[' -f2 | cut -d']' -f1 > {params.logPerFolder}/removed_phenotypes_loci.txt;"
        #"tr '\n' ',' < {params.logPerFolder}/$(find {params.logPerFolder}/*exclude_phenotypes.txt -not -empty -ls | rev | cut -d'/' -f1 | rev) > {params.logPerFolder}/removed_phenotypes_list.txt;"
       # "paste {params.logPerFolder}/removed_phenotypes_chunks.txt {params.logPerFolder}/removed_phenotypes_loci.txt {params.logPerFolder}/removed_phenotypes_list.txt > {params.logPerFolder}/Removed_Phenotypes.txt;"
       # "rm {params.logPerFolder}/removed_phenotypes_chunks.txt;"
       # "rm {params.logPerFolder}/removed_phenotypes_loci.txt;"

## TENSORQTL -----------------------------------------------------------------------

# tensorQTL requires genotypes in PLINK format
# convert using plink2
# test VCF has multiallelic SNPs, hence the forcing with max-alleles
# in real data we've previously excluded multiallelics
# should go back to at some point
rule VCFtoPLINK:
    input:
        VCFstem + ".vcf.gz"
    output:
        VCFstem + ".fam"
    shell:
        "ml plink2; "
        "plink2 --make-bed "
            "--output-chr chrM "
        "--max-alleles 2 "
            "--vcf {input} "
            "--out {VCFstem} "

rule tensorQTL_cis:
    input:
        VCFstem + ".fam",
        phenotypes = phenotype_matrix,
        #expression = prefix + ".expression.bed.gz",
                covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt"
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.cis_qtl.txt.gz"
    params:
        num_peer = "{PEER_N}"
    shell:
        "python3 -m tensorqtl {VCFstem} {input.phenotypes} "
        "{outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} " #" ${prefix} "
            " --covariates {input.covariates} "
            " --mode cis "

rule tensorQTL_cis_nominal:
    input:
        VCFstem + ".fam",
        phenotypes = phenotype_matrix,
        #expression = prefix + ".expression.bed.gz",
        covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt"
    output:
        outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.cis_nominal_qtl.txt.gz"
    params:
        num_peer = "{PEER_N}"
    shell:
        "python3 -m tensorqtl {VCFstem} {input.phenotypes} "
        "{outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} " #" ${prefix} "
        "--covariates {input.covariates} "
        " --mode cis_nominal "
        

        "rm {params.logPerFolder}/removed_phenotypes_list.txt;"
