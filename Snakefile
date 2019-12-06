
#VCF = "test/test_all_chr.vcf.gz"
#GTF = "test/test.gtf" # cannot be gzipped
#dataCode = "test"
#sampleKey  = "test/test_sample_key.txt" # has to be "sample_id", "participant_id"
#genotypePCs = "test/test_genotype_PCs.txt" # rows are PCs, columns are samples
#countMatrixRData = 

VCF = config["VCF"]
GTF = config["GTF"]
dataCode = config["dataCode"]
sampleKey = config["sampleKey"]
genotypePCs = config["genotypePCs"]
countMatrixRData = config["countMatrixRData"]
# covariate file?

# derived variables
outFolder = "results/" + dataCode + "/"
prefix = outFolder + dataCode

# these will be created
counts_gct_file = "test/test_counts.gct"
tpm_gct_file = "test/test_tpm.gct"


# hardcoded variables
PEER_values = [20] # a list so can have a range of different values
chunk_number = 2 # at least as many chunks as there are chromosomes
chunk_range = range(1,chunk_number + 1)

QTLtools = "/hpc/packages/minerva-centos7/qtltools/1.2/bin/QTLtools"

shell.prefix('export PS1="";source activate QTL-pipeline; module load qtltools/1.2;')

rule all:
	input: 
		expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt", PEER_N = PEER_values),
		expand(outFolder + "peer{PEER_N}/" + dataCode +'_peer{PEER_N}_chunk{CHUNK}.nominals.txt', PEER_N = PEER_values, CHUNK = chunk_range)

rule collapseGTF:
	input: 
		GTF
	output:
		outFolder + "collapsed.gtf"
	params: 
		script = "scripts/collapse_annotation.py"
	shell: 
		"python {params.script} {input} {output} "
		
rule createGCTFiles:
	input:
		counts = countMatrixRData,
		key = sampleKey,
		gtf = outFolder + "collapsed.gtf"
	output:
		counts_gct_file = prefix + "_counts.gct",
		tpm_gct_file = prefix + "_tpm.gct"
	script:
		"scripts/create_GCT_files.R"

rule VCF_chr_list:
	input: 
		VCF
	output:
		outFolder + "vcf_chr_list.txt"
	shell:
		"for i in {input}; do tabix -l $i; done > {output}"

rule prepareExpression:
	input:
		vcf_chr_list = outFolder + "vcf_chr_list.txt",
		tpm_gct = prefix + "_tpm.gct",
		counts_gct =  prefix + "_counts.gct",
		gtf = outFolder + "collapsed.gtf",
		sampleKey = sampleKey
	output:
		prefix + ".expression.bed.gz"
	params:
		script = "scripts/eqtl_prepare_expression.py"
	shell:
		" python {params.script} {input.tpm_gct} {input.counts_gct} {input.gtf} "
		" {input.sampleKey} {input.vcf_chr_list} {prefix} "
		" --tpm_threshold 0.1 "
		" --count_threshold 6 "
		" --sample_frac_threshold 0.2 "
		" --normalization_method tmm "


rule runPEER:
	input: 
		prefix + ".expression.bed.gz"
	params:
		script = "scripts/run_PEER.R",
		num_peer = "{PEER_N}"
	output:
		outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt"
	shell:
		"ml R/3.6.0; "
		"Rscript {params.script} {input} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} {params.num_peer}"

rule combineCovariates:
	input:
		geno =	genotypePCs,
		peer =	outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt" 
	output:
		outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt"
	params: 
		num_peer = "{PEER_N}",
		script = "scripts/combine_covariates.py"
	shell:
		"python {params.script} {input.peer} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} "
    		" --genotype_pcs {input.geno} "
#    		" --add_covariates {add_covariates} "

rule QTLtools_nominal:
	input:
		expression = prefix + ".expression.bed.gz",
		covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt",
		vcf = VCF
	output:
		outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk{CHUNK}.nominals.txt"
	params:
		pval_threshold = 0.01,
		chunk_num = "{CHUNK}",
		chunk_max = chunk_number
	shell:
		"{QTLtools} cis --vcf {input.vcf} --bed {input.expression} --cov {input.covariates} "
		" --nominal {params.pval_threshold} "
		" --out {output}"
		" --chunk {params.chunk_num} {params.chunk_max} "
		" --normal "

rule QTLtools_permutation:
        input:
                expression = prefix + ".expression.bed.gz",
                covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt",
                vcf = VCF
        output:
                outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}_chunk{CHUNK}.permutations.txt"
        params:
                permutations = 1000,
		chunk_num = "{CHUNK}",
		chunk_max = chunk_number
        shell:
                "{QTLtools} cis --vcf {input.vcf} --bed {input.expression} --cov {input.covariates} "
                " --permute {params.permutations} "
                " --out {output}"
		" --chunk {params.chunk_num} {params.chunk_max} "
		" --normal "

rule summariseResults:
	input:
		expand(outFolder + "peer{PEER_N}/" + dataCode + '_peer{PEER_N}_chunk{CHUNK}.permutations.txt', PEER_N = PEER_values, CHUNK = chunk_range)
	output:
		full = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.full.txt.gz",
		sig = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt"
	params:
		files = outFolder + "peer{PEER_N}/" + dataCode + "*" + "_peer{PEER_N}*permutations.txt",
		script = "scripts/runFDR_cis.R",
		file_prefix = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes"
	shell:
		"ml R/3.6.0;"
		"cat {params.files} | gzip -c > {output.full};"
		"Rscript {params.script} {output.full} 0.05 {params.file_prefix}"
