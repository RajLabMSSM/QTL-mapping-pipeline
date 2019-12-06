VCF = "test/test_chr1.vcf.gz"
VCFstem = "test/test_chr" # + CHR + ".vcf.gz"
GTF = "test/test_chr1.gtf" # cannot be gzipped
dataCode = "test"
counts_gct_file = "test/test_counts.gct"
tpm_gct_file = "test/test_tpm.gct"
sample_lookup_file  = "test/test_sample_key.txt" # has to be "sample_id", "participant_id"
outFolder = "results/" + dataCode + "/"
genotype_PCs = "test/test_genotype_PCs.txt" # rows are PCs, columns are samples

prefix = outFolder + dataCode

PEER_values = [1,5,10]
chromosomes = [1]

chunk_number = 1

tabix = "/hpc/packages/minerva-centos7/htslib/1.9/bin/tabix"
#python = "/hpc/packages/minerva-centos7/python/3.7.3/bin/python"
python = "python"

shell.prefix('ml tabix; ml qtltools/1.2; ml python/3.7.3;')

chunk_range = range(1,chunk_number + 1)

rule all:
	input: 
		expand(outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}" + "_results.genes.significant.txt", PEER_N = PEER_values),
		#expand(prefix + '_chr{CHR}_peer{PEER_N}_chunk{CHUNK}.permutations.txt', PEER_N = PEER_values, CHR = chromosomes, CHUNK = chunk_range),
		expand(outFolder + "peer{PEER_N}/" + dataCode +'_chr{CHR}_peer{PEER_N}_chunk{CHUNK}.nominals.txt', PEER_N = PEER_values, CHR = chromosomes, CHUNK = chunk_range)

rule collapseGTF:
	input: 
		GTF
	output:
		outFolder + "collapsed.gtf"
	params: 
		script = "scripts/collapse_annotation.py"
	shell: 
		"{python} {params.script} {input} {output} "
		

rule VCF_chr_list:
	input: 
		VCFstem + "{CHR}.vcf.gz"
	output:
		outFolder + "vcf_chr_list.txt"
	shell:
		"{tabix} -l {input} > {output}"

rule prepareExpression:
	input:
		vcf_chr_list = outFolder + "vcf_chr_list.txt",
		tpm_gct = tpm_gct_file,
		counts_gct = counts_gct_file,
		gtf = outFolder + "collapsed.gtf",
		sample_lookup = sample_lookup_file
	output:
		prefix + ".expression.bed.gz"
	params:
		script = "scripts/eqtl_prepare_expression.py"
	shell:
		" {python} {params.script} {input.tpm_gct} {input.counts_gct} {input.gtf} "
		" {input.sample_lookup} {input.vcf_chr_list} {prefix} "
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
		geno =	genotype_PCs,
		peer =	outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.PEER_covariates.txt" 
	output:
		outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt"
	params: 
		num_peer = "{PEER_N}",
		script = "scripts/combine_covariates.py"
	shell:
		"{python} {params.script} {input.peer} {outFolder}peer{params.num_peer}/{dataCode}_peer{params.num_peer} "
    		" --genotype_pcs {input.geno} "
#    		" --add_covariates {add_covariates} "

rule QTLtools_nominal:
	input:
		expression = prefix + ".expression.bed.gz",
		covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt",
		vcf = VCFstem + "{CHR}.vcf.gz"
	output:
		outFolder + "peer{PEER_N}/" + dataCode + "_chr{CHR}_peer{PEER_N}_chunk{CHUNK}.nominals.txt"
	params:
		pval_threshold = 0.01,
		chunk_num = "{CHUNK}",
		chunk_max = chunk_number
	shell:
		"ml qtltools/1.2;"
		"QTLtools cis --vcf {input.vcf} --bed {input.expression} --cov {input.covariates} "
		" --nominal {params.pval_threshold} "
		" --out {output}"
		" --chunk {params.chunk_num} {params.chunk_max} "


rule QTLtools_permutation:
        input:
                expression = prefix + ".expression.bed.gz",
                covariates = outFolder + "peer{PEER_N}/" + dataCode + "_peer{PEER_N}.combined_covariates.txt",
                vcf = VCF
        output:
                outFolder + "peer{PEER_N}/" + dataCode + "_chr{CHR}_peer{PEER_N}_chunk{CHUNK}.permutations.txt"
        params:
                permutations = 1000,
                region = "chr{CHR}",
		chunk_num = "{CHUNK}",
		chunk_max = chunk_number
        shell:
                "ml qtltools/1.2;"
                "QTLtools cis --vcf {input.vcf} --bed {input.expression} --cov {input.covariates} "
                " --permute {params.permutations} "
                " --out {output}"
		" --chunk {params.chunk_num} {params.chunk_max} "


rule summariseResults:
	input:
		expand(outFolder + "peer{PEER_N}/" + dataCode + '_chr{CHR}_peer{PEER_N}_chunk{CHUNK}.permutations.txt', PEER_N = PEER_values, CHR = chromosomes, CHUNK = chunk_range)
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
