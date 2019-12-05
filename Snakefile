VCF = "test/test_chr1.vcf.gz"
GTF = "test/test_chr1.gtf" # cannot be gzipped
dataCode = "test"
counts_gct_file = "test/test_counts.gct"
tpm_gct_file = "test/test_tpm.gct"
sample_lookup_file = "test/NYGC_ALS_CervicalSpinalCord_sample_key_columns_reversed.txt"

outFolder = "results/" + dataCode + "/"

prefix = outFolder + dataCode

PEER_values = [1, 5, 10, 15]

rule all:
	input: expand(prefix + '_{PEER_N}.PEER_covariates.txt',  PEER_N = PEER_values)  

rule collapseGTF:
	input: 
		GTF
	output:
		outFolder + "collapsed.gtf"
	params: 
		script = "scripts/collapse_annotation.py"
	shell: 
		"python {params.script} {input} {output} "
		

rule VCF_chr_list:
	input: 
		VCF
	output:
		outFolder + "vcf_chr_list.txt"
	shell:
		"ml tabix; tabix -l {input} > {output}"

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
		#" ml python/3.7.3;"
		" python {params.script} {input.tpm_gct} {input.counts_gct} {input.gtf} "
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
		prefix + "_{PEER_N}.PEER_covariates.txt"
	shell:
		"ml R/3.6.0; "
		"Rscript {params.script} {input} {prefix}_{params.num_peer} {params.num_peer}"

