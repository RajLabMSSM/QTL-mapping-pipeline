

CHR = "chr22"

rule all:
	input: 
		"summary.txt"

SAMPLES = ["HG00381"]
outFolder = "test/match_bam_to_genotypes/"

VCF = outFolder + "genotypes.chr22.vcf.gz"

rule all:
	input:
		outFolder + "summary.txt"

# for each bam file 
# extract single chr - 22
rule extractBAMChr:
	input:
		outFolder + "{sample}.bam"
	output:
		outFolder + "{sample}." + CHR + ".bam",
		outFolder + "{sample}." + CHR + ".bam.bai"
	shell:
		"ml samtools/1.9;"
		"samtools view -bh {input} {CHR} > {output}"
		"samtools index {input} "

rule matchBAM2VCF:
	input:
		bam = outFolder + "{sample}." + CHR + ".bam",
		vcf = VCF
	output:
		outFolder + "{sample}." + CHR + "bamstat.txt"
	shell:
		"ml qtltools/1.2;"
		"QTLtools mbv --bam {input.bam} --vcf {input.vcf} --filter-mapping-quality 150 --out {output}"

rule summariseResults:
	input:
		expand(outFolder + "{sample}." + CHR + "bamstat.txt", sample = SAMPLES)
	output:
		outFolder + "summary.txt"
	shell:
		"for i in {input}; do;"
		"cat $i | sort -k9nr,10nr | head -1 | awk -v i=$i '\{print i, $0\}'  ;"
		"done > {output}"
