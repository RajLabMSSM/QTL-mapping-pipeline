import os
import glob

#SAMPLES = ["HG00381"]
#outFolder = "test/match_bam_to_genotypes/"
#VCF = outFolder + "genotypes.chr22.vcf.gz"
shell.prefix('export PS1="";source activate QTL-pipeline; module load qtltools/1.2;')

VCF = config["VCF"] # full path to gzipped, tabixed VCF of same chromosome as CHR
inFolder = config["inFolder"]
outFolder = config["outFolder"]

bamSuffix = config["bamSuffix"]
CHR = config["CHR"]


# find sample names from BAM files in inFolder
SAMPLES = [os.path.basename(x).strip(bamSuffix) for x in glob.glob(inFolder + "/*" + bamSuffix)]

print(SAMPLES)

rule all:
	input:
		outFolder + "summary.txt"

# for each bam file 
# extract single chr - 22

rule indexBam:
	input:
		inFolder + "{sample}.bam"
	output:
		inFolder + "{sample}.bam.bai"
	shell:
		"ml samtools/1.9;"
		"samtools index {input}"

rule extractBAMChr:
	input:
		bam = inFolder + "{sample}.bam",
		bai = inFolder + "{sample}.bam.bai"
	output:
		bam = outFolder + "{sample}." + CHR + ".bam",
		bai = outFolder + "{sample}." + CHR + ".bam.bai"
	shell:
		"ml samtools/1.9;"
		"samtools view -bh {input.bam} {CHR} > {output.bam};"
		"samtools index {output.bam} "

rule matchBAM2VCF:
	input:
		bam = outFolder + "{sample}." + CHR + ".bam",
		vcf = VCF
	output:
		outFolder + "{sample}." + CHR + ".bamstat.txt"
	shell:
		#"ml qtltools/1.2;"
		"QTLtools mbv --bam {input.bam} --vcf {input.vcf} --filter-mapping-quality 150 --out {output}"

rule summariseResults:
	input:
		files = expand(outFolder + "{sample}." + CHR + ".bamstat.txt", sample = SAMPLES)
	output:
		outFolder + "summary.txt"
	shell:
		"set +o pipefail;"
		"for i in {input.files};"
		"do cat $i | sort -k9nr,10nr | head -1 | awk -v i=$i \'{{print i, $0}}\'  ;"
		"done > {output};"
