library(readr)
library(dplyr,warn.conflicts = FALSE)

# using the provided gene matrix and TPM matrix
# extract just the genes in the test GTF file
# write to a GCT format file
# which has the stupid header and columns

library(optparse)

option_list <- list(
    make_option(c('--counts'), help = 'RData count matrix and TPM matrix', default = ""),
    make_option(c('--gtf'), help = 'GENCODE GTF file', default = ""),
    make_option(c('--key'), help = 'sample key', default = ""),
    make_option(c('--outFileCounts'), help = "output a counts.gct file", default = ""),
    make_option(c('--outFileTPM'), help = "output a tpm.gct file", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


countMatrixRData <- opt$counts
gtf <- opt$gtf
sample_key <- opt$key

counts_gct_file <- opt$outFileCounts
tpm_gct_file <- opt$outFileTPM

#countMatrixRData <- snakemake@input[["counts"]]
#gtf <- snakemake@input[["gtf"]]
#sample_key <- snakemake@input[["key"]]

# output files
#counts_gct_file <- snakemake@output[["counts_gct_file"]]
#tpm_gct_file <- snakemake@output[["tpm_gct_file"]]

message(" * loading count matrix")
load(countMatrixRData)

message(" * loading GTF")
gtf <- rtracklayer::readGFF(gtf)

samples <- read_tsv(sample_key, col_types = "cc")

genes <- unique( gtf$gene_id )

counts <- genes_counts[ genes, samples$sample_id]
tpm <- genes_tpm[ genes, samples$sample_id]

# make first 2 columns for GCT 

gct_columns <- data.frame( 
	Name = row.names(counts),
	Description = row.names(counts),
	stringsAsFactors = FALSE 
	)

# first 2 rows are pointless but have to be in there 
# pipeline code skips them anyway
# writeLines them first then append the table after with write_tsv (append = TRUE)
gct_header <- c("#1.2", paste(dim(counts), collapse = "\t") )

message(" * writing GCT files")
counts_gct <- cbind(gct_columns, counts)
tpm_gct <- cbind(gct_columns, tpm)

writeLines(gct_header, con = counts_gct_file)
writeLines(gct_header, con = tpm_gct_file)

write_tsv(counts_gct, counts_gct_file, append = TRUE, col_names = TRUE)
write_tsv(tpm_gct, tpm_gct_file, append = TRUE, col_names = TRUE)
