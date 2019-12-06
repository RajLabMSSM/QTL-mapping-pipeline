library(readr)
library(dplyr,warn.conflicts = FALSE)

# using the provided gene matrix and TPM matrix
# extract just the genes in the test GTF file
# write to a GCT format file
# which has the stupid header and columns

message(" * loading count matrix")
load("/sc/orga/projects/als-omics/NYGC_ALS/data/oct_2019_gene_matrix.RData")

message(" * loading GTF")
gtf <- rtracklayer::readGFF("test.gtf")

samples <- read_tsv("test_sample_key.txt", col_types = "cc")

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

writeLines(gct_header, con = "test_counts.gct")
writeLines(gct_header, con = "test_tpm.gct" )

write_tsv(counts_gct, "test_counts.gct", append = TRUE, col_names = TRUE)
write_tsv(tpm_gct, "test_tpm.gct", append = TRUE, col_names = TRUE)
