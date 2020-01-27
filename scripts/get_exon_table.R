# Jack Humphrey

# reads in gtf with  rtracklayer
# finds exons
# writes out table
suppressMessages(library(rtracklayer, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE ))
suppressMessages(library(readr, quietly=TRUE ))
suppressMessages(library(optparse, quietly=TRUE ))

# get args
option_list <- list(
    make_option(c('--gtf'), help='the GTF file being used', default = "example"),
    make_option(c('--outFolder'), help='', default = "results/test/")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

gtf_file <- opt$gtf
outFolder <- opt$outFolder

#gtf_file <- "test/test.gtf"

gtf_id <- basename(gtf_file)
outFileName <- paste0(gtf_id, ".exons.txt.gz")
outFile <- file.path(outFolder, outFileName)

message( paste0(" * Reading GTF file from ", gtf_file ) )
gtf <- import.gff(gtf_file)

exons <- gtf[ gtf$type == "exon" ]

#Must have columns: chr, start, end, strand, gene_id[, gene_name]

exons_table <- data.frame(
	chr = seqnames(exons),
	start = start(exons), 
	end = end(exons), 
	strand = strand(exons), 
	gene_id = exons$gene_id, 
	gene_name = exons$gene_name,
	stringsAsFactors = FALSE 
)

exons_table <- distinct(exons_table)

message(paste0(" * Writing exon table to ", outFile) ) 
# readr::write_tsv automatically detects gzip extension and compresses output
write_tsv(exons_table, path = outFile)
#write.table(exons_table, out_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
