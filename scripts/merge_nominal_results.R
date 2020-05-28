# merge together all nominal associations into a single file
# match in chr and pos from VCF
# sort, bgzip and tabix for random access

library(tidyverse)
library(arrow)
library(optparse)

option_list <- list(
    make_option(c('--vcf', '-v' ), help='The full path to the VCF used in the QTL analysis', default = ""),
    make_option(c('-o', '--out_folder'), help = "the full path to where the input parquet files are and the output should be written", default = ""),
    make_option(c('--prefix', '-p'), help = "the dataset name used") 
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

vcf <- opt$vcf
prefix <- opt$prefix
out_folder <- opt$out_folder


#VCF <- "/sc/arion/projects/als-omics/WGS_QC/NYGC_Freeze02_European_Feb2020/WGS_QC_Pipeline/NYGC_Freeze02_European_Feb2020/output/chrAll_QCFinished_MAF0.01.anno.vcf.gz"

#out_folder <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/Hippocampus_expression/peer15/"
#prefix <- "Hippocampus_expression_peer15_gene"


outFile <- paste0( out_folder, prefix, ".cis_qtl_nominal_tabixed.tsv" )

stopifnot(file.exists(vcf) )

vcf_df <- tempfile()

# read in SNP positions from VCF
cmd <- paste0("ml bcftools; bcftools query -f '%CHROM\t%POS\t%ID\n' ", vcf," > ", vcf_df )
message(" * extracting SNP IDs from VCF")
system(cmd)

message(" * reading in SNP positions from VCF")
snp_df <- read_tsv(vcf_df, col_names = c("chr", "pos", "variant_id"), col_types = c("cnc") )

system(paste0("rm ", vcf_df ))


# find all parquet files for a QTL analysis
all_parquet <- paste0(out_folder, prefix, ".cis_qtl_pairs.chr",1:22,".parquet")

for( file in all_parquet){
    if( !file.exists(file) ){
        stop( file, "does not exist" )
    }
}
stopifnot( all(file.exists(all_parquet) ))

# read in each file with read_parquet
# join into one huge table
all_res <- 
    map_df( all_parquet, ~{
        message( " * reading in ", .x )
        df <- read_parquet(.x)
        df <- dplyr::left_join(df, snp_df, by = "variant_id")
        # sort by position
        df <- df[order(df$pos), ]
        # any weirdness - remove
        df <- df[!is.na(df$pos) ,]
    return(df)
})

message(" * writing to ", outFile)
write_tsv(all_res, path = outFile)

# bgzip
bgzip_cmd <- paste0(" ml bcftools; bgzip -f ", outFile)
message( " * ", bgzip_cmd)
system(bgzip_cmd)

# tabix
tabix_cmd <- paste0(" ml bcftools; tabix -S 1 -s 10 -b 11 -e 11 ", outFile, ".gz" )
message(" * ", tabix_cmd)
system(tabix_cmd)

