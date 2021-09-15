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

# for testing
#vcf <- "/sc/arion/projects/als-omics/WGS_QC/NYGC_Freeze02_European_Feb2020/WGS_QC_Pipeline/NYGC_Freeze02_European_Feb2020/output/chrAll_QCFinished_MAF0.01.anno.vcf.gz"
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
all_parquet <- paste0(out_folder, prefix, ".cis_qtl_pairs.chr",gsub("chr","",unique(snp_df$chr)),".parquet")

for( file in all_parquet){
    if( !file.exists(file) ){
        stop( file, "does not exist" )
    }
}
stopifnot( all(file.exists(all_parquet) ))

header_file <- paste0(outFile, "_header")

# read in each file with read_parquet
# join on SNP info 
# sort by position
# write out each section
all_res <- 
    walk( all_parquet, ~{
        message( " * reading in ", .x )
        df <- read_parquet(.x,  compression = "uncompressed")
        df <- dplyr::left_join(df, snp_df, by = "variant_id")
        # sort by position
        df <- df[order(df$pos), ]
        # any weirdness - remove
        df <- df[!is.na(df$pos) ,]
        chr <- stringr::str_split_fixed(.x, ".cis_qtl_pairs.", n =2)[,2]
        chr <- gsub(".parquet", "", chr)
        tempFile <- paste0(outFile, "_", chr)
        message( "* writing to: ", tempFile )
        # only write out header for chr1
        if( chr == "chr1" ){
            # write column names as separate file
            write_tsv(df[0,], path = header_file, col_names = TRUE)
            write_tsv(df, path = tempFile, col_names = FALSE) 
        }else{
            write_tsv(df, path = tempFile, col_names = FALSE)
        }
        rm(df)
        gc()        
})

# use cat to concatenate together 
message(" * concatenating to ", outFile)

all_temp <- paste0(outFile, "_chr", gsub("chr","",unique(snp_df$chr)), collapse = " ")
all_temp <- paste0(header_file, " ", all_temp)
cmd <- paste( "cat", all_temp, " > ", outFile )
message( cmd )  
system(cmd)
# remove temp files
cmd <- paste( "rm ", all_temp)
#system(cmd)
#write_tsv(all_res, path = outFile)

# sorting
message( "* sorting " )
sort_cmd <- paste0("cat ", outFile, " | (sed -u 1q; sort -k10,10 -k11,11n) > ", outFile, ".sorted")
message( " * ", sort_cmd)
system(sort_cmd)

# bgzip
message( "* bgzipping " )
bgzip_cmd <- paste0(" ml bcftools; bgzip -f -c ", outFile, ".sorted > ", outFile, ".gz")
message( " * ", bgzip_cmd)
system(bgzip_cmd)

# tabix
message(" * tabixing " )
tabix_cmd <- paste0(" ml bcftools; tabix -S 1 -s 10 -b 11 -e 11 ", outFile, ".gz" )
message(" * ", tabix_cmd)
system(tabix_cmd)

clean_cmd <- paste0("rm ", outFile, ".sorted" )

