#script for using existing genotype PCs, sample key, and covariates file from eQTL/sQTL to run edQTLs
#Winston Cuddleston

library(tidyverse)
library(optparse)

option_list <- list(make_option(c('--inputDir'), help = 'the path to the QTL inputs directory', default = './'),
                    make_option(c('--jacusaDir'), help = 'the path to the Jacusa project directory with ratio file', default = ''),
                    make_option(c('--rat'), help = 'name of the editing ratio matrix', default = 'all_sites_pileup_editing.tsv.gz'),
                    make_option(c('--previous_key'), help = '', default = ''),
                    make_option(c('--edQTL_key'), help = 'name of the sample key for edQTL mapping', default = ''),
                    make_option(c('--previous_PCs'), help = '', default = ''),
                    make_option(c('--edQTL_PCs'), help = 'name of the genotype PCs file for edQTL mapping', default = ''),
                    make_option(c('--metadata'), help = 'full path and name of the metadata file containing ADAR expression', default = ''),
                    make_option(c('--previous_covars'), help = '', default = ''),
                    make_option(c('--edQTL_covars'), help = 'name of the covariates file for edQTL mapping', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

inputDir <- opt$inputDir
jacusaDir <- opt$jacusaDir
rat <- opt$rat
previous_key <- opt$previous_key
edQTL_key <- opt$edQTL_key
previous_PCs <- opt$previous_PCs
edQTL_PCs <- opt$edQTL_PCs
previous_covars <- opt$previous_covars
edQTL_covars <- opt$edQTL_covars
metadata <- opt$metadata


#covars file for QTL mapping has to match genotype ID not RNA-seq ID

rat <- read_tsv(paste0(jacusaDir,rat))
previous_key <- read_tsv(paste0(inputDir,previous_key))
previous_PCs <- read_tsv(paste0(inputDir,previous_PCs))
previous_covars <- read_tsv(paste0(inputDir,previous_covars))
metadata <- readRDS(metadata)

#make sure that genotype PCs and sample key only contain samples that were not dropped during editing pipeline 
rat <- rat[colnames(rat) %in% previous_key$sample_id]
key <- previous_key[previous_key$sample_id %in% colnames(rat),]
a <- key$participant_id
b <- append(a, "ID")
covars <- previous_covars[,which(names(previous_covars) %in% b)]

PCs <- previous_PCs %>% 
  column_to_rownames(., var = "ID") %>%
  select(any_of(key$participant_id)) %>%
  rownames_to_column(., var = "ID")

ADARs <- metadata[,c("sample_id", "participant_id.x", "ADAR1", "ADAR2")]

x <- merge(ADARs, key, by = "sample_id")

covars_t <- x %>%
  pivot_longer(names_to = "ID", values_to = "value", cols = !participant_id) %>%
  pivot_wider(names_from = "participant_id", values_from = "value")

fin <- rbind(covars_t,covars)
FIN <- fin[-c(1,2),]
FIN[,-1] <- sapply(FIN[,-1],as.double)

write_tsv(covars, file = paste0(inputDir,edQTL_covars))
write_tsv(PCs, file = paste0(inputDir,edQTL_PCs))
write_tsv(key, file = paste0(inputDir,edQTL_key))
