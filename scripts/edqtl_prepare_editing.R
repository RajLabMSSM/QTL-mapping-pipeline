library(optparse)
library(dplyr)
library(tibble)
library(tidyverse)
library(preprocessCore)

option_list <- list(make_option(c('--vcf_chr_list'), help = "list of chromosomes to use", default = ''),
                    make_option(c('--rat'), help = "name of raw editing ratio output from Jacusa", default = ''),
                    make_option(c('--anno'), help = "name of annotation file for editing sites from Jacusa", default = ''),
                    make_option(c('--keyIn'), help = "full path and name of sample key already filtered post-jacusa pipeline", default = ''),
                    make_option(c('--pheno'), help = "phenotype matrix file in BED format", default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

chrs <- opt$vcf_chr_list
rat <- opt$rat
anno <- opt$anno
sampleKey <- opt$keyIn
phenotype_matrix <- opt$pheno

chrs <- read.delim(chrs, header = F, sep = "\t")
rat <- read_tsv(rat)
anno <- read_tsv(anno)
key <- read_tsv(sampleKey)

rat <- rat %>%
  column_to_rownames(., var = "ESid") %>%
  select(any_of(key$sample_id)) %>%
  rownames_to_column(., var = "ESid")

keyfix <- data.frame(t(key))
colnames(keyfix) <- keyfix['sample_id',]
keyfix <- keyfix[!(rownames(keyfix) == 'sample_id'),]
ratios <- rat %>% column_to_rownames(., var = "ESid")
ratios <- rbind(keyfix,ratios)
colnames(ratios) <- ratios['participant_id',]
Rat <- ratios[!(rownames(ratios) == 'participant_id'),]
RAT <- lapply(Rat, function(x) as.numeric(as.character(x)))
RAT <- data.frame(RAT, row.names = rownames(Rat), check.names = F)

rat <- RAT %>% rownames_to_column(., var = "ESid")
rat[is.na(rat)] <- 0
AGrat <- rat[grep("A:G", rat$ESid),]
TCrat <- rat[grep("T:C", rat$ESid),]
CTrat <- rat[grep("C:T", rat$ESid),]
GArat <- rat[grep("G:A", rat$ESid),]
pheno <- rbind(AGrat,TCrat,CTrat,GArat)

norm <- pheno %>%
  remove_rownames(.) %>%
  column_to_rownames(., var = "ESid") %>%
  scale(., center = T, scale = T) %>%
  preprocessCore::normalize.quantiles(., keep.names = T) %>%
  data.frame(., check.names = F) %>%
  rownames_to_column(., var = "ESid")

geneIDs <- subset(anno, select = c("ESid2", "ensembl_id", "strand"))
normAnno <- merge(geneIDs, norm, by.x = "ESid2", by.y = "ESid") %>% rename(., "ESid" = "ESid2")

bed <- normAnno %>%
  mutate(gene_id = paste(ESid, ensembl_id, sep = ":")) %>%
  mutate(group_id = gene_id) %>%
  separate(ESid, c('Chr','end', 'ref', 'alt')) %>%
  mutate(., chr = gsub("chr", "", Chr)) %>%
  arrange(., chr, end) %>%
  mutate(start = as.numeric(end) - 1) %>%
  relocate(start, .before = 'end') %>%
  relocate(gene_id, .after = 'end') %>%
  relocate(group_id, .after = 'gene_id') %>%
  relocate(strand, .after = 'group_id') %>%
  select(-c(ref, alt, ensembl_id, chr)) %>%
  rename(., '#Chr'='Chr')
  
colnames(bed) <- gsub("\\.","-",colnames(bed))
final_bed <- bed[bed$`#Chr` %in% chrs$V1,]

write_tsv(final_bed, phenotype_matrix)
