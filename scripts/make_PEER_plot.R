# Plot number of QTLs found per number of PEER factors used 
# Jack Humphrey
library(tidyverse)
library(optparse)

# arguments
## out_folder - where the various results are kept
## mode - whether eQTL or sQTL
## data_code - the name of the dataset

option_list <- list(
    make_option(c('--prefix'), help='the name of the dataset being used', default = "example"),
    make_option(c('--dir'), help='path to directory containing subfolders with cis_qtl.txt.gz results', default = "results/example"),
    make_option(c('--FDR'), help ="the threshold for significance", default = 0.05)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

data_code <- opt$prefix
out_folder <- opt$dir
fdr <- opt$FDR

out_file <- paste0(out_folder, "/", data_code, "_N_per_PEER_plot.pdf")

# get all cis_qtl outputs 
sig_res <- list.files(out_folder, pattern = "cis_qtl.txt.gz", full.names = TRUE, recursive = TRUE)

# read in and filter at 5% FDR
sig_res_files <- map(sig_res, ~{read_tsv(.x ) })

# get number of tested features and number sig at FDR
sig_res_n_qtl <- map_dbl(sig_res_files, nrow)
sig_res_n_sig <- map_dbl(sig_res_files, ~{.x %>% filter(qval < fdr) %>% nrow()  })

sig_res_info <- gsub(".cis_qtl.txt.gz", "", basename(sig_res))


# get number of genes with at least one QTL
# for sQTL results grouped by cluster
sig_res_n_gene <- map_dbl(sig_res_files, ~{ 
    .x <- filter(.x, qval < fdr)
    if("group_id" %in% names(.x) ){ 
    .x$gene <- str_split_fixed(.x$group_id, ":", 2)[,1]
    length(unique(.x$gene) ) 
    }else{
      return(nrow(.x))
    } })

# make table
sig_res_df <- 
  stringr::str_split_fixed(sig_res_info, pattern = "_peer", n = 2) %>%
  as.data.frame( stringsAsFactors = FALSE ) %>%
  rename(dataset = V1, n_peer = V2) %>%
  mutate(n_peer = as.numeric(str_split_fixed(n_peer, "_", 2)[,1]) ) %>%
  mutate(n_features = sig_res_n_qtl) %>%
  mutate(n_sig_QTL = sig_res_n_sig) %>%
  mutate(n_genes = sig_res_n_gene) %>%  
mutate(file = basename(sig_res) )

print(as.data.frame(sig_res_df))

out_table <- paste0(out_folder, "/", data_code, "_N_per_PEER.tsv")

write_tsv(sig_res_df, file = out_table)


ymax <- max(sig_res_df$n_genes) + (max(sig_res_df$n_genes) * 0.1) 
xbreaks <- unique(sig_res_df$n_peer)

n_QTL_per_PEER_plot <-
  sig_res_df %>%
  ggplot(aes( y = n_genes, x = n_peer, group = dataset)) + 
  geom_point(aes(colour = dataset)) +
  geom_line(aes(colour = dataset), show.legend=FALSE) +
  geom_text(aes(label = n_genes), vjust = 0, nudge_y = 0.5, show.legend=FALSE, size = 3 ) +
  scale_x_continuous( breaks = xbreaks, labels = xbreaks ) + 
  xlab("Number of PEER factors") +
  ylab("Number of QTL genes") +
  ylim(0, ymax) +
  theme(legend.position = "bottom") +
  scale_colour_brewer(type = "qual") +
  labs(title = data_code, subtitle = paste0("FDR = ", fdr) ) +
  theme_bw()


ggsave(plot = n_QTL_per_PEER_plot, file = out_file)
