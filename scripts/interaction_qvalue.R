library(tidyverse)
library(qvalue)
library(optparse)


option_list <- list(
    make_option(c('--input'), help='', default = "example"),
        make_option(c('--output'), help='', default = "example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


input <- opt$input
output <- opt$output
df <- read_tsv(input)

do_qvalue <- function(x){ qvalue(x)$qvalues}

df$qval_g <- do_qvalue(df$pval_g)
df$qval_i <- do_qvalue(df$pval_i)
df$qval_gi <- do_qvalue(df$pval_gi)

write_tsv(df, output)
