#' ---
#' title: Comparing eQTL running with lm different inputs
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule comp_lm from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias/Snakefile


library(data.table)
library(ggplot2)

lm1 <- snakemake@input[['lm1']]
lm2 <- snakemake@input[['lm2']]


lm <- lapply(list(lm1, lm2), function(i) rbindlist(lapply(i,fread)))

lm <-Reduce(function(a,b) merge(a,b,by.x=c('Gene_id', 'tag'),
                                by.y=c('Gene_id','SNP')),
            lm)

p <- ggplot(lm, aes(log2.aFC_mean.x, log2.aFC_mean.y)) +
    geom_point(shape=1) +
    xlab("Inputs for lm") +
    ylab("Inputs for Deseq2")
print(p)
