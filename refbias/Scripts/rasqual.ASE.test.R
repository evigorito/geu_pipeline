#' ---
#' title: Comparing rasqual output with/without "ASE" info for homozygous calls
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule Btrecase_analysis from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias/Snakefile

## Load R packages and functions:
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
##source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")

## open rasqual outputs

rasq.files <- snakemake@input[['rasq']]
 
rasq.header <- snakemake@params[['rasqheader']]


## rasq.files <- list.files(path='/mrc-bsu/scratch/ev250/bin/rasqual/data/output', full.names=T)
## rasq.header <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/output", "rasqual.header.txt", full.names=T)

inp <- unique(gsub(".*\\.([a-z]*)\\..*", "\\1", basename(rasq.files)))



rasq <- lapply(inp,
               function(i) {
                   f <- grep(i , rasq.files, value=T)
                   tmp <- rbindlist(lapply(f, function(j) format_rasqual(j, header=rasq.header, top.hits="no")))
                   tmp[, ASE:=i]
                   })
                    

rasq.w <- Reduce(function(a,b) merge(a,b, by=names(rasq[[1]])[1:7], suffixes=paste0(".", inp)),
                 rasq)

## plot 

effects <- ggplot(rasq.w, aes(Fold_change.control, Fold_change.mod)) + geom_point() + geom_abline() + ggtitle("Effect estimates control  vs modified input")

pval <- ggplot(rasq.w, aes(-log(p.control, 10), -log(p.mod, 10))) + geom_point() + geom_abline() + ggtitle("P-val controls vs modified input")

plot_grid(effects, pval, ncol=1)

