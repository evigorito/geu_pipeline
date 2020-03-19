
library(parallel)
suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

#' Extracts parameters from stan output across genes with signif associations
#'
#' @param dir full path to directory with stan outputs
#' @param par character vector with parameters to extract
#' @param out full path to write results
#' @export
#' @return saves a data table to file cols: Gene_id, model, par, rSNP
#' ex.par()

ex.par <- function(dir,par,out){
    ## get genes and tags run
    gt <- comb.files(path=dir, pattern = "ENSG[0-9]+\\.stan.summary.txt")

    ## select significant assoc only

    gt <- gt[log2_aFC_null == "no",]
    ## get posterior files
    post.f <- list.files(path=dir, pattern = "ENSG[0-9]+\\.stan.post.RDS", full.names=T)
    names(post.f) <- unlist(lapply(post.f, function(i) gsub("\\..*","", basename(i))))

    ## select post files with genes with signif assoc only

    post.f <- post.f[unique(gt$Gene_id)]

    ## sample param estimates
    post <- lapply(names(post.f), function(i) {
    tags  <- gt[ Gene_id == i, tag]
    tmp <- readRDS(post.f[i])
    l <- list()
    for(j in 1:length(tmp)){             
        if(j == 1){
            w <- which(names(tmp[[j]]) %in% tags)
            tmp2 <- tmp[[j]][w]
        } else {
            w <- which(names(tmp[[j]][[1]]) %in% tags)
            tmp2 <- tmp[[j]][[1]][w]
        }
        e <- rbindlist(mclapply(names(tmp2), function(k) {
            r <- rstan::extract(tmp2[[k]], pars=par )
            ## make data table including distance to origin
            edt <- data.table(bj=r$bj, rSNP=k)
            return(edt)
        }))
        l[[j]] <- e
    }
    names(l) <- names(tmp)
    return(rbindlist(l, idcol="model"))
    })
    names(post) <- names(post.f)
    write.table(post, out, row.names=F)
    
}
    

ex.par(dir=snakemake@params[['GT_dir']], snakemake@params[['par']], snakemake@output[['out']])
