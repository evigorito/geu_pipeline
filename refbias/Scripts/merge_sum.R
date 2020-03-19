## Load R packages and functions:
library(data.table)

source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R", verbose=FALSE)
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R', verbose=FALSE)


#' Prepare summary of associations  and compute r2 between rsnps per gene
#'
#' @param stan.dir path to directory with input stan summaries and tags
#' @param prefix character witht the prefix used when creating stan summaries, defaults to NULL
#' @param tags whether SNPs were grouped and tagged or not, defaults to tags
#' @param eSNPs full name to file with fSNPs
#' @param stan_sum full name to file to save stan summary
#' @param legend full name for legend file to extract EAF
#' @return saves stan summary combining all genes
#' @export
#' merge_stan()

merge_stan <- function( stan.dir, prefix=NULL, tags="yes", eSNPs,  stan_sum, legend){
        
    stan_all <- comb.files(path=stan.dir, pattern=paste0(prefix,".noGT.stan.summary.txt"))

    fisher <- comb.files(path=stan.dir, pattern=paste0(prefix, ".fsnps.het.fisher.test.txt"))

    ## get fSNPs from fisher output and add EAF
    
    fsnps <- fisher[,.(gene_id,fsnp)]

    if(tags == "yes"){

        tags <- comb.files(path=stan.dir, pattern=paste0(prefix, ".noGT.eqtl.tags.lookup.txt"))
   
        tags.run <- merge(tags, unique(stan_all[,.(Gene_id,tag)]), by= c("Gene_id","tag"))

        eaf <- rbindlist(lapply(unique(fsnps$fsnp), function(i) snp.eaf(legend,i)))

        fsnps <- merge(fsnps, eaf, by.x="fsnp", by.y= "snp")
        
        setnames(fsnps, c("fsnp","eaf"), c("fSNP.id","fSNP.EAF"))

        fsnps.run <- merge(fsnps,tags.run,  by.y=c("Gene_id", "SNP"), by.x=c("gene_id","fSNP.id"), all.x=T)
        
        fsnps.run <- fsnps.run[!is.na(tag),]

        ## add fSNPs to stan_all

        stan_all <- merge(stan_all, fsnps.run, by.x=c("Gene_id", "tag"), by.y=c("gene_id", "tag"), all.x=T)
        
        stan_all[!is.na(fSNP.id) , tag.fsnp.op:="no"][(fSNP.EAF < 0.5 & tag.EAF > 0.5) | (fSNP.EAF > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]
        
    } else {

        fsnps[ , fSNP:= "yes"]
        stan_all <- merge(stan_all, fsnps, by.x=c("Gene_id", "SNP"), by.y=c("gene_id", "fsnp"), all.x=T)
    }
        
    
    ## save stan_all
    write.table(stan_all, stan_sum, row.names=F)  

}



eSNPs <- snakemake@input[['eSNPs']]
legend <- snakemake@input[['legfile']]
stan.dir <- snakemake@params[['out_stan']]
prefix <- snakemake@params[['prefix']]
tags <- snakemake@params[['tags']]
stan_sum <- snakemake@output[['comb']]

if(!exists("prefix")){
    prefix <- NULL
}

if(!exists("tags")){
    tags <- "yes"
}


merge_stan(stan.dir, prefix, tags, eSNPs, stan_sum, legend)

