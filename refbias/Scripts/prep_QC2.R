source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")


#' Prepares list of exonic SNPs that passed Fisher exact test testing the frequency of heterozigocity between the sample and referece panel to perform genotyping error QC2
#'
#' @param dir full path to dir with Fisher exact test output
#' @param pattern pattern to match files in dir
#' @param snps full file name with exonic SNPs pre-filtering
#' @param pval min pval to filter fSNPS, defaults to NULL
#' @param out file name to save filtered results
#' @export
#' @keywords inputs Fisher QC
#' @return saves file to use in btrecase.noGT functions
#' prep.QC()

prep.QC <- function(dir, pattern, snps, pval=NULL, out){
    
    fisher <- comb.files(path=dir, pattern=pattern)
    ## exclude fSNPs below cut-off
    if(!is.null(pval)) fisher <- fisher[pvalue>=pval,]
    
    fsnp <- fread(snps)
    fsnp[, id:=paste(POS,REF,ALT, sep=":")]
       
    temp <- merge(fsnp, fisher, by.x=c("id", "gene_id"), by.y=c("fsnp", "gene_id"))
    temp <- temp[, names(fsnp), with=F]
    temp[, id:=NULL]
    
    write.table(temp, file=out, row.names=F)

}

pval=snakemake@params[['pfsnp']]

if(!exists("pval")){
    pval <- NULL
} else {
    pval <- as.numeric(pval)
}

prep.QC(dir=snakemake@params[['FisherDir']],
        pattern=snakemake@params[['pattern']],
        snps=snakemake@input[['ueSNPS']],
        pval=pval,
        out=snakemake@output[['out']]
        )
    
