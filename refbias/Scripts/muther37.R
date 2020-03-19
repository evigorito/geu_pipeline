library(rtracklayer)
library(liftOver)
library(data.table)

#' Format MuTHER data file: uplift built36 to optional built and add reference/alt allele to each SNP. Input format provides col A1 as "ref/alt allele" but only provides 1 allele.
#'
#' based on https://bioconductor.org/packages/devel/workflows/vignettes/liftOver/inst/doc/liftov.html
#' @param chain full name to chain file to use in UCSC format
#' @param muther full name to file with muther results, one file per chromosome
#' @param chrom chromosome of muther input file
#' @param out full name to save formatted file
#' @export
#' @return saves file

lift <- function(chain, muther,chrom, out){
    ## if chain ends with gz, unzip
    if(grep(".gz$", chain)){
        system(paste("gzip -d ", chain))
        chain <- sub(".gz$", "", chain)
    
    }
    
    ch <- import.chain(chain)

    ## make Granges object from muther
    ## start with dt
    dt <- fread(cmd=paste("zcat",  muther))
    ## add chr col in UCSC format
    dt[, chr:=paste0("chr", chrom)]
    
    gr <- makeGRangesFromDataFrame(dt,
                                   keep.extra.columns=T,
                                   start.field="SNP_Coor",
                                   end.field="SNP_Coor")
                                   
    ##liftOver
    seqlevelsStyle(gr) = "UCSC"
    gr19 <- liftOver(gr, ch)
    gr19 <- unlist(gr19)
    genome(gr19) <- "hg19"

    ## convert to back to dt

    dt <- data.table(data.frame(gr19))
    dt[, end:=NULL]
    setnames(dt, "start", "POS")

    ## save 
    write.table(dt, row.names=F, file=out)
    
}

chain <- snakemake@input[['chain']]
muther <- snakemake@input[['muther']]
out <- snakemake@output[['mutherF']]
chrom <- snakemake@params[['chrom']]

lift(chain, muther,chrom, out)
