source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")

#' Remove ASE from het inds in rasqual input to check if results are the same when genotypes are fix
#'
#' @param data name of vcf with rasqual input
#' @param out name to save output file
#' @keywords rasqual QC
#' @export
#' @return saves file with modified input to make vcf
#' rasq.mod()

rasq.mod <- function(data, out){
    dt <- vcf_w(data)
    samples <- gsub("_GT", "", grep("_GT", names(dt), value=T))

    
    ## replace counts for "0|0" or "1|1" with 0,0
    tmp <- lapply(samples, function(i) {
        g <- paste0(i,"_GT")
        a <- paste0(i, "_AS")
        dt[(get(g) == "0|0" | get(g) == "1|1") & get(a) != "0,0", eval(a):="0,0"]
        dt[, eval(g):=paste(get(g),get(a), sep=":")]
    })

    ## format as vcf: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO and FORMAT followed by samples

    dt[, QUAL:=100][,FILTER:="PASS"][,INFO:="RSQ=1"][,FORMAT:="GT:AS"]

    dt[, c("id", paste0(samples, "_AS")):=NULL]

    ##order and save

    dt <- setcolorder(dt, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", paste0(samples, "_GT")))
    
    write.table(dt, grep("mod", out, value=T), row.names=F,col.names=F,quote=F,sep="\t")



    ## save file control with RSQ=1 to make sure I am running rasqual n the same settings
    
    dt <- vcf_w(data)
    samples <- gsub("_GT", "", grep("_GT", names(dt), value=T))

    ## make GT:ASE cols
    tmp <- lapply(samples, function(i) {
        g <- paste0(i,"_GT")
        a <- paste0(i, "_AS")
        dt[, eval(g):=paste(get(g),get(a), sep=":")]
    })

    ## format as vcf: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO and FORMAT followed by samples

    dt[, QUAL:=100][,FILTER:="PASS"][,INFO:="RSQ=1"][,FORMAT:="GT:AS"]

    dt[, c("id", paste0(samples, "_AS")):=NULL]

    ##order and save

    dt <- setcolorder(dt, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", paste0(samples, "_GT")))
    
    write.table(dt, grep("control", out, value=T) ,row.names=F,col.names=F,quote=F,sep="\t")


}


rasq.mod(snakemake@input[['data']], snakemake@output[['data2']])
## data <- "/mrc-bsu/scratch/ev250/bin/rasqual/data/chr11.gz"

