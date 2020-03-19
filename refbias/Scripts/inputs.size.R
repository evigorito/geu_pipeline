library(data.table)

#' prepares inputs for running subsamples of 5 and 25 GBR individuals
#'
#' @param counts.f name of file with filtered raw counts to extract individuals
#' @param libsize.f name of rds file with library size covariate for samples
#' @param size number of samples to subset. Take nested random samples of indicated sizes.
#' @param out path to dir to write output files
#' @export
#' @return saves files to out dir
#' inp.sub()

inp.sub <- function(counts.f, lisize.f, size, out){
    counts <- fread(counts.f)
    lisize <- readRDS(lisize.f)
    ## sort size decreasing=T to start from bigger to smaller
    size <- sort(size, decreasing=T)
    ## nested sampling
    samp <- list()
    i=1
    x=rownames(lisize)
    while(i <= length(size)){
        samp[[i]] <- sort(sample(x, size[i]))
        x <- samp[[i]]
        i <- i+1
    }

    ## extract samp from counts and lisize so they are in the same order
    counts.sub <- lapply(samp, function(i) counts[, c("gene_id", i), with=FALSE])
    lisize.sub <- lapply(samp, function(i) lisize[i, , drop=F])

    ## save files
    for(s in 1:length(size)){
        i=size[s]
        write.table(counts.sub[[s]], paste0(out, "/b37_filtered.raw_counts.", i ,"inds.txt"), row.names=F)
        saveRDS(lisize.sub[[s]], paste0(out, "/library.size.",i ,"inds.rds"))
        writeLines(samp[[s]], paste0(out,"/samples.",i ,"inds.txt"))
    }
    
}

counts.f <- snakemake@input[['counts']]
lisize.f <- snakemake@input[['libsize']]
size <- as.numeric(snakemake@params[['sample']])
out <- snakemake@params[['out']]

## counts.f <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt"
## lisize.f <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.comp.rds"
## size=c(5,25)
## out <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/inputs/sample5_25"

inp.sub(counts.f, lisize.f, size, out)
