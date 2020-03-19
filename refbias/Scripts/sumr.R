library(tictoc)
tic("Start running")
## Load R packages and functions:
library(data.table)
library(ggplot2)
library(parallel)
mc.cores = snakemake@threads

source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R", verbose=FALSE)
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R', verbose=FALSE)
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

#' Prepare summary of associations combining Normal and psoriatic skin and compute r between rsnps per gene
#'
#' @param path to dir with skin summaries
#' @param gene_coord text file with cols gene_id and chrom, other cols allowed
#' @param lefile character vector with full names to reference panel legend files per chromosome
#' @param hapfile character vector with full names to reference panel hap files per chromosome
#' @param sample full name for sample file linked to legend and haps to select samples, defaults to NULL
#' @param group character with group to select (col in samples, ID, POP, GROUP,SEX), defaults to NULL
#' @param maf value to filter snps, defaults to 0.05 which was used to select rsnps to run
#' @param r full name to file to save list of matrices with r between rsnps per gene
#' @return save list of matrices with r between rsnps per gene
#' @export
#' sumr()

sumr <- function(stan.dir, gene_coord,  lefile, hapfile,sample=NULL, group=NULL, maf=0.05, r){

    sumall <- rbindlist(lapply(stan.dir, function(i) comb.files(path=i, pattern=paste0("\\.ENSG[0-9]+.*stan.summary.txt"))), fill=T)
    sumall <- unique(sumall[, .(Gene_id, tag)])

    geneStEnd <- fread(gene_coord)
    sumall <- merge(sumall, geneStEnd[, .(gene_id, chr)], by.x="Gene_id", by.y="gene_id")

    ## get r, add tag.pos to sum1t
    sumall[, tag.pos:=as.numeric(gsub(":.*", "", tag))]
    setkey(sumall, Gene_id, tag.pos)

    tagr <- mclapply(unique(sumall$Gene_id), function(i){
        dt <- sumall[Gene_id==i,]
        ##file1 <- grep(paste0("chr",unique(dt[,chrom]),".legend.gz"), lefiles,value=T)
        ##file2 <- grep(paste0("chr",unique(dt[,chrom]),".hap.gz"), hapfiles,value=T)
        cw <- c(min(dt$tag.pos), max(dt$tag.pos))
        mat <- haps.range(file1=lefile,file=hapfile,cw=cw, maf=maf)
        mat <- mat[unique(dt$tag),]
        
        if(!is.null(sample) & !is.null(group)){
            samples <- fread(sample)
            ## get row numbers to keep
            keep <- samples[GROUP == group,which=T]
            ## for each sample there are 2 cols in haps file, so I need keep*2 and keep*2-1
            keep.cols <- sort(c(keep*2, keep*2-1))
            mat <- mat[,keep.cols, drop=F]
        }
        
        r <- cor2(t(mat)) 
        return(r)
    },
    mc.cores=mc.cores)
    names(tagr) <- unique(sumall$Gene_id)

    ## save tagr
    saveRDS(tagr, r)    

}

stan.dir <- snakemake@params[['btrecase']]
gene_coord <- snakemake@input[['gene_coord']]
lefile <- unlist(snakemake@input[['leRef']])
hapfile <- unlist(snakemake@input[['hapRef']])
sample <- snakemake@input[['sample']]
r <- snakemake@output[['r']]
group <- snakemake@params[['group']]

sumr(stan.dir, gene_coord, lefile, hapfile, sample, group, maf=0, r)

## gene_coord="/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"
## stan.dir <- "/mrc-bsu/scratch/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna"
## lefiles <- list.files(path="/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3",
##                       pattern="1000GP_Phase3_chr[0-9]*.\\.legend.gz",
##                       full.names=T)

## hapfiles <- list.files(path="/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3",
##                       pattern="1000GP_Phase3_chr[0-9]*.\\.hap.gz",
##                       full.names=T)

## sample <- "/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3.sample"

## r='/mrc-bsu/scratch/ev250/psoriasis/Btrecase/objects/fisher001EURr.rds'
## group="EUR"
toc()
