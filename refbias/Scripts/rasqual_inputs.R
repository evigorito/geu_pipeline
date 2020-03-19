

library(rasqualTools)
library(data.table)
library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)

source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")


##' # Prepare inputs for rasqual
##'
##' @param counts full file name to data table with raw counts, firt col gene_id, followed by samples
##' @param cov name for file with covariate matrix: row names samples and cols covariates
##' @param ebg name for file with a genomic ranges list, each element exons for a gene
##' @param snps vcf file to extract snps from
##' @param cis numeric with length of cis-window, defaults to 5*10^5
##' @param out path to dir to save rasqual inputs
##' @param gc.correct, whether to correct lib size by gc content, defaults to YES, otherwise it is centred and scaled and suffix for files is "libsize_centred_scaled"
##' @return saves files with rasqual inputs
##' @export
##'
##' rasq.in()

rasq.in <- function(counts,cov,ebg,snps, cis=5e5, out, gc.correct=c("rasqual_gc","libsize_centred_scaled") ){

    ## load objects
    counts <- fread(counts)
   
    ##covs <- readRDS(cov)
    ebg <- readRDS(ebg)
    snps <- vcf_w(snps)
    gc.correct <- match.arg(gc.correct)

    ## transform counts in matrix as required by rasqualTools

    counts_matrix <- as.matrix(counts[,2:ncol(counts)], rownames=counts$gene_id)
    saveRasqualMatrices(list(counts=counts_matrix), out , file_suffix = "rasqual")

    ## covariates

    if(gc.correct == "rasqual_gc"){ ## use rasqual function
    ## need gene GC content
        ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
        
        gc <- data.table(getBM(attributes = c("ensembl_gene_id",  "transcript_length", "percentage_gene_gc_content"),
                               filters = "ensembl_gene_id",
                               values = counts$gene_id,
                               mart = ensembl))
        
        setkey(gc, ensembl_gene_id, transcript_length)

        ## get gc and length for the longest transcript
        gc <- gc[, .SD[.N], ensembl_gene_id]
        setnames(gc, names(gc), c("gene_id", "length", "percentage_gc_content"))

        size_factors =rasqualCalculateSampleOffsets(counts_matrix, gc ,gc_correct=TRUE)
       

    } else { ## use log(lib.size) centred and scaled, rasqual normalises total counts diving by the mean

        lib_size = scale(log(colSums(counts_matrix)), center=TRUE,scale=TRUE)
        size_factors = matrix(rep(lib_size, nrow(counts_matrix)), nrow = nrow(counts_matrix), byrow = TRUE)
        rownames(size_factors) = rownames(counts_matrix)
        colnames(size_factors) = colnames(counts_matrix)        
    }
     saveRasqualMatrices(list(covariates=size_factors),  out, file_suffix = gc.correct)
    

    ## Calculate the number of SNPs overlapping each gene

    ## Need a data table with gene_id, chr, strand exon_starts and exon_ends cols as in rasqual instructions

    ebg <- ebg[counts$gene_id]
    gene_id <- names(ebg)
    chrom <- unlist(lapply( seqnames(ebg), function(i) as.character(i@values)))
    st <- as.data.table(start(ebg))[, paste(value, collapse=","), group_name]
    end <- as.data.table(end(ebg))[, paste(value, collapse=","), group_name]
    strand <- unlist(lapply( strand(ebg), function(i) as.character(i@values)))

    dt <- data.table(gene_id=gene_id, chr=chrom, strand=strand, exon_starts=st$V1, exon_ends=end$V1)

    ## Select cols in snps and named thyem as required by rasqualtools

    snps <- snps[,.(CHROM, POS)]
    snps[, snp_id:=paste(CHROM, POS, sep=":")]
    setnames(snps, c("CHROM", "POS"), c("chr", "pos"))
    
    snp_counts = countSnpsOverlapingExons(dt, snps, cis_window = cis)
    snp_select <- dplyr::select(snp_counts, gene_id, feature_snp_count, cis_snp_count)

    ## Add dt to snp_select

    snp_select <- as.data.table(merge(snp_select, dt, by="gene_id"))

    ## Select genes with rsnps
    snp_select <- snp_select[ cis_snp_count>0, ]

    write.table(snp_select,file=paste0(out, '/snps_perGene.txt'), row.names=F, quote=F)

    }

gc.correct <- snakemake@params[['gc_correct']]
if(!exists("gc.correct")){
    gc.correct <- "rasqual_gc"
}


rasq.in(counts=snakemake@input[['counts']],
        cov=snakemake@input[['libsize']],
        ebg=snakemake@input[['ebg']],
        snps=snakemake@input[['snps']],
        cis=as.numeric(snakemake@params[['cis']]),
        out=snakemake@params[['out']],
        gc.correct=gc.correct)
