library(data.table)
library(rasqualTools)
library(biomaRt)

#' get library size corrected by gc content as matrix with rownames genes and colnames samples saved in rds format
#' 
#' @param counts.f list or character vector with names of raw count files
#' @param rds list or character vector with names of rds files to save
#' @export
#' @return saves matrix as rds file
#' counts2gc_cov()

counts2gc_cov <- function(counts.f, rds){
    tmp <- mapply(function(a,b) {
        counts <- fread(a)
        ## convert to matrix
        counts_matrix <- as.matrix(counts[,2:ncol(counts)], rownames=counts$gene_id)
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
       
        saveRDS(size_factors, b)
    },
    a=counts.f,
    b=rds,
    SIMPLIFY=F)
}

counts2gc_cov(counts.f=snakemake@input[['counts']], rds=snakemake@output[['gc_cov']])
    
