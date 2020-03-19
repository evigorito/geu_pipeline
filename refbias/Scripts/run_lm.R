source("/home/ev250/Bayesian_inf/trecase/Functions/lm.GT.R")
library(tictoc)

tic("Start running")

Gene=snakemake@wildcards[['gene']]
if(is.null(Gene)){
    Gene <- snakemake@params[['gene']]
}

prefix <- snakemake@params[['prefix']]
if(!exists("prefix")){
    prefix <- NULL
}

snps <- snakemake@params[['snps']]

if(is.null(snps)) {
    tmp <- fread(snakemake@input[['tags2run']])
    snps <- tmp[Gene_id == Gene, tag]
} else {
    if(!is.na(as.numeric(snps))){
    snps <- as.numeric(snps)
    }
}



tag <- snakemake@params[['tag']]

if(!is.na(as.numeric(tag))){
    tag <- as.numeric(tag)
}

nhets <- snakemake@params[['nhets']]

if(!is.null(nhets)){
    nhets <- as.numeric(nhets)
} else {
    nhets <- 5
}



pval <- snakemake@params[['pval']]
if(!exists("pval")){
    pval <- NULL
}

vcf <- snakemake@input[['vcf']]
if(!exists("vcf")){
    vcf <- NULL
}

gt <- snakemake@input[['gt']]
if(!exists("gt")){
    gt <- NULL
}




lm.gt(gene=Gene,
      chr=as.numeric(snakemake@params[['chrom']]),
      snps=snps,
      counts.f=snakemake@input[['counts']],
      covariates=snakemake@input[['libsize']],
      gene.coord=snakemake@input[['genecoord']],
      vcf=vcf,
      gt=gt,
      tag.threshold=tag,
      out=snakemake@params[['out']],
      prefix=prefix,
      pval=pval,
      nhets=nhets
      )


toc()
