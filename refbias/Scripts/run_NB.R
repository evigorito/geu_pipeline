prefix=snakemake@params[['prefix']]

if(!exists("prefix")){
    prefix <- NULL
}

source("/home/ev250/Bayesian_inf/trecase/Functions/NegBinom.R")    

print("start")
neg.binom(gene=snakemake@wildcards[['gene']],
          chr=as.numeric(snakemake@params[['chrom']]),
          snps=as.numeric(snakemake@params[['snps']]),
          counts.f=snakemake@input[['counts']],
          covariates=snakemake@input[['libsize']],
          gene.coord=snakemake@input[['genecoord']],
          vcf=snakemake@input[['vcf']],
          tag.threshold=as.numeric(snakemake@params[['tag']]),
          out=snakemake@params[['out']],
          prefix=prefix
          )



