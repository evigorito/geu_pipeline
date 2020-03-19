source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")

##' Saves a simplified version of full vcf file, removing homs and/or missing genotypes for all samples

dna <- vcf_w(snakemake@input[['vcf_GT']])

saveRDS(dna, snakemake@output[['out']])
