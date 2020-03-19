## Making arguments compatible across workflows calling same script

library(tictoc)

prefix=snakemake@params[['prefix']]

if(!exists("prefix")){
    prefix <- NULL
}

ufsnps=snakemake@input[['ueSNPS']]
if(!exists("ufsnps")){
    ufsnps <- NULL
}

if(!is.null(snakemake@wildcards[['rbias']])){
    if(snakemake@wildcards[['rbias']] == "norefbias"){
        AI <- NULL
    } else {
        AI <- snakemake@input[['AI']]
    }
    
}else {
    AI <- snakemake@input[['AI']]
}

prior=snakemake@params[['prior']]
if(!exists("prior")) {
    prior <- NULL
} else {
    k=length(prior)/3 ## number of gaussians
    s <- seq(1,length(prior),k)
    l <- lapply(1:3, function(i) as.numeric(prior[s[i]: (s[i]+k-1)]))
    names(l) <- c("mean", "sd", "mix")
    prior <- l
    
}

exfsnp=snakemake@params[['pfsnp']]
if(!exists("exfsnp")){
    exfsnp <- NULL
} else {
    exfsnp <- as.numeric(snakemake@params[['pfsnp']])
}

negonly=snakemake@input[['model_negonly']]
if(!exists("negonly")){
    negonly <- NULL
}

sample=snakemake@input[['sample']]
if(!exists("sample")){
    sample <- NULL
}

source("/home/ev250/Bayesian_inf/trecase/Functions/Btrecase.noGT.refpanelBias.R")


tic("Start running")
btrecase.nogt.rna.refbias(gene=snakemake@wildcards[['gene']],
                  chr=as.numeric(snakemake@params[['chrom']]),
                  snps=as.numeric(snakemake@params[['snps']]),
                  counts.f=snakemake@input[['counts']],
                  covariates=snakemake@input[['libsize']],
                  e.snps=snakemake@input[['eSNPs']],
                  u.esnps=ufsnps,
                  gene.coord=snakemake@input[['genecoord']],
                  vcf=snakemake@input[['vcf']],
                  sample.file=sample,
                  le.file=snakemake@input[['leRef']],
                  h.file=snakemake@input[['hapRef']],
                  population=snakemake@params[['pop']],
                  maf=as.numeric(snakemake@params[['maf']]),
                  min.ase=as.numeric(snakemake@params[['minAse']]),
                  min.ase.snp=as.numeric(snakemake@params[['minAseSnp']]),
                  min.ase.n=as.numeric(snakemake@params[['minAseN']]),
                  tag.threshold=as.numeric(snakemake@params[['tag']]),
                  info=as.numeric(snakemake@params[['info']]),
                  out=snakemake@params[['out']],
                  prefix=prefix,
                  prob=snakemake@params[['prob']],
                  model=snakemake@input[['model']],
                  model.negonly=negonly,
                  prior=prior,
                  ex.fsnp=exfsnp,
                  AI_estimate=snakemake@input[['AI']],
                  pretotalReads=as.numeric(snakemake@params[['pretotalReads']])
                  )

toc()

## e.snps='/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/fSNP/chr11.fSNP.genes.txt'
## gene="ENSG00000233645"
## gene='ENSG00000002330'
## counts.f <-  "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/Psoriasis_skin.txt"
## covariates <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/Psoriasis_skin_gc_lib_size.rds"
## vcf <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/GT/chr11.ASE.Psoriasis_skin.vcf.gz"
## chr <- 11
## gene.coord <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt"

## sample="/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3.sample"

## le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr11.legend.gz'

## h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr11.hap.gz'

## model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan'

## population="EUR"
## maf=0.05
## nhets=5
## min.ase=5
## min.ase.snp=5
## min.ase.n=5
## tag.threshold=.9
## q.test="no"
## info=0.3
## snps=5*10^5
## prob=0.99
## ex.fsnp=10^-4
## out='/mrc-bsu/scratch/ev250/psoriasis/Btrecase/output99/ENSG00000002330'
## prefix = "ENSG00000002330.Psoriasis_skin"

## test=btrecase.nogt.rna (gene, chr, snps=5*10^5,counts.f,covariates,e.snps,gene.coord,vcf,le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), maf=0.05, min.ase=5,min.ase.snp=5,min.ase.n=5,tag.threshold=.9,q.test="no", info=0.3, out, prefix=prefix, model="/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan", prob, ex.fsnp)
