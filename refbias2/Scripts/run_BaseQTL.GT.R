## Making arguments compatible across workflows calling same script

library(tictoc)
library(data.table)

prefix=snakemake@params[['prefix']]

if(!exists("prefix")){
    prefix <- NULL
}

nhets=snakemake@params[['nhets']]

if(!exists("nhets")){
    prefix <- NULL
}


min.ase.het=snakemake@params[['minasehet']]

if(!exists("min.ase.het")){
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

Gene=snakemake@wildcards[['gene']]
if(is.null(Gene)){
    Gene <- snakemake@params[['gene']]
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
print(head(snps))

tag <- snakemake@params[['tag']]

if(!is.na(as.numeric(tag))){
    tag <- as.numeric(tag)
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



source("/home/ev250/Bayesian_inf/trecase/Functions/BaseQTL.GT.refbiasCorrec.R")


tic("Start running")

    baseqtl.gt.refbias(gene=Gene,
                        chr=as.numeric(snakemake@params[['chrom']]),
                        snps=snps,
                        counts.f=snakemake@input[['counts']],
                        covariates=snakemake@input[['libsize']],
                        e.snps=snakemake@input[['eSNPs']],
                        u.esnps=ufsnps,
                        gene.coord=snakemake@input[['genecoord']],
                        vcf=snakemake@input[['vcf']],
                        le.file=snakemake@input[['leRef']],
                        h.file=snakemake@input[['hapRef']],
                        population=snakemake@params[['pop']],
                        tag.threshold=tag,
                        out=snakemake@params[['out']],
                        prefix=prefix,
                        model=snakemake@params[['model']],
                        stan.model=snakemake@input[['smodel']],
                       prior=prior,
                       nhets=nhets,
                       min.ase.het=min.ase.het,
                        stan.trec=snakemake@input[['trecmodel']],
                        ex.fsnp=snakemake@params[['pfsnp']],
                        AI_estimate=AI,
                        pretotalReads=as.numeric(snakemake@params[['pretotalReads']])
                        )

toc()
