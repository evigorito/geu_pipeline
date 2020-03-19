#' ---
#' title: Comparing fSNP GT QC measures
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule Btrecase_analysis from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias/Snakefile



library(cowplot)
source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")

dna <- vcf_w(snakemake@input[['dna']])
rna <- name(snakemake@input[['rna']])
fsnps=snakemake@input[['eSNPS']]
fisher <- snakemake@input[['fisher_het']]
HW <- snakemake@input[['HW']]

rna_stan=snakemake@params[['rna_stan']]
dna_dir=snakemake@params[['dna_dir']]
dna_fsnps=snakemake@params[['dna_fsnps']]
depth=as.numeric(snakemake@params[['depth']])


## dna <- vcf_w("/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/DNA/chr22.ASE.allsamples.vcf.gz")
## rna <- name("/mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.txt")
## fsnps <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/inputs/fSNP/chr22.fSNP.unique.genes.txt"
## fisher <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/QCfSNP/fisher.fSNPs.txt"
## HW <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/QCfSNP/HWE.fSNPs.txt"

## rna_stan <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/RNA'
## dna_dir <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/GT'
## dna_fsnps <- "rbias.ENSG[0-9]+.*.\\.GT.fsnps.with.counts.rds"
## depth=10



rna[, id:=paste(POS, REF, ALT, sep=":")]

## get fSNPs run in inference in DNA
fsnp.i <- lapply(list.files(dna_dir, pattern=dna_fsnps, full.names=T), readRDS)
fsnp.i <- unique(unlist(fsnp.i))

##remove .n .m ending
fsnp.i <- unique(unlist(lapply(fsnp.i, function(i) gsub(".[n:m]$", "", i))))

## link fSNPs to genes
fsnp <- fread(fsnps)
fsnp[,id:=paste(POS, REF,ALT, sep=":")]

## select genes run in RNA
rna.genes <- list.files(path=rna_stan, pattern="^rbias\\.ENSG[0-9]+.*stan.summary.txt")
rna.genes <- gsub(".*(ENSG[0-9]+).*", "\\1", rna.genes)

## select fSNPs used in GT inference
fsnp.i <- fsnp[id %in% fsnp.i,][gene_id %in% rna.genes,]


## merge dna/rna with fsnp.i to get genes and fsnps
rna <- merge(rna, fsnp.i[,.(gene_id, id)], by="id")
dna <- merge(dna, fsnp.i[,.(gene_id, id)], by="id")

## add Fisher pvalue, HWE pvalue and Allele feq pvalue to rna

fish <- fread(fisher)
HW <- fread(HW)
print(head(names(rna)))
print(names(fish))
rna <- merge(rna, fish[, .(gene_id, fsnp, pvalue)] , by.x=c("gene_id", "id"), by.y=c("gene_id", "fsnp"))
rna <- merge(rna, HW, by=c("gene_id", "id"))

## Get GT errors by pvalue and min depth
rna.dna.dp <- convert2(dna, rna, chr=22)

err.fish <-  err2.var(col="pvalue", x=depth, DT=rna.dna.dp)
err.HW <- err2.var(col="HW.p", x=depth, DT=rna.dna.dp)
err.A0 <- err2.var(col="RP.p", x=depth, DT=rna.dna.dp)

l <- mapply( function(i,j) {
    ggplot(data=i) +
    geom_point(mapping=aes(y=1-prop.conc, x=-log(pval, 10)), shape=1) +
    labs(title= j,
        y="Proportion of errors",
        x=expression(paste(-log[10], "(p-value)")))
},
i=list(err.fish , err.HW, err.A0),
j=c("Heterozygous frequency to RP", "Hardy-Weinberg equilibrum", "Allele frequency relative RP" ),
SIMPLIFY=F)


#+ fig.width= 6, fig.height=7.5
plot_grid(plotlist=l, ncol=1, labels="auto")

## a corresponds to the test I have done so far.
## b corresponds to Hardy-Weinberg equilibrum and
## c compares the proprtion of reference allele in the sample with the reference panel.
## I selected fSNPs called with depth >=10.
