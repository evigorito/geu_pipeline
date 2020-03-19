#' ---
#' title: Comparing eQTL with or without reference bias correction for selected gene-SNP associations
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule Btrecase_analysis from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias/Snakefile

## Load R packages and functions:
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(cowplot)


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/aux_btrecase.R")


##' # Extract and format data for analysis

## Get stan summaries

## Look at genes run with GT with ref bias correction

gt <- comb.files(path=snakemake@params[['GT_dir']], pattern = "ENSG[0-9]+\\.stan.summary.txt")

## gt <- comb.files(path="/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2", pattern = "ENSG[0-9]+\\.stan.summary.txt")


##  Look at genes run with GT without ref bias correction
gt_old <- comb.files(path=snakemake@params[['GT_old']], pattern="ENSG[0-9]+\\.stan.summary.txt")
## gt_old <- comb.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT/BtrecaseV2', pattern="ENSG[0-9]+\\.stan.summary.txt")

## REmove associations with extreme values for Rhat (poorly run MCMC models)

## both models struggled with some SNPS in ENSG00000211666
summary(gt$log2_aFC_Rhat)
summary(gt_old$log2_aFC_Rhat)

gt <- gt[log2_aFC_Rhat < 1.1 & log2_aFC_Rhat > 0.9,]
gt_old <- gt_old[log2_aFC_Rhat < 1.1 & log2_aFC_Rhat > 0.9,]

gt.comb <- merge(gt, gt_old, by=c("Gene_id", "tag"), suffixes=c(".post", ".pre"))


## select associations run with ASE with correction (lower number of fSNPs)

gt.ase <- gt.comb[model.post == "trec-ase",]

gt.ase <- add.signif(gt.ase, x1="log2_aFC_null.pre", x2="log2_aFC_null.post", col=c("Without","With") )


## get stan inputs for genes with ASE info

## f.in <- list.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2',pattern=".GT.stan1.input.rds", full.names=T)


f.in <- list.files(path=snakemake@params[['GT_dir']],pattern=".GT.stan1.input.rds", full.names=T)
stan.in <- lapply(f.in, readRDS)

names(stan.in) <- sapply(f.in, function(i) gsub("\\..*", "", gsub(".*GT2/", "",i)))

## convert format to ease analysis

inp <- lapply(stan.in, function(i) lapply(i, in.neg.beta.prob.eff2))

## start with GT, look for significant eQTL in close LD with a fSNP. Either the fSNP is eQTL or tagged at 0.9

## Get fSNPs and their target gene

vcf.dna <- snakemake@input[['rds_GT']]
## vcf.dna='/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/inputs/chr22.ASE.allsamples.remHom.rds'

dna <- readRDS(vcf.dna)


## Get fsnps and their host gene 
## fsnps <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt")

## unique fSNPs
##ufsnps <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/inputs/fSNP/chr22.fSNP.unique.genes.txt")
##ufsnps[, id:=paste(POS,REF,ALT, sep=":")]

fsnps <- fread(snakemake@input[['eSNPs']])

## select fsnps from p.dna
fsnps[, id:=paste(POS,REF,ALT, sep=":")]


## for each fSNP in dna I need the target gene

gene.fsnp <- merge(dna[, c("CHROM","POS", "REF","ALT","id"), with=F], fsnps, by=c("CHROM","POS", "REF","ALT","id"))
    

## for each fsnp run, get the tagging SNP and EAF

gt.tags <- comb.files(path=snakemake@params[['GT_dir']],pattern="eqtl.tags.lookup.txt")

## gt.tags <- comb.files(path = '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2', pattern="eqtl.tags.lookup.txt")

gt.tags.run <- merge(gt.tags, gt.ase[,.(Gene_id,tag)], by= c("Gene_id","tag"))

tag.fsnp <-  merge(gene.fsnp, gt.tags.run, by.x=c("gene_id", "id"), by.y =c("Gene_id", "SNP") )
## add EAF for SNP col to check direction of effects with tag
##DT <- snp.eaf(file1='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz', snps = tag.fsnp$id)

DT <- snp.eaf(file1=snakemake@input[['leRef']], snps = tag.fsnp$id)
tag.fsnp <- merge(tag.fsnp, DT, by.x="id", by.y="snp")
setnames(tag.fsnp, c("id","eaf"), c("fSNP.id","EAF.fsnp"))
    
## use long format  and add 2 cols: fSNP and op.dir. fSNP is "yes" if tag is tagging an fSNP and op.dir tells whether the EAF of tag and fSNP are in the same direction

long <- function(x,y, xn, yn){
    ## make long format from 2 data tables with extra column identifying source of x and y, coded in xn (x name) and yn (y name). Use the same associations in both data tables to ease comparison
    x[, source:=xn]
    y[, source:=yn]
    x <- merge(x, y[, .(Gene_id, tag)], by=c("Gene_id","tag"))
    y <- merge(y, x[, .(Gene_id, tag)], by=c("Gene_id","tag"))
    rbind(x,y, fill =T)
}

## select only trec-ase model, ASE info
gt.l <- long(gt[model=="trec-ase",], gt_old[model=="trec-ase",], "GT-post", "GT-pre")

## add fsnps

gt.l <- merge(gt.l, tag.fsnp[,c("gene_id", "tag", "fSNP.id","EAF.fsnp") ,with=F], by.x=c("Gene_id", "tag"), by.y=c("gene_id", "tag"), all.x=T)
gt.l[, fSNP:="yes"][is.na(EAF.fsnp), fSNP:="no"]
gt.l[fSNP=="yes" , tag.fsnp.op:="no"]
gt.l[(EAF.fsnp < 0.5 & tag.EAF > 0.5) | (EAF.fsnp > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]


## add trec

## trec <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.stan.summary.txt")

trec <- comb.files(path=snakemake@params[['trec']], pattern="trec.stan.summary.txt")

## get trec with same tags as gt.l run with ASE info
trec.sub <- merge(trec, unique(gt.l[,.(Gene_id, tag, tag.EAF)]), by=c("Gene_id", "tag"))
trec.sub[, source:="NegBinomial-GT"]

## trec.sub has same tags as gt.l, compare same associations
trec.sub <- merge(trec.sub, tag.fsnp[,c("gene_id", "tag", "fSNP.id","EAF.fsnp") ,with=F], by.x=c("Gene_id", "tag"), by.y=c("gene_id", "tag"), all=T)

trec.sub[, fSNP:="yes"][is.na(EAF.fsnp), fSNP:="no"]
trec.sub[fSNP=="yes" , tag.fsnp.op:="no"][(EAF.fsnp < 0.5 & tag.EAF > 0.5) & (EAF.fsnp > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]


## combine with gt.l
gt.trec <- rbind(gt.l, trec.sub, fill=T)
gt.trec[, Signif:="No significant"][log2_aFC_null=="no", Signif:="Significant"]


## Some fSNPs are in different direction as tag SNP, create new column to adjust eQTL estimates

gt.trec[, log2mean.eaf:=log2_aFC_mean][tag.fsnp.op=="yes", log2mean.eaf:= -log2_aFC_mean]
gt.trec[,source:=factor(source, levels=c( "GT-pre" ,   "GT-post" , "NegBinomial-GT"))]



##################################################################################################################################

##' # Look at gene-snp association ENSG00000211685 23261726:C:A which I previously suspected reference panel bias (out.chr22.lm.GT.noGT.rna.v3.R line 1413)

##' This is a complex locus even if few fSNPs. The fSNPs are close to each other and some are not observed. Some of the reads must be overlapping more than 1 fSNP so helps to see how the estimates are affected.

## fsnps
fsnps.g <- gene.fsnp[ gene_id == 'ENSG00000211685',id]

## snps to look for 
g.snps <- c('23261726:C:A', fsnps.g)

## get tags
tags.g <- gt.tags.run[SNP %in% g.snps & Gene_id == 'ENSG00000211685',]

## fSNP "23264842:G:A" wasnt run as tag due to "rsnp with less than 5 het ind."

ref.test <- stan.in[['ENSG00000211685']][['23261726:C:A']]

dt.ref <- inp.qc.gt(ref.test)

## Comapre eQTL estimates using different models and ref bias correction
gt.trec[Gene_id == 'ENSG00000211685' & tag == '23261726:C:A',]

## very strong effect with ASE only without ref panel correction
## grep 23261726:C:A /mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/noGTrsnp/ASE.only/ENSG00000211685*

## Look at estimates for reference panel bias estimates
## ai.est <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/pre_post_AI.txt')

ai.est <- fread(snakemake@input[['AI']])

## get ai for fSNPs
ai.est[POS >= 23264785 & POS<=23264864,]


gt.as <- dna[id%in% g.snps,]

## recode GT in gt.as

rec.rs <- rec_mytrecase_rSNPs(x=gt.as$POS, y=gt.as)

##' # Look in more detail at the estimates, the estimates are not as estreme as the data and they were estimated with large number of counts

## Look at the variability of the ref panel bias estimates across samples
 
## dir.ai <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/"
dir.ai <- snakemake@params[['dir_ai']]

samp <- unique(gsub("_GT", "", grep("_GT", names(gt.as), value=T)))

l.ai <- lapply(samp, function(i) {
    tmp <- fread(paste0(dir.ai,"/",  i, "/Aligned.sortedByCoord.out.post_remapping_AI.txt"))
    tmp[ , id:=paste(POS,REF,ALT, sep=":")]
    tmp <- tmp[id %in% g.snps,]
    tmp <- merge(tmp, rec.rs[ , c("id", paste0(i, "_GT")), with=F], by="id")
    names(tmp)[ncol(tmp)] <- "GT"
    return(tmp)
})

names(l.ai) <- samp

ai.samp <- rbindlist(l.ai, idcol="sample")
ai.samp[, rb.est:=NALT/(NALT+NREF)]
ai.samp[,rb.sd:= sqrt(1/NALT + 1/NREF)]

      
## fSNP '23264842:G:A' doesn't have ref bias estimate due to low maf

## get observed ASE per sample per fSNP, for homozygous SNPs NaN

as <- grep("_AS",names(gt.as), value=T)
tmp <- gt.as[,as,with=F]

 ## get counts for hap2 (n)
hap2 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], ","), `[[`, 2))))
 ## get counts for hap1+2 (m)
hap12 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], ","), function(j) sum(as.numeric(j))))))
ase.obs <- hap2/hap12
colnames(ase.obs) <- samp
rownames(ase.obs) <- g.snps


dt <- data.table(id=rep(rownames(ase.obs), ncol(ase.obs)), sample=rep(samp,each=nrow(ase.obs)), ai.obs=c(ase.obs))

## combine ase obs with ref panel bias estimates per sample

comb <- merge(ai.samp, dt, by=c("sample", "id"))

## need to correct ase.obs, when g==-1, 1-ase.obs
comb[ ,ai.corr.obs:=ai.obs][GT==-1, ai.corr.obs:=1-ai.obs]


## collect mean effect and sd for each fSNP

mean.eff <- ai.samp[, mean(rb.est,na.rm=T), id]
sd <- merge(ai.samp[, 1/sum(NREF), id], ai.samp[, 1/sum(NALT), id], by="id")
sd[, V1:=V1.x+V1.y]
sd[,sd:=sqrt(V1)]

##upper limit
sd[,V1:=2*sd + mean.eff$V1]
mean.eff <- rbind(mean.eff, sd[,.(id,V1)])
##lower
sd[,V1:=-2*sd + mean.eff$V1[1:3]]
mean.eff <- rbind(mean.eff, sd[,.(id,V1)])

## When calculating ASE I can only use het individuals. Check whether estimates for ref panel bias are affected by the genotype of the sample.
## I have assumed that reads from het individuals are a random sample from total reads


p1 <- ggplot(comb, aes(x=log(NREF+NALT), y=rb.est, color=as.factor(abs(GT)))) + geom_point(shape=1) +
    geom_hline(yintercept=0.5, colour="blue")+   
    facet_grid(.~id) +
    ylab("Reference panel bias \n estimate per sample")+
    xlab("log(Total re-aligned counts to fSNP)") +
    geom_hline(data = mean.eff, aes(yintercept=V1), linetype = "dashed", size=0.5)+
    ggtitle("Estimates for reference panel bias per sample")

p2 <- ggplot(comb, aes(x=log(NREF+NALT), y=ai.corr.obs, color=as.factor(abs(GT)))) + geom_point(shape=1) +
    geom_hline(yintercept=0.5, colour="blue")+
    ylab("Observed allelic imbalance \n at each fSNP per sample")+
    xlab("log(Total re-aligned counts to fSNP)") +
    facet_grid(.~id) +
    geom_hline(data = mean.eff, aes(yintercept=V1), linetype = "dashed", size=0.5)+
    ggtitle("Observed allelic imbalance per sample")


#' Upper plot corresponds to the reads used for calculating estimates for reference panel bias. The horizontal lines correspond to the ref panel bias estimate calculated by polling all reads  plus/mins 2*SD (as used when modelling eQTL effects). The lower plot corrresponds to the observed ASE in het samples.

#''23261726:C:A' is not listed as fSNP though it has ASE counts associated, probably due to variations in annotation files.

plot_grid(p1,p2 ,ncol=1)


## look at the fSNPs, proxy for fSNP 23264785:C:A

ref.test <- stan.in[['ENSG00000211685']][['23258438:A:G']]
dt.ref <- inp.qc.gt(ref.test)

## look at fSNP, proxi for fSNP  23264864:T:A
ref.test <- stan.in[['ENSG00000211685']][['23268677:G:A']]
dt.ref <- inp.qc.gt(ref.test)


########################################################################################################################

## Look for association of fSNPs where eQTL pre <0 and Signif and with trec or ref bias correction is not

gt.ase <- merge(gt.ase, gt.l[fSNP=="yes", .(Gene_id,tag, fSNP, tag.fsnp.op)], by=c("Gene_id", "tag"))

gt.ase.f <- gt.ase[fSNP=="yes" & log2_aFC_mean.pre <0 & log2_aFC_null.pre == "no" & log2_aFC_null.post == "yes",]


## Most clear example : ENSG00000100030 22190163:C:A
## Compare estimates
unique(gt.ase.f[,.(Gene_id, tag, log2_aFC_mean.pre, log2_aFC_mean.post, log2_aFC_null.pre , log2_aFC_null.post)])
trec[Gene_id  %in% gt.ase.f$Gene_id & tag %in% gt.ase.f$tag,.(Gene_id, tag, log2_aFC_mean, log2_aFC_null)]


ref.test <- stan.in[['ENSG00000100030']][['22190163:C:A']]
dt.ref <- inp.qc.gt(ref.test)


## Opposite direction

gt.ase.f <- gt.ase[fSNP=="yes" & log2_aFC_mean.post >0 & log2_aFC_null.post == "no" & log2_aFC_null.pre == "yes",]

## Compare estimates
unique(gt.ase.f[,.(Gene_id, tag, log2_aFC_mean.pre, log2_aFC_mean.post, log2_aFC_null.pre, log2_aFC_null.post)])
trec[Gene_id  %in% gt.ase.f$Gene_id & tag %in% gt.ase.f$tag,.(Gene_id, tag, log2_aFC_mean, log2_aFC_null)]


g <- 'ENSG00000189060'
s <- '38138284:A:T'

ref.test <- stan.in[[g]][[s]]
dt.ref <- inp.qc.gt(ref.test)

g <-  'ENSG00000183473'
s <- '37580334:G:A'

ref.test <- stan.in[[g]][[s]]
dt.ref <- inp.qc.gt(ref.test)


##' # Conclusions

##' The reference panel bias correction seems to be working as expected, I can extend the analysis but it appears that the reads mapped to hets are similar to the reads mapped to hom.
##' For  ENSG00000211685- 23261726:C:A the direction of the eQTL estimate from Trec is different from Btrecase but applying ref panel bias correction get closer to trec. it doesnt seem to be only due to reference panel bias, probably some discrepancy between total counts and mapped counts, though it is only one exon gene.

##' I selected a couple of examples in which the ref panel bias correction makes the eQTL estimates more in line with trec, either by being or not significant. For simplicity in interpretation I have chosen fSNPs.


