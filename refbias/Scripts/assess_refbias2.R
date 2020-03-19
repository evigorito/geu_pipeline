#' ---
#' title: Comparing eQTL with or without reference bias correction
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


##' # Compare eQTL estimates with GT with and without reference bias correction

## Get stan summaries

## Look at genes run with GT with ref bias correction

gt <- comb.files(path=snakemake@params[['GT_dir']], pattern = "ENSG[0-9]+\\.stan.summary.txt")

## gt <- comb.files(path="/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2", pattern = "ENSG[0-9]+\\.stan.summary.txt")


##  Look at genes run with GT without ref bias correction
gt_old <- comb.files(path=snakemake@params[['GT_old']], pattern="ENSG[0-9]+\\.stan.summary.txt")
##gt_old <- comb.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT/BtrecaseV2', pattern="ENSG[0-9]+\\.stan.summary.txt")

## REmove associations with extreme values for Rhat (poorly run MCMC models)

## both models struggled with some SNPS in ENSG00000211666
summary(gt$log2_aFC_Rhat)
summary(gt_old$log2_aFC_Rhat)

gt <- gt[log2_aFC_Rhat < 1.1 & log2_aFC_Rhat > 0.9,]
gt_old <- gt_old[log2_aFC_Rhat < 1.1 & log2_aFC_Rhat > 0.9,]

gt.comb <- merge(gt, gt_old, by=c("Gene_id", "tag"), suffixes=c(".post", ".pre"))

gt.all <- merge(gt, gt_old, by=c("Gene_id", "tag"), suffixes=c(".post", ".pre"), all.y=T)

gt.all <- add.signif(gt.all, x1="log2_aFC_null.pre", x2="log2_aFC_null.post", col=c("Without","With") )

cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%", "log2_aFC_null")

## look at changes from trec to trec-ase, following gene gain association when only trec was run 

btrecase.plot(dt=gt.all[Signif != "None",] , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without reference panel bias correction"
              )


## select associations run with ASE with correction (lower number of fSNPs)

gt.ase <- gt.comb[model.post == "trec-ase",]

gt.ase <- add.signif(gt.ase, x1="log2_aFC_null.pre", x2="log2_aFC_null.post", col=c("Without","With") )


## Count significant associations with or without reference panel bias correction
table2 <- gt.ase[,.N,Signif]
names(table2) <- c("Signif", "SNPs")

table2[Signif=="None", color:="#999999"][Signif=="Without", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="With", color:="#0072B2"]

table2[, Signif:=factor(Signif, levels=c("None","Without", "With","Both"))]

setkey(table2,Signif)

p.ase <- btrecase.plot(dt=gt.ase[Signif != "None",] , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without \n reference panel bias correction"
              )+
    annotation_custom(tableGrob(table2[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table2$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)

ggsave(file=paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/plots_tables/GTRefBiasChr22ASE.png"),
     p.ase, width=4,
           height=4)


## btrecase.plot(dt=gt.ase , x1=paste0(cols, ".pre"),
##               x2=paste0(cols, ".post"),
##               xl="eQTL effect without correction",
##               yl="eQTL effect with corrrection",
##               col=c("Without", "With"),
##               title="eQTL estimates with or without reference panel bias correction"
##               )+
##     annotation_custom(tableGrob(table2[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table2$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)


##' # Compare allelic imbalance estimates for haplotypes per individual with and without random effect

## get stan inputs for genes with ASE info

## f.in <- list.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2',pattern=".GT.stan1.input.rds", full.names=T)


f.in <- list.files(path=snakemake@params[['GT_dir']],pattern=".GT.stan1.input.rds", full.names=T)
stan.in <- lapply(f.in, readRDS)

names(stan.in) <- sapply(f.in, function(i) gsub("\\..*", "", gsub(".*GT2/", "",i)))

## convert format to ease analysis

inp <- lapply(stan.in, function(i) lapply(i, in.neg.beta.prob.eff2))

## Collect AI estimate and sd for each ind/gene-snp stratifying by significant associations
a0.sig <- list()
sd0.sig <- list()

a0.nosig <- list()
sd0.nosig <- list()

for(g in names(inp)){ ## gene level
    sig <- gt[Gene_id==g & log2_aFC_null=="no" & model=="trec-ase", tag]
    a0.sig[[g]] <- Reduce("c", lapply(inp[[g]][names(inp[[g]])%in%sig], '[[', "ai0"))
    a0.nosig[[g]] <- Reduce("c", lapply(inp[[g]][!names(inp[[g]])%in%sig], '[[', "ai0"))
    sd0.sig[[g]] <- Reduce("c", lapply(inp[[g]][names(inp[[g]])%in%sig], '[[', "sdai0"))
    sd0.nosig[[g]] <- Reduce("c", lapply(inp[[g]][!names(inp[[g]])%in%sig], '[[', "sdai0"))
}

ai.all <- lapply(list(a0.sig, a0.nosig, sd0.sig, sd0.nosig), function(i) Reduce("c", i))
names(ai.all) <- c("a0.sig", "a0.nosig", "sd0.sig", "sd0.nosig")

## sample from a0 using hap sd, by sig

ai.r <- mapply(function(a,b) rnorm(length(a), a, b), SIMPLIFY=FALSE, a=ai.all[grep("a0", names(ai.all), value=T)], b=ai.all[grep("sd0", names(ai.all), value=T)])

dt <- data.table(AI_fix=c(ai.all$a0.sig, ai.all$a0.nosig),AI_random=c( ai.r$a0.sig, ai.r$a0.nosig), SD=c(ai.all$sd0.sig, ai.all$sd0.nosig), Significant=rep(c("Signif", "Not signif"), sapply(ai.r, length)), Gene_id = unlist(mapply(function(a,b) rep(a, length(b)), SIMPLIFY=FALSE, a=c(names(a0.sig), names(a0.nosig)), b=c(a0.sig, a0.nosig))))


## plot equal number of haplotye reference panel bias estimates from significant and no significant gene-snp associations

## need to sample from dt as otherwise too many points

n=100000
cap <- paste("Sampling", n, "associations per signif. level")

ggplot(dt[,.SD[sample(.N,n)], by=Significant] , aes(x=AI_fix, y=SD)) +
    geom_point(alpha=1/5, size=1, shape=1) +
    scale_y_continuous(name="hap SD", breaks=seq(0,1,0.25), limits=c(0, 0.5)) +
    scale_x_continuous(name="Allelic imbalance per haplotype per ind (log)") +
    facet_grid(.~Significant) +
    labs(title="Pooled haplotypic allelic imbalance estimates",
         subtitle="Stratified by eQTL significance",
         caption=cap)


p <- ggplot(dt[,.SD[sample(.N,n)], by=Significant] , aes(x=AI_fix, y=AI_random)) +
    geom_point(alpha=1/5, size=1, shape=1) +
    facet_grid(.~Significant) +
    labs(title="Pooled haplotypic allelic imbalance estimates",
         subtitle="Stratified by eQTL significance",
         caption=cap)

## plot eQTL effect in log scale (currently log2)

gt[,log_aFC:=log2_aFC_mean/log2(exp(1))]


p2 <- ggplot(gt[model=="trec-ase" ,], aes(x=log_aFC, color=log2_aFC_null)) + geom_density() +
    scale_x_continuous("log eQTL-effect", seq(-2.5,2,0.5)) + ggtitle("eQTL effect") + labs(color = "aFC_null")

plot_grid(p, p2, ncol=1)


## Look at the eQTL estimates with and whithout correction for extreme hap allelic imbalance

genes.ex <- unique(dt[Significant=="Signif" & AI_fix >abs(.4), Gene_id])

## Look at all haps for genes.ex to see how variable the estimates are within gene

ggplot(dt[Gene_id %in% genes.ex ,] , aes(x=AI_fix, y=AI_random, color=Gene_id)) +
    geom_point(alpha=1/5, size=1, shape=1) +
    facet_grid(.~Significant) +
    labs(title="Pooled haplotypic allelic imbalance estimates",
         subtitle="Stratified by eQTL significance",
         caption=cap)

## Look at the change in eQTL effect at those genes

btrecase.plot(dt=gt.ase[Gene_id %in% genes.ex & Signif != "None" ,] , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without reference panel bias correction"
              )


##############################################################################################

##' # Look at effect size of cis-SNP stratified by being fSNP

## start with GT, look for significant eQTL in close LD with a fSNP. Either the fSNP is eQTL or tagged at 0.9

## Get fSNPs and their target gene

vcf.dna <- snakemake@input[['rds_GT']]
## vcf.dna='/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/inputs/chr22.ASE.allsamples.remHom.rds'

dna <- readRDS(vcf.dna) 

## Get fsnps and their host gene 
##fsnps <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt")

## unique fSNPs
##ufsnps <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/inputs/fSNP/chr22.fSNP.unique.genes.txt")

fsnps <- fread(snakemake@input[['eSNPs']])

## select fsnps from p.dna
fsnps[, id:=paste(POS,REF,ALT, sep=":")]
ufsnps[, id:=paste(POS,REF,ALT, sep=":")]

## for each fSNP in dna I need the target gene

gene.fsnp <- merge(dna[, c("CHROM","POS", "REF","ALT","id"), with=F], fsnps, by=c("CHROM","POS", "REF","ALT","id"))
    

## for each fsnp run, get the tagging SNP and EAF

gt.tags <- comb.files(path=snakemake@params[['GT_dir']],pattern="eqtl.tags.lookup.txt")

## gt.tags <- comb.files(path = '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2', pattern="eqtl.tags.lookup.txt")

gt.tags.run <- merge(gt.tags, gt.ase[,.(Gene_id,tag)], by= c("Gene_id","tag"))

tag.fsnp <-  merge(gene.fsnp, gt.tags.run, by.x=c("gene_id", "id"), by.y =c("Gene_id", "SNP") )
## add EAF for SNP col to check direction of effects with tag
DT <- snp.eaf(file1='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz', snps = tag.fsnp$id)
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

## prepare tables to add to plot

## get all combinations
g <- data.table(expand.grid(unique(gt.trec$Signif), levels(gt.trec$source)))



tables3 <- lapply(1:nrow(g), function(i) {
    dt <- gt.trec[Signif==g[i,Var1] & source == g[i,Var2] ,.N,.(fSNP,sign(log2mean.eaf))]
    dt[, sign:=as.character(sign)][sign=="-1", sign:="<0"][sign=="1", sign:=">0"]
    setkey(dt, fSNP,sign)
    setnames(dt, c("sign", "N"), c("Effect", "SNPs"))
    return(dt)
})

pal <- brewer.pal(n = 8, name = "Set1")[c(4,3)]

gl <- lapply(tables3, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 7,padding = unit(c(1.5, 1,1,1), "mm"), core=list(fg_params=list(col=c( rep(pal, each=2), rep("black", 8)))))))   

gt.trec <- gt.trec[order(match(source, levels(g$Var2))),]

## make sure names in dt are the same as in gt.trec (Signif and source)
dt <- data.table(Signif=g$Var1, source=g$Var2, grob = gl )

ref.p <- ggplot(gt.trec, aes(log2mean.eaf, color= fSNP)) +
    xlab("eQTL-effect") +
    ylab("Density") +
    geom_density()+
    scale_color_manual(values=pal) +
    theme_bw() +
    geom_vline(xintercept= 0, linetype = "dashed",  color = "grey")+
    facet_grid(rows= vars(source), cols=vars(Signif), scales = "free_x")+
    theme(strip.background=element_rect(fill="white"))+
    geom_custom(data=dt, aes(grob=grob), x = 0.75, y = 0.75)

ggsave(file="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/plots_tables/GTRefBiasChr22fSNPs.png",
       ref.p, width=15, height=13, uni="cm")

###############################################################################################

##' # Compare eQTL effects with GT (Trec and full model)  and external databases

## Start by full model and trec

gt.trec.w <- merge(gt[model=="trec-ase",],trec, by=c("Gene_id", "tag"), suffixes=c(".full",".trec"))

gt.trec.w <- add.signif(gt.trec.w, x1="log2_aFC_null.full", x2="log2_aFC_null.trec", col=c("full","trec") )

## Count significant associations with or without reference panel bias correction
table3 <- gt.trec.w[,.N,Signif]
names(table3) <- c("Signif", "SNPs")

table3[Signif=="None", color:="#999999"][Signif=="full", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="trec", color:="#0072B2"]

table3[, Signif:=factor(Signif, levels=c("None","full", "trec","Both"))]

setkey(table3,Signif)

btrecase.plot(dt=gt.trec.w[Signif != "None",] , x1=paste0(cols, ".full"),
              x2=paste0(cols, ".trec"),
              xl="eQTL effect full model",
              yl="eQTL effect trec model",
              col=c("full", "trec"),
              title="eQTL estimates with full or trec model"
              )+
    annotation_custom(tableGrob(table3[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table3$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)


## Apply multiple testing correction to Trec

trec.mult <- rej.recode(a=0.6, b=trec, c=0.9)

ci.mult(0.6, b=trec.mult)


## Look at Gtex ebv 

##ebv <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz")
ebv <- fread(cmd=paste("zcat", snakemake@input[['gtex']]))

##ebv.sig <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz")

ebv.sig <- fread(cmd=paste("zcat", snakemake@input[['sigGtex']]))

## Array Express and Chris file ##### AE: All the associations below false discovery rate 5%

## sig.qtl <- fread('/home/ev250/rds/rds-cew54-wallace-share/Data/GEUVADIS/Geuvadis_sig_eqtl') ## Chris file

sig.qtl <- fread(snakemake@input[['geu']])

## select  chrom 22 and format compatible with gt

sig.qtl <- sig.qtl[CHROM==22,][,Gene_id:=gsub("\\..*","", Gene)]
sig.qtl[, id := paste(POS,REF,ALT, sep=":")][, "null":="no"]
##geu <- geu[CHR_GENE==22,][,Gene_id:=gsub("\\.*","", GENE_ID)]

ebv[, Gene_id:=gsub("\\..*","",gene_id)]
ebv <- ebv[Gene_id %in% gt$Gene_id,]
ebv[, SNP:=gsub("^22_|_b37$", "", variant_id)][, SNP:=gsub("_", ":",SNP)]

ebv.sig <- ebv.sig[variant_id %in% ebv$variant_id,][, null:="no"]
ebv <- merge(ebv,ebv.sig[,.(gene_id, variant_id, null)], by=c("gene_id", "variant_id"), all.x=T)
ebv[is.na(null), null:="yes"]


## merge with trec.mult, previous version was with trec

trec.ebv <- merge(trec.mult, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), all.x=T)
trec.ebv <- add.signif(trec.ebv, x1="null.rej", x2="null", col=c("trec", "Gtex-ebv"))

btrecase.plot(dt=trec.ebv[Signif != "None"] , x1=c(cols[1:3], "null.rej"),
              x2=c(rep("slope",3), "null"),
              xl="eQTL effect Trec (log2)",
              yl="eQTL effect Gtex-EBV (slope)",
              col=c("Trec-mult", "Gtex-EBV"),
              title="eQTL estimates by Trec with multiple testing correction \n or Gtex-EBV"
              ) 

trec.sig <- merge(trec.mult,sig.qtl,by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
trec.sig[is.na(null), null:="yes"]
trec.sig <- add.signif(trec.sig, x1="null.rej", x2="null", col=c("Trec", "GEUVADIS"))


table3 <- trec.sig[,.N,Signif]
names(table3) <- c("Signif", "SNPs")

table3[Signif=="None", color:="#999999"][Signif=="Trec", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="GEUVADIS", color:="#0072B2"]

table3[, Signif:=factor(Signif, levels=c("None","Trec", "GEUVADIS","Both"))]

setkey(table3,Signif)


btrecase.plot(dt=trec.sig[Signif != "None"] , x1=paste0(cols),
              x2=c(rep("Beta",3), "null"),
              xl="eQTL effect Trec (log2)",
              yl="eQTL effect GEUVADIS (Beta)",
              col=c("Trec", "GEUVADIS"),
              title="eQTL estimates by Trec 85 samples or whole GEU"
              ) +
    annotation_custom(tableGrob(table3[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table3$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2) +
    labs(subtitle="Only significant GEUVADIS estimates available")


## merge with gt

gt.ebv <- merge(gt,ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), all.x=T)
gt.ebv <- add.signif(gt.ebv, x1="log2_aFC_null", x2="null", col=c("full", "Gtex-ebv"))

btrecase.plot(dt=gt.ebv[Signif != "None"] , x1=paste0(cols),
              x2=c(rep("slope",3), "null"),
              xl="eQTL effect BI-WI (log2)",
              yl="eQTL effect Gtex-EBV (slope)",
              col=c("BI-WI", "Gtex-EBV"),
              title="eQTL estimates by full model or Gtex-EBV"
              ) 


gt.sig <- merge(gt,sig.qtl,by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
gt.sig[is.na(null), null:="yes"]
gt.sig <- add.signif(gt.sig, x1="log2_aFC_null", x2="null", col=c("BI-WI", "GEUVADIS"))



table3 <- gt.sig[,.N,Signif]
names(table3) <- c("Signif", "SNPs")

table3[Signif=="None", color:="#999999"][Signif=="BI-WI", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="GEUVADIS", color:="#0072B2"]

table3[, Signif:=factor(Signif, levels=c("None","BI-WI", "GEUVADIS","Both"))]

setkey(table3,Signif)


btrecase.plot(dt=gt.sig[Signif != "None"] , x1=paste0(cols),
              x2=c(rep("Beta",3), "null"),
              xl="eQTL effect BI-WI (log2)",
              yl="eQTL effect GEUVADIS (Beta)",
              col=c("Trec", "GEUVADIS"),
              title="eQTL estimates by full model  85 samples or whole GEU"
              ) +
    annotation_custom(tableGrob(table3[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table3$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2) +
    labs(subtitle="Only significant GEUVADIS estimates available")





################################################################################################
##' # Rreference panel bias correction without GT

## eQTl estimates without corrections
old <- fread(snakemake@input[['stan_before']])

##old <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/RNAnoGTrsnp/rna99.summary.txt')

## old is rna99 without reference bias correction run on the same SNPs

## open corrected summaries

new <- fread(snakemake@input[['ref_bias']])

##new <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/results/rna99.BiasCorrec2.summary.txt")

## merge old and new

all <- merge(new, old, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), suffixes=c(".post", ".pre") )

all <- add.signif(all, x1="log2_aFC_null.pre", x2="log2_aFC_null.post", col=c("Without","With") )

cols <- c("log2_aFC_mean", "log2_aFC_0.5%", "log2_aFC_99.5%", "log2_aFC_null")

## Count significant associations with or without reference panel bias correction
table1 <- all[,.N,Signif]
names(table1) <- c("Signif", "SNPs")

table1[Signif=="None", color:="#999999"][Signif=="Without", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="With", color:="#0072B2"]

table1[, Signif:=factor(Signif, levels=c("None","Without", "With","Both"))]

setkey(table1,Signif)

btrecase.plot(dt=all[Signif != "None"] , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without reference panel bias correction"
              )+
    annotation_custom(tableGrob(table1[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table1$color, rep("black", 4)))))), xmin=-1.0, xmax=-.5, ymin=0.5, ymax=1)
                          



##################################################################################################################################

## Look at gene-snp association ENSG00000211685 23261726:C:A which I previously suspected reference panel bias (out.chr22.lm.GT.noGT.rna.v3.R line 1413)


ref.test <- stan.in[['ENSG00000211685']][['23261726:C:A']]

dt.ref <- inp.qc.gt(ref.test)

gt.no.corec <- readRDS('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT/ENSG00000211685.GT.stan1.input.rds')

noref.test <- gt.no.corec[['23999991:A:T']]

inp.qc.gt(noref.test)

gt.ase[Gene_id == 'ENSG00000211685' & tag == '23261726:C:A',]

gt.trec[Gene_id == 'ENSG00000211685' & tag == '23261726:C:A',]

## very strong effect with ASE only without ref panel correction
## grep 23261726:C:A /mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/noGTrsnp/ASE.only/ENSG00000211685*

ai.est <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/pre_post_AI.txt')

## get ai for fSNPs
ai.est[POS>= 23264785 & POS<=23264864,]



g.snps <- c('23261726:C:A', '23264785:C:A','23264842:G:A', '23264864:T:A')

gt.as <- dna[id%in% g.snps,]

## recode GT in gt.ase

rec.rs <- rec_mytrecase_rSNPs(x=gt.as$POS, y=gt.ase)


## Compare reference panel bias estimates per samples

dir.ai <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/"

samp <- unique(gsub("_GT", "", grep("_GT", names(gt.as), value=T)))

l.ai <- lapply(samp, function(i) {
    tmp <- fread(paste0(dir.ai, i, "/Aligned.sortedByCoord.out.post_remapping_AI.txt"))
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


plot_grid(p1,p2 ,ncol=1)


## look at fSNP, proxi for fSNP 23264785:C:A

gt.tags.run[SNP %in% g.snps & Gene_id == 'ENSG00000211685',]

ref.test <- stan.in[['ENSG00000211685']][['23258438:A:G']]
dt.ref <- inp.qc.gt(ref.test)

## look at fSNP, proxi for fSNP  23264864:T:A
ref.test <- stan.in[['ENSG00000211685']][['23268677:G:A']]
dt.ref <- inp.qc.gt(ref.test)


##' # Conclusions
##'
##' 





######### Ignore for the time being #########

##### We need to adjust for multiple testing ######

## post.f <- list.files(path=snakemake@params[['GT_dir']], pattern = "ENSG[0-9]+\\.stan.post.RDS", full.names=T)

## ## post.f <- list.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT2', pattern = "ENSG[0-9]+\\.stan.post.RDS", full.names=T)

## ## select signif associations
## ## calculate % of posterior out of region of rejection
## ## caterpillar plot for the extremes
## ## after ranking compare with p-value GEU and/or gtex

## gt.Sig <- gt[log2_aFC_null == "no",]

## post.fsub <- unlist(lapply(unique(gt.Sig$Gene_id), function(i) grep(i,post.f, value=T)))
## names(post.fsub)  <- unique(gt.Sig$Gene_id)

## post <- lapply(names(post.fsub)[1:10], function(i) {
##     tags  <- gt.Sig[ Gene_id == i, tag]
##     tmp <- readRDS(post.fsub[i])
##     l <- list()
##     for(j in 1:length(tmp)){             
##         if(j == 1){
##             w <- which(names(tmp[[j]]) %in% tags)
##             tmp2 <- tmp[[j]][w]
##         } else {
##             w <- which(names(tmp[[j]][[1]]) %in% tags)
##             tmp2 <- tmp[[j]][[1]][w]
##         }
##         e <- rbindlist(mclapply(names(tmp2), function(k) {
##             r <- rstan::extract(tmp2[[k]], pars="bj" )
##             ## make data table including distance to origin
##             edt <- data.table(bj=r$bj, rSNP=k)
##             return(edt)
##         }))
##         l[[j]] <- e
##     }
##     names(l) <- names(tmp)
##     return(rbindlist(l, idcol="model"))
## })


## names(post) <- names(post.fsub)[1:10]
        
## post <- rbindlist(post, idcol="Gene_id")            

## ## 2- calculate % of posterior out of region of rejection

## post.out <- per.post(rej=0.2, post, "bj")

## ## caterpillar plot for the extremes

## gt.post <- merge(post.out, gt.Sig, by.x=c("Gene_id", "rSNP"), by.y=c("Gene_id", "tag"))

## gt.post <- gt.post[order(log2_aFC_mean, Gene_id)]

## gt.post[ , index:=1:nrow(gt.post)]

## ggplot(gt.post, aes(index, log2_aFC_mean, ymin=`log2_aFC_2.5%`, ymax=`log2_aFC_97.5%`)) +
##     geom_pointrange(aes(colour=factor(Gene_id))) +coord_flip() + scale_colour_discrete(name="Gene_id") + geom_hline(yintercept= c(-0.1, 0.1), linetype="dashed")

