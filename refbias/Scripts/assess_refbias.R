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


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")



##' # Compare eQTL estimates with and without reference bias correction

## open summaries with no corrections
old <- fread(snakemake@input[['stan_before']])

##old <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/RNAnoGTrsnp/rna99.summary.txt')

## old is rna99 without refrence bias correction run on the same SNPs

## open corrected summaries

new <- fread(snakemake@input[['ref_bias']])

##new <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/results/rna99.BiasCorrec.summary.txt")

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
                          
btrecase.plot(dt=all , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("without", "with"),
              title="eQTL estimates with or without reference panel bias correction"
              )




ai <- fread(snakemake@input[['AI']])
#ai <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/pre_post_AI.txt')

ai[, id:=paste(POS,REF,ALT, sep = ":")]

fsnps.AI <- ai[id %in% unique(all$fSNP.id),]

## plot AI for fsnps run in model

hist(fsnps.AI$AI_post)
summary(fsnps.AI$AI_post)

## select fsnps with the strongest imbalance and look at the difference in eQTL estimates

fsnps.AI[AI_post<0.4 ,]

unique(all[min_AI <0.4& !is.na(fSNP.id),.(Gene_id, n.fsnps.post, n.fsnps.pre, log2_aFC_mean.post, log2_aFC_mean.pre)])




## Look at genes run with GT

gt <- comb.files(path=snakemake@params[['GT_dir']], pattern = "ENSG[0-9]+\\.stan.summary.txt")

## gt <- comb.files(path="/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT", pattern = "ENSG[0-9]+\\.stan.summary.txt")

gt_old <- comb.files(path=snakemake@params[['GT_old']], pattern="ENSG[0-9]+\\.stan.summary.txt")

gt.comb <- merge(gt, gt_old, by=c("Gene_id", "tag"), suffixes=c(".post", ".pre"))

## select associations run with ASE in post (lower number of fSNPs)

gt.ase <- gt.comb[model.post == "trec-ase",]

hist(gt.ase$min_AI)

gt.ase <- add.signif(gt.ase, x1="log2_aFC_null.pre", x2="log2_aFC_null.post", col=c("Without","With") )

cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%", "log2_aFC_null")

## Count significant associations with or without reference panel bias correction
table2 <- gt.ase[,.N,Signif]
names(table2) <- c("Signif", "SNPs")

table2[Signif=="None", color:="#999999"][Signif=="Without", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="With", color:="#0072B2"]

table2[, Signif:=factor(Signif, levels=c("None","Without", "With","Both"))]

setkey(table2,Signif)

btrecase.plot(dt=gt.ase[Signif != "None"] , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without reference panel bias correction"
              )+
    annotation_custom(tableGrob(table2[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table2$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)


btrecase.plot(dt=gt.ase , x1=paste0(cols, ".pre"),
              x2=paste0(cols, ".post"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without reference panel bias correction"
              )+
    annotation_custom(tableGrob(table2[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table2$color, rep("black", 4)))))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)

## Look at extreme examples

sig.post <- gt.ase[Signif=="With",]
sig.post[, dif := abs(`log2_aFC_mean.pre` - `log2_aFC_mean.post`)/abs(`log2_aFC_mean.post`)]
sig.post[dif >=0.5, bias:="big"][dif<0.5, bias:="small"]


## In this case the 95% CI has been reduced after correction plus an effect on the estimate
sig.post[log2_aFC_mean.pre < -1.8,]


## look at effect size of rsnp stratified by being fSNP

## start with gt, look for significant eQTL in close LD with a fSNP. Either the fSNP is eQTL or tagged at 0.9

## Get fSNPs and their target gene

vcf.dna <- snakemake@input[['vcf_GT']]
##vcf.dna='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz'

dna <- vcf_w(vcf.dna) 

## Get fsnps and their host gene 
##fsnps <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt")

fsnps <- fread(snakemake@input[['eSNPs']])

## select fsnps from p.dna
fsnps[, id:=paste(POS,REF,ALT, sep=":")]


## for each fSNP in dna I need the target gene

gene.fsnp <- merge(dna[, c("CHROM","POS", "REF","ALT","id"), with=F], fsnps, by=c("CHROM","POS", "REF","ALT","id"))
    

## for each fsnp run, get the tagging SNP and EAF

gt.tags <- comb.files(path=snakemake@params[['GT_dir']],pattern="eqtl.tags.lookup.txt")

## gt.tags <- comb.files(path = '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/GT', pattern="eqtl.tags.lookup.txt")

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

## prepare tables to add to plot

## get all combinations
g <- data.table(expand.grid(unique(gt.trec$Signif), unique(gt.trec$source)))[order(-Var2),]


tables3 <- lapply(1:nrow(g), function(i) {
    dt <- gt.trec[Signif==g[i,Var1] & source == g[i,Var2] ,.N,.(fSNP,sign(log2mean.eaf))]
    dt[, sign:=as.character(sign)][sign=="-1", sign:="<0"][sign=="1", sign:=">0"]
    setkey(dt, fSNP,sign)
    setnames(dt, c("sign", "N"), c("Effect", "SNPs"))
    return(dt)
})

pal <- brewer.pal(n = 8, name = "Set1")[c(4,3)]

gl <- lapply(tables3, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 7,padding = unit(c(1.5, 1,1,1), "mm"), core=list(fg_params=list(col=c( rep(pal, each=2), rep("black", 8)))))))   

gt.trec <- gt.trec[order(match(source, unique(g$Var2))),]

## make sure names in dt are the same as in gt.trec (Signif and source)
dt <- data.table(Signif=g$Var1, source=g$Var2, grob = gl )

ggplot(gt.trec, aes(log2mean.eaf, color= fSNP)) +
    xlab("eQTL-effect") +
    ylab("Density") +
    geom_density()+
    scale_color_manual(values=pal) +
    theme_bw() +
    geom_vline(xintercept= 0, linetype = "dashed",  color = "grey")+
    facet_grid(rows= vars(source), cols=vars(Signif), scales = "free_x")+
    theme(strip.background=element_rect(fill="white"))+
    geom_custom(data=dt, aes(grob=grob), x = 0.75, y = 0.75)





## all[is.na(fSNP.id), fSNP:="no"][!is.na(fSNP.id), fSNP:= "yes"]


## ## add trec as control: no reference panel bias expected

## trec <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.stan.summary.txt")

## ## tags in trec and all may be different

## ##trec.tags <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.tags.lookup.txt")

## trec.tags <- comb.files(path=snakemake@params[['trec']], pattern="trec.tags.lookup.txt")

## trec.tags.run <- merge(trec.tags, trec[,.(Gene_id,tag)], by= c("Gene_id","tag"))

## ## get tags run in rna99


## ## rna.tags <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/rna99/',pattern="eqtl.tags.lookup.txt")

## rna.tags <- comb.files(path=snakemake@params[['Imp_new']],pattern="eqtl.tags.lookup.txt")

## rna.tags.run <- merge(rna.tags, all[,.(Gene_id,tag)], by= c("Gene_id","tag"))

## ## combine trec and rna (all), allowing for different tags
## rna.trec <- comb.sum(x=trec, y=gt.ase, z=all, s=c(".trec", ".gt", ".rna"), trec.tags.run, rna.tags.run)

## ## make into long format for plotting

## ## first select relevant columns

## pre.cols <- c("log2_aFC_mean",  "log2_aFC_null", "Gene_id", "tag.rna", "tag.gt", "fSNP.rna")
## suf.cols <- c("trec", "post.gt","pre.gt","post.rna","pre.rna")


## rna.trec.l <- lapply(suf.cols, function(i){
##     dt <- rna.trec[, unlist(lapply(pre.cols, function(i) grep(i, names(rna.trec)))), with=F]
##     tmp <- dt[, c(grep(i, names(dt), value=T), pre.cols[3:length(pre.cols)]), with=F]
##     return(tmp)
## })
## names(rna.trec.l) <- suf.cols

## ## from list to long and formatting
## rna.trec.l  <- rbindlist(lapply(rna.trec.l, function(i) setNames(i, pre.cols )), idcol="Source")
## rna.trec.l[, Signif:="No significant"][log2_aFC_null=="no", Signif:="Significant"][is.na(fSNP.rna), fSNP.rna:="no"]
## setnames(rna.trec.l, "fSNP.rna", "fSNP")

## ## make tables for plot
## g <- data.table(expand.grid(unique(rna.trec.l$Signif), unique(rna.trec.l$Source)))[order(Var2),]

## tables3 <- lapply(1:nrow(g), function(i) {
##     dt <- rna.trec.l[Signif==g[i,Var1] & Source == g[i,Var2] ,.N,.(fSNP,sign(log2_aFC_mean))]    
##     dt[, sign:=as.character(sign)][sign=="-1", sign:="<0"][sign=="1", sign:=">0"]
##     setkey(dt, fSNP,sign)
##     setnames(dt, c("sign", "N"), c("Effect", "SNPs"))
##     return(dt)
## })

## pal <- brewer.pal(n = 8, name = "Set1")[c(4,3)]

## gl <- lapply(tables3, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 7,padding = unit(c(1.5, 1,1,1), "mm"), core=list(fg_params=list(col=c( rep(pal, each=2), rep("black", 8)))))))   

## dt <- data.table(Signif=g$Var1, source=g$Var2, grob = gl )


## lapply(rna.trec.l, function(i) summary

## ref.bias <- ggplot(rna.trec.l, aes(log2_aFC_mean, color= fSNP)) +
##     xlab("eQTL-effect") +
##     ylab("Density") +
##     geom_density()+
##     scale_color_manual(values=pal) +
##     theme_bw() +
##     geom_vline(xintercept= 0, linetype = "dashed",  color = "grey")+
##     facet_grid(rows= vars(Source), cols=vars(Signif), scales = "free_x")+
##     theme(strip.background=element_rect(fill="white"))+
##     geom_custom(data=dt, aes(grob=grob), x = 0.75, y = 0.75)
