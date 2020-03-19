#' ---
#' title: Set threshold for making significant calls with Bayesian model
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
library(grid)


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")

##' # Compare associations using  Bayesian and frequentist output for associations

##' Use frequentist 5%FDR threshold


## trec <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.stan.summary.txt")
## freq <- comb.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/NB", "eqtl.NegBinom.txt")



trec <- comb.files(path=snakemake@params[['trec']], pattern="trec.stan.summary.txt")
freq <- comb.files(path=snakemake@params[['freq']], pattern="eqtl.NegBinom.txt")

## add BH correction to freq

setkey(freq, `Pr(>|z|)`)
freq[,p.adj:= p.adjust(`Pr(>|z|)`,method = "BH")]
setkey(freq, p.adj)

##  express estimate in log2 scale in the same units as log(a_FC)
freq[, log2_est:=2*Estimate/log(2)]

## merge trec and freq and save

freq.trec <- merge(trec, freq, by=c("Gene_id", "tag"))

saveRDS(freq.trec, "/home/ev250/newshare/elena/trecase/objects/NB.freq.trec.chr22.rds")



## Add null column for freq based on 5% FDR
freq[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]

## Bayesian trec: using normal approximation calculate proportion of posterior out of rejection zone. If the mean of the posterior is within the rejection zone (-r,r) I set the posterior out of the rejection zone as 0% as I dont want to call any of those associations significant. If the rejection zone is narrow I could have a high % of the posterior out of the zone. The I coded a variable, null.rej "yes" if the % of the posterior out of the rejection zone is below a threshold I define (post.level) and "no" otherwise.

rej<-seq(0,0.8,0.2)
post.level=0.9

post.dt <- rej.recode(a=rej, b=trec,c=post.level)

setkey(post.dt, Gene_id,tag)

## compare multiple testing correction with calling singnificance using 95% CI

ci.mult(rej,post.dt)

## merge post.l and freq

all <- merge(post.dt,freq, by=c("Gene_id", "tag"))
setkey(all, p.adj)

all <- add.signif(all, x1="null.rej", x2="null.fdr5", col=c("trec","freq") )

cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%")

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

##' #Plot thresholds for multiple testing correction in Bayesian negative binomial  and compare with frequentist model at 5%FDR
##'
##'

table <- all[,.N,.(Signif, rej.level)]
setnames(table, "N", "SNPs")

table[Signif=="None", color:="#999999"][Signif=="trec", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="freq", color:="#0072B2"]

table[, Signif:=factor(Signif, levels=c("None","trec", "freq","Both"))]

setkey(table,Signif)

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)


#+ fig.width= 10, fig.height=21
btrecase.plot(dt=all[Signif != "None" ,] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_est",3), "null.fdr5"),
              s=nrow(all[Signif != "None" ,]),
              xl="eQTL effect Bayesian trec",
              yl="eQTL effect freqentist",
              col=c("trec", "freq"),
              title="NB eQTL estimates based on 5% FDR or\n 90% posterior our of the indicated rejection zone",
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
    theme(strip.background=element_rect(fill="white"),strip.text.y = element_text(size = 14) )+
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.8, y = 0.3)

##' # Extend analysis using a range of FDR and a range of %of posterior and calculate Jaccard index (intersection over union)


## freq: add fdr at 0.1, 0.01, 0.001
## post.per: 0.5-0.9

fdr <- c(0.001, 0.01, 0.1)

for (i in fdr){
    freq[,eval(paste0("null.fdr",i*100)):= "yes"][p.adj<=i, paste0("null.fdr",i*100):="no"]
}

all2 <- merge(post.dt,freq, by=c("Gene_id", "tag"))
setkey(all2, p.adj)

## Prepare table for Jaccard index

post.per <- seq(0.5,0.9, 0.1)

JI <- rbindlist(lapply(post.per, function(i) {
    rbindlist(lapply(c(fdr, 0.05), function(j){
        tmp <- all2[,.N , .(get(paste0("null.fdr", j*100)), rej.level, post.out >=i)]
        tmp <- tmp[!(post.out == FALSE & get =="yes"),]
        setkey(tmp, rej.level)
        num <- tmp[get == "no" & post.out ==TRUE, N, rej.level]
        den <- tmp[, .(tot=sum(N)), rej.level]
        ind <- merge(num, den, by="rej.level")
        ind[, JI:=N/tot]
        ind[, fdr:= j]
        ind[, post.per:=i]
    }))
    }))
    
 #+ fig.width= 25, fig.height=10   
ggplot(JI, aes(fdr*100, JI, color=as.factor(post.per))) + geom_line() + geom_point() + facet_grid(.~rej.level) + ggtitle("Jaccard index by FDR(%), \n rejection zone and Posterior in rej. zone (%)") + xlab("FDR") + theme(axis.text.x = element_text(angle = 45)) + scale_x_log10(breaks = fdr*100) + labs(color = "Post (%)")

   

JI[JI==max(JI),]

## add CI freq cols for plot

all2[ , freq.cilow:=log2_est-1.96*`Std. Error`][, freq.cih:=log2_est+1.96*`Std. Error`]

all2 <- add.signif(all2, x1="null.rej", x2="null.fdr10", col=c("Btrec","Freq") )


sig.cols <- c("None"="#999999","Btrec"="yellow3", "Both"="#D55E00", "Freq"="#0072B2")
tab <- tab2bplot(all2[rej.level==0.4,], colors=sig.cols)


#+ fig.width= 6, fig.height=5
btrecase.plot(dt=all2[rej.level==0.4 & !(post.out<=0.9 & null.fdr10 =="yes")  ,] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_est",3), "null.fdr10"),
              s=nrow(all2),
              xl="eQTL effect Bayesian",
              yl="eQTL effect freqentist",
              col=c("Btrec", "Freq"),
              title="NB freq with 10% FDR and Bayesian based \n on 90% posterior out rejection zone = 0.4"
              ) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=8, ymax=15) +
    annotation_custom(grobTree(textGrob(paste("Jaccard index =", round(max(JI$JI),1)), x=0.25, y=0.75, gp=gpar(fontsize=10))))


## Explore outlier in eQTl freq vs Bayesian

outl <- all[Signif =="freq" & log2_est == min(log2_est) & rej.level==0,]

## no ASE, need to prepare input

counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'
counts.g <-  fread(cmd=paste("grep -e gene_id -e ",outl$Gene_id,counts.f), header=TRUE)
counts.g <- counts.g[,2:ncol(counts.g),with=F] ## removes gene_id
genotypes <- readRDS("/rds/project/cew54/rds-cew54-wallace-share/elena/trecase/objects/chr22.GT.geuvardis.rds")
genotypes[,id:=paste(POS, REF,ALT, sep=":")]

genotypes <- genotypes[id == outl$tag,]
genotypes <- genotypes[,sapply(names(counts.g), function(i) grep(i, names(genotypes), value=T)), with=F]
setnames(genotypes, names(genotypes), gsub("_GT$", "", names(genotypes)))

## merge counts.g and genotypes

c.g <- data.table(t(rbind(genotypes, counts.g)), keep.rownames=T)
names(c.g) <- c("sample", "gt", "counts")
c.g[, Genotype:=as.factor(abs(gt))]

##' # Inputs for extreme eQTL in freq but null in Bayesian
c.g[,.N, Genotype]
c.g[, .(`counts summary`=summary(counts)), Genotype]
ggplot(as.data.table(c.g), aes(y=counts, x=Genotype, color=Genotype)) + geom_jitter(width = 0.1, shape=1) + theme(legend.position = "none") + ggtitle("Inputs for extreme eQTL in freq vs Bayesian NB")

glm.nb(counts ~ abs(gt), data=c.g)

##########################################################################################################
## REpeat analysis using dseq2 NB model (run by Wei-Yu)

## dseq <- rbindlist(lapply(list.files("/mrc-bsu/scratch/wyl37/ElenaData/RunNBmodel", pattern="RunNBmodelbatch[0-9]+_chr22.nbmodelRes.csv", full.names=T), fread))

dseq <- rbindlist(lapply(snakemake@input[['dseq']], fread))

dseq[, SNP:=NULL]

## add log2_aFC to dseq

dseq[, log2_aFC:=log2FoldChange*2]

## add BH correction to dseq

setkey(dseq, pvalue)
dseq[,p.adj:= p.adjust(pvalue,method = "BH")]
setkey(dseq, p.adj)

## Add null column for dseq based on 5% FDR and  exclude "ENSG00000211664"
dseq[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]
dseq <- dseq[Gene_id != "ENSG00000211664",]

## Compare dseq with freq,

d.f <- merge(dseq ,freq, by=c("Gene_id", "tag"), suffixes=c(".dseq", ".freq"))
d.f <- add.signif(d.f, x1="null.fdr5.dseq", x2="null.fdr5.freq", col=c("dseq","freq") )

d.f.tab <- tab2bplot(d.f, colors=setNames( c("#999999", "yellow3", "#0072B2", "#D55E00"), c("None", "dseq","freq", "Both")))

btrecase.plot(dt=d.f[Signif != "None",] , x1=c(rep("log2_aFC",3),"null.fdr5.dseq") ,
              x2=c(rep("log2_est",3), "null.fdr5.freq"),
              xl="eQTL effect NB-frequentist",
              yl="eQTL effect NB-DSeq2",
              col=c("dseq", "freq"),
              title="eQTL estimates using NB-freq or NB-DSeq2 model\n at 5% FDR"
              ) +
    annotation_custom(tableGrob(d.f.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(d.f.tab$color, rep("black", 4)))))), xmin=.8, xmax=15, ymin=-15, ymax=-5)


###########################################################################################################
###############  Compare Bayesian NB with Dseq2 NB model ###############

## merge post.l and dseq

btrec.dsq <- merge(post.dt,dseq, by=c("Gene_id", "tag"))
setkey(btrec.dsq, p.adj)

btrec.dsq <- add.signif(btrec.dsq, x1="null.rej", x2="null.fdr5", col=c("trec","dseq") )

##' #Plot thresholds for multiple testing correction in Bayesian negative binomial  and compare with NB-Dseq2  model at 5%FDR


table <- btrec.dsq[,.N,.(Signif, rej.level)]
setnames(table, "N", "SNPs")

table[Signif=="None", color:="#999999"][Signif=="trec", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="dseq", color:="#0072B2"]

table[, Signif:=factor(Signif, levels=c("None","trec", "dseq","Both"))]

setkey(table,Signif)

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)


#+ fig.width= 10, fig.height=21
btrecase.plot(dt=btrec.dsq[Signif != "None" ,] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_aFC",3), "null.fdr5"),
              s=nrow(btrec.dsq[Signif != "None" ,]),
              xl="eQTL effect Bayesian NB",
              yl="eQTL effect Dseq2",
              col=c("trec", "freq"),
              title="NB eQTL estimates based on 5% FDR or\n 90% posterior out of the indicated rejection zone",
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
    theme(strip.background=element_rect(fill="white"), strip.text.y = element_text(size = 14))+
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.1, y = 0.8)


##' # Extend analysis using a range of FDR and a range of %of posterior and calculate Jaccard index (intersection over union)

for (i in fdr){
    dseq[,eval(paste0("null.fdr",i*100)):= "yes"][p.adj<=i, paste0("null.fdr",i*100):="no"]
}

btrec.dsq2 <- merge(post.dt,dseq, by=c("Gene_id", "tag"))
setkey(btrec.dsq2, p.adj)

## Prepare table for Jaccard index (Both/Total excluding None)

JI <- rbindlist(lapply(post.per, function(i) {
    rbindlist(lapply(c(fdr, 0.05), function(j){
        tmp <- btrec.dsq2[,.N , .(get(paste0("null.fdr", j*100)), rej.level, post.out >=i)]
        tmp <- tmp[!(post.out == FALSE & get =="yes"),]
        setkey(tmp, rej.level)
        num <- tmp[get == "no" & post.out ==TRUE, N, rej.level]
        den <- tmp[, .(tot=sum(N)), rej.level]
        ind <- merge(num, den, by="rej.level")
        ind[, JI:=N/tot]
        ind[, fdr:= j]
        ind[, post.per:=i]
    }))
    }))
    
 #+ fig.width= 25, fig.height=10   
ggplot(JI, aes(fdr*100, JI, color=as.factor(post.per))) + geom_line() + geom_point() + facet_grid(.~rej.level) + ggtitle("Jaccard index by FDR, \n rejection zone and posterior out of rej. zone (%)") + xlab("FDR (%)") + theme(axis.text.x = element_text(angle = 45)) + scale_x_log10(breaks = fdr*100) + labs(color = "Post (%)")

   

JI[JI==max(JI),]


btrec.dsq2 <- add.signif(btrec.dsq2, x1="null.rej", x2="null.fdr10", col=c("Btrec","Dseq2") )


sig.cols <- c("None"="#999999","Btrec"="yellow3", "Both"="#D55E00", "Dseq2"="#0072B2")
tab <- tab2bplot(btrec.dsq2[rej.level==0.4,], colors=sig.cols)


#+ fig.width= 6, fig.height=5
btrecase.plot(dt=btrec.dsq2[rej.level==0.4 & !(post.out<=0.9 & null.fdr10 =="yes")  ,] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_aFC",3), "null.fdr10"),
              xl="eQTL effect Bayesian trec",
              yl="eQTL effect Dseq2",
              col=c("Btrec", "Dseq2"),
              title="NB Dseq2 with 10% FDR and Bayesian based\n on 90% posterior out of rejection zone = 0.4"
              ) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=8, ymax=15) +
    annotation_custom(grobTree(textGrob(paste("Jaccard index =", round(max(JI$JI),2)), x=0.2, y=0.75, gp=gpar(fontsize=10))))


###############################################################################################################################


##### Extend analysis to Btrec vs Btrecase-GT and Btrecase-noGT

## Bayesian models: using normal approximation calculate proportion of posterior out of rejection zone. If the mean of the posterior is within the rejection zone (-r,r) I set the posterior out of the rejection zone as 0% as I dont want to call any of those associations significant. If the rejection zone is narrow I could have a high % of the posterior out of the zone. The I coded a variable, null.rej "yes" if the % of the posterior out of the rejection zone is below a threshold I define (post.level) and "no" otherwise.

## Get btrecase-GT

## btrecase <- comb.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/GT', "^refbias\\.ENSG[0-9]+.*stan.summary.txt")

btrecase <- comb.files(snakemake@params[['btrecase']], "^refbias\\.ENSG[0-9]+.*stan.summary.txt")


## Get rna with tags matched to gt
## gt.rna <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/results/refbias.gt.rna.txt')

gt.rna <- fread(snakemake@input[['gt_rna_sum']])

## reduce rejection zone
rej<-seq(0,0.5,0.1)
post.level=0.9

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

cols <- mapply(function(x,i,j) sapply(x,function(k) grep(paste0(k, i, "$"), names(j), value=T)),              
               i=c("", "", ".rna"),
               j=list(trec, btrecase, gt.rna),
               MoreArgs=list(x=c("log2_aFC_mean" ,"log2_aFC_sd", "tag")),
               SIMPLIFY=F,
               USE.NAMES=F
               )


post.btrecase <- mapply(function(a,b,cols,c,t) {
    dt <- rej.recode(a, b,cols, c)
    setkeyv(dt, c("Gene_id",t))  
    return(dt)},
    b=list(nb=trec, gt=btrecase, rna=gt.rna),
    cols= lapply(cols, function(i) i[1:2]),
    t=sapply(cols, function(i) i[3]),
    MoreArgs=list(a=rej, c=post.level),
    SIMPLIFY=F)
                  
   
## Merge btrec with post.btrecase$gt

trec.gt <- merge(post.btrecase$nb, post.btrecase$gt, by=c("Gene_id", "tag", "rej.level"), suffixes=c(".trec", ".gt"))
trec.gt <- add.signif(trec.gt, x1="null.rej.trec", x2="null.rej.gt", col=c("trec","btrecase") )

## plot hist for ASE assoc comparing trec-ase with trec

##' #Plot thresholds for multiple testing correction in Bayesian negative binomial (trec) and Bayesian NB-ASE (btrecase) only considering associations with ASE information in btrecase

trec.gtase <- trec.gt[model.gt=="trec-ase",]

table <- tab2bplot(dt=trec.gtase, var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", btrecase="#0072B2", Both= "#D55E00"))

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)


## So far I was comparing effects using 95% CI:

btrecase <- add.null(dt=btrecase)
trec.refbias <- merge(trec, btrecase[model =="trec-ase",], by=c("Gene_id", "tag"), suffixes=c(".trec",".trec-ase"))
trec.refbias <- add.signif(trec.refbias, x2="null.95", x1= "log2_aFC_null" , col=c("trec-ase", "trec"))
 
trec.rb.tab <- tab2bplot(trec.refbias, colors=c(None="#999999", `trec-ase`="yellow3", trec="#0072B2", Both= "#D55E00"))

#+ fig.width= 6, fig.height=5
btrecase.plot(dt=trec.refbias[Signif != "None",] ,
              x2=c(paste0(c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%"), ".trec-ase"), "null.95"),
              x1=c(paste0(c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%"), ".trec"),  "log2_aFC_null" ),
              yl="eQTL effect Bayesian trec-ase",
              xl="eQTL effect Bayesian trec",
              col=c("trec-ase", "trec"),
              title="eQTL estimates using 95% CI"
              ) + 
    annotation_custom(tableGrob(trec.rb.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(trec.rb.tab$color, rep("black", 4)))))), xmin=-2, xmax=-1, ymin=1, ymax=2.5)



p <- btrecase.plot(dt=trec.gtase[Signif != "None" ,] , x1=c(rep("log2_aFC_mean.trec",3),"null.rej.trec") ,
                   x2=c(rep("log2_aFC_mean.gt",3), "null.rej.gt"),
                   s=50000,
                   xl="eQTL effect Bayesian trec",
                   yl="eQTL effect Bayesian trec-ase",
                   col=c("trec", "btrecase"),
                   title="eQTL estimates based on 90% posterior out of the indicated\n rejection zone for associations with ASE",
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
                   ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
    theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14)) +
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)

#+ fig.width= 10, fig.height=21
p

#' # Need to change rejection zone as eQTL estimates with ASE are of small magnitude than trec only

## reduce rejection zone
rej<-seq(0,0.25,0.05)
post.level=0.9

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

cols <- mapply(function(x,i,j) sapply(x,function(k) grep(paste0(k, i, "$"), names(j), value=T)),              
               i=c("", "", ".rna"),
               j=list(trec, btrecase, gt.rna[model.gt == "trec-ase",]),
               MoreArgs=list(x=c("log2_aFC_mean" ,"log2_aFC_sd", "tag")),
               SIMPLIFY=F,
               USE.NAMES=F
               )


post.btrecase <- mapply(function(a,b,cols,c,t) {
    dt <- rej.recode(a, b,cols, c)
    setkeyv(dt, c("Gene_id",t))  
    return(dt)},
    b=list(nb=trec, gt=btrecase, rna=gt.rna),
    cols= lapply(cols, function(i) i[1:2]),
    t=sapply(cols, function(i) i[3]),
    MoreArgs=list(a=rej, c=post.level),
    SIMPLIFY=F)
                  
   
## Merge post.btrec$nb with post.btrecase$gt

trec.gt <- merge(post.btrecase$nb, post.btrecase$gt[model == "trec-ase",], by=c("Gene_id", "tag", "rej.level"), suffixes=c(".trec", ".gt"))
trec.gt <- add.signif(trec.gt, x1="null.rej.trec", x2="null.rej.gt", col=c("trec","btrecase") )

table <- tab2bplot(dt=trec.gt, var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", btrecase="#0072B2", Both= "#D55E00"))

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)


p <- btrecase.plot(dt=trec.gt[Signif != "None" ,] , x1=c(rep("log2_aFC_mean.trec",3),"null.rej.trec") ,
                   x2=c(rep("log2_aFC_mean.gt",3), "null.rej.gt"),
                   xl="eQTL effect Bayesian trec",
                   s=nrow(trec.gt),
                   yl="eQTL effect Btrecase",
                   col=c("trec", "btrecase"),
                   title="eQTL estimates based on 90% posterior out of the indicated\n rejection zone for associations with ASE",
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
                   ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
    theme(strip.background=element_rect(fill="white"),strip.text.y = element_text(size = 14))+
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)

##+ fig.width= 12, fig.height=21
p


##' #Plot thresholds for multiple testing correction in Bayesian trecase with and without GT


gt.rna2 <- merge(post.btrecase$gt, post.btrecase$rna, by.x=c("Gene_id", "tag", "rej.level"), by.y=c("Gene_id", "tag.gt", "rej.level"), suffixes=c(".gt", ".rna"))
                
gt.rna2 <- add.signif(gt.rna2, x1="null.rej.gt", x2="null.rej.rna", col=c("Ob-GT","Hidden-GT") )

table <- tab2bplot(dt=gt.rna2, var= c("Signif", "rej.level"), colors=c(None="#999999", `Ob-GT`="yellow3", `Hidden-GT`="#0072B2", Both= "#D55E00"))

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)

## remove op.dir for plotting without CI

gt.rna2[,op.dir:=NULL]

## Plot associations using 99%CI to compare:

gt.rna<- add.signif(gt.rna, x1="null.99.gt", x2="null.99.rna", col=c("obs-GT","hidden-GT"))
gt.rna.tab <- tab2bplot(dt=gt.rna, colors= c(None="#999999", `obs-GT`="yellow3", `hidden-GT`="#0072B2", Both= "#D55E00"))

btrecase.plot(dt=gt.rna[Signif != "None"  ,],
                        x1=c('log2_aFC_mean.gt', 'log2_aFC_2.5%.gt','log2_aFC_97.5%.gt', paste0('null.',99,'.gt')),
                        x2= paste0(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', paste0('null.', 99)), ".rna"),
                        xl='eQTL-effect (observed GT)', yl='eQTL-effect (hidden GT)',
                        col=c("obs-GT","hidden-GT"),axis.title=12, axis.text=10,
                        legend.title=12, legend.text=10,
                        legend.symbol=4, point.size=3 ,
                        title=paste0("Observed (",99 ," %CI) vs hidden genotypes (",99," %CI)") , title.size=12) +
    annotation_custom(tableGrob(gt.rna.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(gt.rna.tab$color, rep("black", 4)))))), xmin=-1.3, xmax=-.5, ymin=0.5, ymax=1)


p <- btrecase.plot(dt=gt.rna2[Signif != "None" ,] , x1=c(rep("log2_aFC_mean",3),"null.rej.gt") ,
                   x2=c(rep("log2_aFC_mean.rna",3), "null.rej.rna"),
                   s=nrow(gt.rna2[Signif != "None" ,]), 
                   xl="eQTL effect Obs-GT",
                   yl="eQTL effect Hidden-GT",
                   col=c("Ob-GT", "Hidden-GT"),
                   title=paste0("eQTL estimates based on ", post.level*100 , "% posterior\n out of the indicated rejection zone"),
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
                   ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
    theme(strip.background=element_rect(fill="white")  ,strip.text.y = element_text(size = 14))+
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)



#+ fig.width= 12, fig.height=21
p



##' # Compare with % posterior = 0.5 keeping same rejection zone

post.level  <- 0.5

post.btrecase <- mapply(function(a,b,cols,c,t) {
    dt <- rej.recode(a, b,cols, c)
    setkeyv(dt, c("Gene_id",t))  
    return(dt)},
    b=list(nb=trec, gt=btrecase, rna=gt.rna),
    cols= lapply(cols, function(i) i[1:2]),
    t=sapply(cols, function(i) i[3]),
    MoreArgs=list(a=rej, c=post.level),
    SIMPLIFY=F)
   
gt.rna2 <- merge(post.btrecase$gt, post.btrecase$rna, by.x=c("Gene_id", "tag", "rej.level"), by.y=c("Gene_id", "tag.gt", "rej.level"), suffixes=c(".gt", ".rna"))
                
gt.rna2 <- add.signif(gt.rna2, x1="null.rej.gt", x2="null.rej.rna", col=c("Ob-GT","Hidden-GT") )

table <- tab2bplot(dt=gt.rna2, var= c("Signif", "rej.level"), colors=c(None="#999999", `Ob-GT`="yellow3", `Hidden-GT`="#0072B2", Both= "#D55E00"))

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)

## remove op.dir for plotting without CI

gt.rna2[,op.dir:=NULL]

p <- btrecase.plot(dt=gt.rna2[Signif != "None" ,] , x1=c(rep("log2_aFC_mean",3),"null.rej.gt") ,
                   x2=c(rep("log2_aFC_mean.rna",3), "null.rej.rna"),
                   s=nrow(gt.rna2[Signif != "None" ,]), 
                   xl="eQTL effect Obs-GT",
                   yl="eQTL effect Hidden-GT",
                   col=c("Ob-GT", "Hidden-GT"),
                   title=paste0("eQTL estimates based on ", post.level*100 , "% posterior\n out of the indicated rejection zone"),
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
                   ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
    theme(strip.background=element_rect(fill="white")  ,strip.text.y = element_text(size = 14))+
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)



#+ fig.width= 12, fig.height=21
p





##' # Extend analysis using a range of  % of posterior out of rejection zone and calculate Jaccard index (intersection over union)

post.per <- seq(0.5,0.9, 0.1)

post.cols <- grep("post.out", names(gt.rna2), value=T)

JI <- rbindlist(lapply(post.per, function(i) {
    tmp <- gt.rna2[,.N , .(rej.level, get(post.cols[1]) >=i, get(post.cols[2]) >=i)]
    ## exclude null associations in both
    tmp <- tmp[ !( get == FALSE & get.1 == FALSE),]
    setkey(tmp, rej.level)
    num <- tmp[get == TRUE & get.1 ==TRUE, N, rej.level]
    den <- tmp[, .(tot=sum(N)), rej.level]
    ind <- merge(num, den, by="rej.level")
    ind[, JI:=N/tot]
    ind[, post.per:=i]
    }))

    
 #+ fig.width= 25, fig.height=10   
ggplot(JI, aes(rej.level, JI, color=as.factor(post.per))) +
    geom_line() +
    geom_point() +
    ggtitle("Jaccard index by posterior out of rej. zone (%)") +
    xlab("Rejection level") +
    labs(color = "Post (%)")
