#' ---
#' title: Set threshold for making significant calls with Bayesian model using same prior with/without GT
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
library(biomaRt)
library(mixtools)


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")

##' # Compare associations using  Bayesian and frequentist output for associations

##' Use frequentist 5%FDR threshold


## trec <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.stan.summary.txt")
## freq <- comb.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/NB", "eqtl.NegBinom.txt")


trec <- comb.files(path=snakemake@params[['trec']], pattern="trec.stan.summary.txt")
freq <- comb.files(path=snakemake@params[['freq']], pattern="eqtl.NegBinom.txt")

## add gene distance to trec

## gene.coord <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt")

gene.coord <- fread(snakemake@input[['geneStEnd']]) 

## select genes in chrom 22

gt22 <- gene.coord[chrom==22,]

## add tag distance to gene (closest to start or end)

trec <- gene.d(trec, gt22[, .(gene_id, start,end,chrom)])


## add BH correction to freq

setkey(freq, `Pr(>|z|)`)
freq[,p.adj:= p.adjust(`Pr(>|z|)`,method = "BH")]
setkey(freq, p.adj)

##  express estimate in log2 scale in the same units as log(a_FC)
freq[, log2_est:=2*Estimate/log(2)]

## merge trec and freq 

freq.trec <- merge(trec, freq, by=c("Gene_id", "tag"))

## Add null column for freq based on 5% FDR
freq[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]

## Bayesian trec: using normal approximation calculate proportion of
## posterior out of rejection zone. If the mean of the posterior is
## within the rejection zone (-r,r) I set the posterior out of the
## rejection zone as 0% as I dont want to call any of those
## associations significant. If the rejection zone is narrow I could
## have a high % of the posterior out of the zone. The I coded a
## variable, null.rej "yes" if the % of the posterior out of the
## rejection zone is below a threshold I define (post.level) and "no"
## otherwise.


post.level=0.95
rej<-seq(0,0.8,0.2)

post.dt <- rej.recode(a=rej, b=trec,c=post.level)

setkey(post.dt, Gene_id,tag)

## compare multiple testing correction with calling singnificance using 95% CI

ci.mult(rej,post.dt)

## merge post.dt and freq

all <- merge(post.dt,freq, by=c("Gene_id", "tag"))
setkey(all, p.adj)

all <- add.signif(all, x1="null.rej", x2="null.fdr5", col=c("trec","freq") )

cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%")

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

##' #Plot thresholds for multiple testing correction in Bayesian negative binomial  and compare with frequentist model at 5%FDR


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
              title=paste0("NB eQTL estimates based on 5% FDR or\n " , post.level*100, "% posterior out of the indicated rejection zone"),
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
## post.per: 0.5-0.95, 0.99

fdr <- c(0.001, 0.01, 0.1)

for (i in fdr){
    freq[,eval(paste0("null.fdr",i*100)):= "yes"][p.adj<=i, paste0("null.fdr",i*100):="no"]
}

all2 <- merge(post.dt,freq, by=c("Gene_id", "tag"))
setkey(all2, p.adj)

## Prepare table for Jaccard index

post.per <- c(seq(0.5,0.9, 0.1), 0.95, 0.99)

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

   

maxJI=JI[JI==max(JI),]

## add CI freq cols for plot

all2[ , freq.cilow:=log2_est-1.96*`Std. Error`][, freq.cih:=log2_est+1.96*`Std. Error`]

all2 <- add.signif(all2, x1="null.rej", x2="null.fdr10", col=c("Btrec","Freq") )


sig.cols <- c("None"="#999999","Btrec"="yellow3", "Both"="#D55E00", "Freq"="#0072B2")
tab <- tab2bplot(all2[rej.level==maxJI$rej.level,], colors=sig.cols)


#+ fig.width= 6, fig.height=5
btrecase.plot(dt=all2[rej.level==maxJI$rej.level & Signif !="None" & !(post.out<=maxJI$post.per & eval(paste0("null.fdr",maxJI$fdr)) =="yes")  ,] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_est",3), "null.fdr10"),
              #s=nrow(all2),
              xl="eQTL effect Bayesian",
              yl="eQTL effect freqentist",
              col=c("Btrec", "Freq"),
              title=paste0("NB freq with ",100*maxJI$fdr, "% FDR and Bayesian based \n on ", maxJI$post.per*100,"% posterior out of rejection zone = ", maxJI$rej.level )
              ) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=8, ymax=15) +
    annotation_custom(grobTree(textGrob(paste("Jaccard index =", round(maxJI$JI,1)), x=0.22, y=0.75, gp=gpar(fontsize=10))))


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

## remove pval NA

dseq <- dseq[!is.na(pvalue),]

## add log2_aFC to dseq

dseq[, log2_aFC:=log2FoldChange*2]

## add BH correction to dseq

setkey(dseq, pvalue)
dseq[,p.adj:= p.adjust(pvalue,method = "BH")]
setkey(dseq, p.adj)

## Add null column for dseq based on 5% FDR and  exclude "ENSG00000211664"
dseq[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]
dseq <- dseq[Gene_id != "ENSG00000211664",]

## Add gene distance

dseq <- gene.d(dseq, gt22[, .(gene_id, start,end,chrom)], snp="tag")

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
              col=c("trec", "dseq"),
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

   

maxJI <- JI[JI==max(JI),]


btrec.dsq2 <- add.signif(btrec.dsq2, x1="null.rej", x2="null.fdr10", col=c("Btrec","Dseq2") )


sig.cols <- c("None"="#999999","Btrec"="yellow3", "Both"="#D55E00", "Dseq2"="#0072B2")
tab <- tab2bplot(btrec.dsq2[rej.level==maxJI$rej.level,], colors=sig.cols)


#+ fig.width= 6, fig.height=5
btrecase.plot(dt=btrec.dsq2[rej.level==maxJI$rej.level& Signif !="None" & !(post.out<=maxJI$post.per & eval(paste0("null.fdr",maxJI$fdr)) =="yes")  ,] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_aFC",3), "null.fdr10"),
              xl="eQTL effect Bayesian trec",
              yl="eQTL effect Dseq2",
              col=c("Btrec", "Dseq2"),
              title=paste0("DEseq ",100*maxJI$fdr, "% FDR and Bayesian based \n on ", maxJI$post.per*100,"% posterior out of rejection zone = ", maxJI$rej.level ),
              ) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=5, ymax=8) +
    annotation_custom(grobTree(textGrob(paste("Jaccard index =", round(max(JI$JI),2)), x=0.2, y=0.7, gp=gpar(fontsize=10))))



###############################################################################################################################


##### Extend analysis to Btrec vs Btrecase-GT and Btrecase-noGT

## Bayesian models: using normal approximation calculate proportion of
## posterior out of rejection zone. If the mean of the posterior is
## within the rejection zone (-r,r) I set the posterior out of the
## rejection zone as 0% as I dont want to call any of those
## associations significant. If the rejection zone is narrow I could
## have a high % of the posterior out of the zone. The I coded a
## variable, null.rej "yes" if the % of the posterior out of the
## rejection zone is below a threshold I define (post.level) and "no"
## otherwise.

## Get btrecase-GT

## btrecase <- comb.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/GT', "^refbias\\.ENSG[0-9]+.*stan.summary.txt")

btrecase <- comb.files(snakemake@params[['btrecase']], "^refbias\\.ENSG[0-9]+.*stan.summary.txt")

btrecase <- gene.d(btrecase, gt22[, .(gene_id, start,end,chrom)])


## Get rna with tags matched to gt
## gt.rna <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/results/refbias.gt.rna.p054.txt')

gt.rna <- fread(snakemake@input[['gt_rna_p054sum']])

##gt.rna <- gene.d(gt.rna, gt22[, .(gene_id, start,end,chrom)], snp="tag.gt")

## reduce rejection zone
rej<-seq(0,0.5,0.1)
post.level=0.95

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

cols <- mapply(function(x,i,j) sapply(x,function(k) grep(paste0(k, i, "$"), names(j), value=T)),              
               i=c("", "", ".rna"),
               j=list(trec, btrecase, gt.rna),
               MoreArgs=list(x=c("log2_aFC_mean" ,"log2_aFC_sd", "tag")),
               SIMPLIFY=F,
               USE.NAMES=F
               )


post.btrecase95.1 <- mapply(function(a,b,cols,c,t) {
    dt <- rej.recode(a, b,cols, c)
    setkeyv(dt, c("Gene_id",t))  
    return(dt)},
    b=list(nb=trec, gt=btrecase, rna=gt.rna),
    cols= lapply(cols, function(i) i[1:2]),
    t=sapply(cols, function(i) i[3]),
    MoreArgs=list(a=rej, c=post.level),
    SIMPLIFY=F)
                  
   
## Merge btrec with post.btrecase$gt

trec.gt <- merge(post.btrecase95.1$nb, post.btrecase95.1$gt, by=c("Gene_id", "tag", "rej.level"), suffixes=c(".trec", ".gt"))
trec.gt <- add.signif(trec.gt, x1="null.rej.trec", x2="null.rej.gt", col=c("trec","btrecase") )

## plot hist for ASE assoc comparing trec-ase with trec

##' #Plot thresholds for multiple testing correction in Bayesian negative binomial (trec) and Bayesian NB-ASE (btrecase)  

trec.gtase <- trec.gt[model.gt=="trec-ase",]

table <- tab2bplot(dt=trec.gtase, var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", btrecase="#0072B2", Both= "#D55E00"))

tables <- lapply(rej, function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=rej, grob=gl)


## So far I was comparing effects using 95% CI:

## btrecase <- add.null(dt=btrecase)
## trec.refbias <- merge(trec, btrecase[model =="trec-ase",], by=c("Gene_id", "tag"), suffixes=c(".trec",".trec-ase"))
## trec.refbias <- add.signif(trec.refbias, x1= "log2_aFC_null" , x2="null.95",  col=c("trec", "trec-ase"))
 
## trec.rb.tab <- tab2bplot(trec.refbias, colors=c(None="#999999", `trec-ase`="yellow3", trec="#0072B2", Both= "#D55E00"))

## #+ fig.width= 6, fig.height=5
## btrecase.plot(dt=trec.refbias[Signif != "None",] ,
##               x2=c(paste0(c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%"), ".trec-ase"), "null.95"),
##               x1=c(paste0(c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%"), ".trec"),  "log2_aFC_null" ),
##               yl="eQTL effect Bayesian trec-ase",
##               xl="eQTL effect Bayesian trec",
##               col=c("trec-ase", "trec"),
##               title="eQTL estimates using 95% CI"
##               ) + 
##     annotation_custom(tableGrob(trec.rb.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(trec.rb.tab$color, rep("black", 4)))))), xmin=-2, xmax=-1, ymin=1, ymax=2.5)



p <- btrecase.plot(dt=trec.gtase[Signif != "None" ,] , x1=c(rep("log2_aFC_mean.trec",3),"null.rej.trec") ,
                   x2=c(rep("log2_aFC_mean.gt",3), "null.rej.gt"),
                   s=50000,
                   xl="eQTL effect Bayesian trec",
                   yl="eQTL effect Bayesian trec-ase",
                   col=c("trec", "btrecase"),
                   title=paste0("eQTL estimates based on ", 100*post.level, "% posterior out of the indicated\n rejection zone for associations with ASE"),
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
post.level=0.95

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

cols <- mapply(function(x,i,j) sapply(x,function(k) grep(paste0(k, i, "$"), names(j), value=T)),              
               i=c("", "", ".rna"),
               j=list(trec, btrecase, gt.rna[model.gt == "trec-ase",]),
               MoreArgs=list(x=c("log2_aFC_mean" ,"log2_aFC_sd", "tag")),
               SIMPLIFY=F,
               USE.NAMES=F
               )


post.btrecase95 <- mapply(function(a,b,cols,c,t) {
    dt <- rej.recode(a, b,cols, c)
    setkeyv(dt, c("Gene_id",t))  
    return(dt)},
    b=list(nb=trec, gt=btrecase, rna=gt.rna),
    cols= lapply(cols, function(i) i[1:2]),
    t=sapply(cols, function(i) i[3]),
    MoreArgs=list(a=rej, c=post.level),
    SIMPLIFY=F)
                  
   
## Merge post.btrec$nb with post.btrecase$gt

trec.gt <- merge(post.btrecase95$nb, post.btrecase95$gt[model == "trec-ase",], by=c("Gene_id", "tag", "rej.level"), suffixes=c(".trec", ".gt"))
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
                    title=paste0("eQTL estimates based on ", post.level*100, "% posterior out of the indicated\n rejection zone for associations with ASE"),
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


gt.rna2 <- merge(post.btrecase95$gt, post.btrecase95$rna, by.x=c("Gene_id", "tag", "rej.level"), by.y=c("Gene_id", "tag.gt", "rej.level"), suffixes=c(".gt", ".rna"))
                
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



##' # Compare with % posterior = 0.99 keeping same rejection zone

post.level  <- 0.99

post.btrecase <- mapply(function(a,b,cols,c,t) {
    dt <- rej.recode(a, b,cols, c)
    setkeyv(dt, c("Gene_id",t))  
    return(dt)},
    b=list(nb=trec, gt=btrecase, rna=gt.rna),
    cols= lapply(cols, function(i) i[1:2]),
    t=sapply(cols, function(i) i[3]),
    MoreArgs=list(a=rej, c=post.level),
    SIMPLIFY=F)
   
gt.rna2 <- merge(post.btrecase$gt, post.btrecase$rna, by.x=c("Gene_id", "tag", "rej.level" ), by.y=c("Gene_id", "tag.gt", "rej.level" ), suffixes=c(".gt", ".rna"))
                
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

post.per <- c(seq(0.5,0.9, 0.1), .05, .99)

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

################################################################################################################################
##' # Alternative analysis based on http://varianceexplained.org/r/bayesian_fdr_baseball/ ############

## select a rule to reject: r value excluded from posterior 95%CI
## define error: posterior probability for true effect in opposite direction of mean estimate
## Plot number of rejections and number of false positives vs r


r <- seq(0,0.8,0.2)
post.level <- 0.95

## compute posterior error probability (1 - post.out when rejection level is 0)
trec.err <- post.dt[rej.level == 0,][, PEP:= 1 -post.out][, .(Gene_id, tag, PEP)]

trec.reject <- rbindlist(lapply(r, function(i) {
    dt <- post.dt[rej.level==i,][post.out>=post.level,]
    dt <- merge(dt, trec.err, by=c("Gene_id", "tag"))
    setkey(dt, PEP)
    dt2 <- data.table(rej.level=i, `Total rejected`=nrow(dt) , `False positives`=sum(dt$PEP))
    dt2[, per.FDR:=round(`False positives`*100/`Total rejected`, 3)]
    return(dt2)
}))

trec.reject.l <- melt(trec.reject[, .(rej.level, `Total rejected`, `False positives`)], id.vars="rej.level", value.name="Associations", variable.name="type")

ggplot(trec.reject.l, aes(rej.level, Associations, color=type)) +
    geom_point() +
    geom_line() +
    xlab("Rejection level") +
    ylab("Number of rejections") +
    theme(legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.95, 0.15), legend.key=element_blank()) +
    ggtitle(paste0("Trec model using PEP = ", 1-post.level) ) +
    annotation_custom(tableGrob(trec.reject[, `False positives`:=round(`False positives`,2)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"))), xmin=.2, xmax=.8, ymin=10000, ymax=20000)



## same with post.level 0.99


post.level <- 0.99

trec.reject <- rbindlist(lapply(r, function(i) {
    dt <- post.dt[rej.level==i,][post.out>=post.level,]
    dt <- merge(dt, trec.err, by=c("Gene_id", "tag"))
    setkey(dt, PEP)
    dt2 <- data.table(rej.level=i, `Total rejected`=nrow(dt) , `False positives`=sum(dt$PEP))
    dt2[, per.FDR:=round(`False positives`*100/`Total rejected`, 3)]
    return(dt2)
}))

trec.reject.l <- melt(trec.reject[, .(rej.level, `Total rejected`, `False positives`)], id.vars="rej.level", value.name="Associations", variable.name="type")

ggplot(trec.reject.l, aes(rej.level, Associations, color=type)) +
    geom_point() +
    geom_line() +
    xlab("Rejection level") +
    ylab("Number of rejections") +
    theme(legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.95, 0.15), legend.key=element_blank()) +
    ggtitle(paste0("Trec model using PEP = ", 1-post.level)) +
    annotation_custom(tableGrob(trec.reject[, `False positives`:=round(`False positives`,2)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"))), xmin=.2, xmax=.8, ymin=3000, ymax=5000)



##########################################################################################################################
## Change definition for posterior error probability : Prob(abs(true
## effect) <= delta), delta=c(0.05,0.1) keep same decision rule

r <- seq(0.2,0.8,0.2)
post.level <- 0.95

trec.err2 <- rej.recode(a=c(0.05, 0.1), b=trec,c=post.level)[, PEP:= 1 -post.out][, .(Gene_id, tag, rej.level, PEP)]

trec.rej2 <- lapply(unique(trec.err2[['rej.level']]), function(j) rbindlist(lapply(r, function(i) {
     dt <- post.dt[rej.level==i,][post.out>=post.level,]
     dt <- merge(dt, trec.err2[rej.level==j,.(Gene_id, tag, PEP)], by=c("Gene_id", "tag"))
     setkey(dt, PEP)
     dt2 <- data.table(rej.level=i, `Total rejected`=nrow(dt) , `False positives`=sum(dt$PEP))
     dt2[, per.FDR:=round(`False positives`*100/`Total rejected`, 3)]
     return(dt2)
})))

names(trec.rej2) <- unique(trec.err2[['rej.level']])

trec.reject.2l <- lapply(trec.rej2, function(i) melt(i[, .(rej.level, `Total rejected`, `False positives`)], id.vars="rej.level", value.name="Associations", variable.name="type"))

p2 <- lapply(seq_along(trec.reject.2l), function(i) ggplot(trec.reject.2l[[i]], aes(rej.level, Associations, color=type)) +
                                   geom_point() +
                                   geom_line() +
                                   xlab("Rejection level") +
                                   ylab("Number of rejections") +
                                   theme(legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.95, 0.15), legend.key=element_blank()) +
                                   labs(title=paste("Trec model with PEP for abs(true effect) <=", names(trec.rej2)[i]),
                                        caption=paste0("Rejection rule based on ",post.level*100, "% CI" )) +
                                   annotation_custom(tableGrob(trec.rej2[[i]][, `False positives`:=round(`False positives`,2)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"))), xmin=.4, xmax=.8, ymin=1500, ymax=2500)
)
#+ fig.width= 10, fig.height=10
plot_grid(plotlist=p2, nrow=2)

################################################################################################################
##' # Compare associations with 95% posterior probability for different rejection levels with Gtex and Geuvadis


## Look at Gtex ebv 

## ebv <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz")
ebv <- fread(cmd=paste("zcat", snakemake@input[['gtex']]))

## ebv.sig <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz")

ebv.sig <- fread(cmd=paste("zcat", snakemake@input[['sigGtex']]))

ebv[, Gene_id:=gsub("\\..*","",gene_id)]
ebv <- ebv[Gene_id %in% unique(post.dt$Gene_id),]
ebv[, SNP:=gsub("^22_|_b37$", "", variant_id)][, SNP:=gsub("_", ":",SNP)]

ebv.sig <- ebv.sig[variant_id %in% ebv$variant_id,][, null:="no"]
ebv <- merge(ebv,ebv.sig[,.(gene_id, variant_id, null)], by=c("gene_id", "variant_id"), all.x=T)
ebv[is.na(null), null:="yes"]



##' # Compare Bayesian trec and DEseq to Gtex EBV cell lines

## merge post.dt with ebv 

trec.ebv <- merge(post.dt, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))

trec.ebv <- add.signif(trec.ebv, x1="null.rej", x2="null", col=c("trec", "Gtex-ebv"))


table <- tab2bplot(dt=trec.ebv[rej.level %in% r,], var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", `Gtex-ebv`="#0072B2", Both= "#D55E00"))

tables.trec <- lapply(r, function(i) {
    table[rej.level==i,]
})



gl <- lapply(tables.trec, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=r, grob=gl)


lab <- paste('r = \u00B1',r)
names(lab) <- r

##' # Bayesian trec to Gtex EBV cell lines

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=trec.ebv[Signif != "None" & rej.level %in% r,] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("slope",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trec",
                   yl="eQTL effect Gtex",
                   col=c("trec", "Gtex-ebv"),
                   title=paste0("Trec estimates based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
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


## DEseq at various FDR to Gtex EBV cell lines at 5% FDR

dseq.l <- rbindlist(lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tmp <- dseq[, null.fdr:="yes"][p.adj<=i, null.fdr:="no"]
    tmp[, Fdr:= i]
    tmp <- tmp[, .(Gene_id, tag, log2FoldChange,log2_aFC, p.adj, null.fdr, Fdr, gene.dist)]
    return(tmp)
    }))

dseq.ebv <- merge(dseq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))

dseq.ebv <-  add.signif(dseq.ebv, x1="null.fdr", x2="null", col=c("DEseq","Gtex-ebv") )

tab <- tab2bplot(dseq.ebv, var= c("Signif", "Fdr"), colors=c("None"="#999999","DEseq"="yellow3", "Both"="#D55E00", "Gtex-ebv"="#0072B2"))

tables.dseq <- lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tab[Fdr==i,]
})


gl <- lapply(tables.dseq, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(Fdr=c(0.1, 0.05, 0.01, 0.001), grob=gl)


lab <- paste('FDR(%) =', c(0.1, 0.05, 0.01, 0.001)*100)
names(lab) <- c(0.1, 0.05, 0.01, 0.001)

##' # DEseq with different FDR relative to Gtex-EBV
#+ fig.width= 12, fig.height=21
btrecase.plot(dt=dseq.ebv[Signif != "None" ,] , x1=c(rep("log2FoldChange",3),"null.fdr") ,
                   x2=c(rep("slope",3), "null"),
                   xl="eQTL effect DEseq",
                   yl="eQTL effect Gtex",
                   col=c("DEseq", "Gtex-ebv"),
                   title="DEseq2 vs Gtex-EBV (5% FDR)",
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
              ) +  facet_grid(Fdr~., labeller=labeller(Fdr=lab))+
    theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14)) +
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)



## Compare Gtex-EBV with Btrecase, all associations

btrecase.ebv <- merge(post.btrecase95.1$gt, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
btrecase.ebv <- add.signif(btrecase.ebv, x1="null.rej", x2="null", col=c("Btrecase", "Gtex-ebv"))


table <- tab2bplot(dt=btrecase.ebv, var= c("Signif", "rej.level"), colors=c(None="#999999", Btrecase="yellow3", `Gtex-ebv`="#0072B2", Both= "#D55E00"))

tables.btrecase <- lapply(unique(table$rej.level), function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables.btrecase, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=unique(table$rej.level), grob=gl)


lab <- paste('r = \u00B1',unique(table$rej.level))
names(lab) <- unique(table$rej.level)

##' # Btrecase (with ASE only) vs Gtex EBV (5% FDR)

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=btrecase.ebv[Signif != "None",] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("slope",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trecase",
                   yl="eQTL effect Gtex",
                   col=c("Btrecase", "Gtex-ebv"),
                   title=paste0("Btrecase based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
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


##' # Summary of associations relative to Gtex-EBV




## Trec: Significant and total associations by rejection level

rbindlist(lapply(tables.trec, format.tab , "Gtex-ebv"))


## DESeq:  Significant and total associations by FDR

rbindlist(lapply(tables.dseq, format.tab, "Gtex-ebv" ))


## Btrecase:  Significant and total associations by rejection level

rbindlist(lapply(tables.btrecase, format.tab , "Gtex-ebv" ))


##########################################################################

## Array Express: All the associations below false discovery rate 5% from
## https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/

## sig.geu <- fread("/mrc-bsu/scratch/ev250/EGEUV1/array_express_eqtl/EUR373.gene.cis.FDR5.all.rs137.txt")
sig.geu <- fread(snakemake@input[['geu_eur']])
sig.geu <- sig.geu[CHR_SNP==22 & CHR_GENE==22,][,Gene_id:=gsub("\\..*","", GENE_ID)]

## get ref and alt allele for each SNP_ID

snpmart <- useEnsembl(biomart="snp", dataset="hsapiens_snp", GRCh="37")

snp <- data.table(getBM(attributes=c('refsnp_id', 'chrom_start', 'allele_1', 'minor_allele'),
             filters = 'snp_filter',
             values = unique(sig.geu[['SNP_ID']]),
             mart=snpmart))

snp[, id:=paste(chrom_start, allele_1, minor_allele, sep=":")]
sig.geu <- merge(sig.geu, snp, by.x=c("SNP_ID", "SNPpos"), by.y=c("refsnp_id","chrom_start"))
sig.geu[, "null":="no"]


## GEUVADIS analysis from Chris
## select  chrom 22 and format compatible with gt
## sig.qtl <- fread('/home/ev250/rds/rds-cew54-wallace-share/Data/GEUVADIS/Geuvadis_sig_eqtl')
sig.qtl <- fread(snakemake@input[['geu_chris']])
sig.qtl <- sig.qtl[CHROM==22,][,Gene_id:=gsub("\\..*","", Gene)]
sig.qtl[, id := paste(POS,REF,ALT, sep=":")][, "null":="no"]
sig.qtl[, aFC:=2*Beta]


##' # Compare GEUVADIS EUR.gene.cis.FDR5 with GEUVADIS Chris for association tested within chromosome 22:

comp <- merge(sig.geu, sig.qtl, by.x=c("GENE_ID", "SNP_ID", "id"), by.y=c("Gene", "SNP", "id"), suffixes=c(".cis.fdr5", ".Chris"), all=T)

comp[,.N, .(null.cis.fdr5, null.Chris)]



##' # Compare Bayesian trec assuming missing associations in whole GEUVADIS signif dataset are not significant

rej <- c(0,r)

trec.sig <- merge(post.dt,sig.qtl,by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
trec.sig[is.na(null), null:="yes"]
trec.sig <- add.signif(trec.sig, x1="null.rej", x2="null", col=c("trec", "GEU-CHRIS"))


table <- tab2bplot(dt=trec.sig[rej.level %in% rej ,], var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", `GEU-CHRIS`="#0072B2", Both= "#D55E00"))

tables.trec.sig <- lapply(rej, function(i) {
    table[rej.level==i,]
})


gl <- lapply(tables.trec.sig, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=sort(rej), grob=gl)

lab <- paste('r = \u00B1',rej)
names(lab) <- rej


#+ fig.width= 12, fig.height=21
btrecase.plot(dt=trec.sig[Signif != "None" & rej.level %in% rej,] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("aFC",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trec",
                   yl="eQTL effect GEUVADIS",
                   col=c("trec", "GEU-CHRIS"),
                   title=paste0("Trec estimates based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
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

## DEseq

dseq.sig <- merge(dseq.l, sig.qtl, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
dseq.sig[is.na(null), null:="yes"]
dseq.sig <-  add.signif(dseq.sig, x1="null.fdr", x2="null", col=c("DEseq","GEU-CHRIS") )

tab <- tab2bplot(dseq.sig, var= c("Signif", "Fdr"), colors=c(None="#999999",DEseq="yellow3", Both="#D55E00", `GEU-CHRIS`="#0072B2"))

tables.dseq.sig <- lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tab[Fdr==i,]
})


gl <- lapply(tables.dseq.sig, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(Fdr=c(0.1, 0.05, 0.01, 0.001), grob=gl)


lab <- paste('FDR(%) =', c(0.1, 0.05, 0.01, 0.001)*100)
names(lab) <- c(0.1, 0.05, 0.01, 0.001)

##' # DEseq with different FDR relative to GEU-CHRIS
#+ fig.width= 12, fig.height=21
btrecase.plot(dt=dseq.sig[Signif != "None" ,] , x1=c(rep("log2FoldChange",3),"null.fdr") ,
              x2=c(rep("Beta",3), "null"),
              xl="eQTL effect DEseq",
              s=nrow(dseq.sig[Signif != "None" ,]),
              yl="eQTL effect GEU_CHRIS",
              col=c("DEseq", "GEU_CHRIS"),
              title="DEseq2 vs GEU-CHRIS (5% FDR)",
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) +  facet_grid(Fdr~., labeller=labeller(Fdr=lab))+
    theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14)) +
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)


## Btrecase, all associations

btrecase.sig <- merge(post.btrecase95.1$gt, sig.qtl, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
btrecase.sig[is.na(null), null:="yes"]
btrecase.sig <- add.signif(btrecase.sig, x1="null.rej", x2="null", col=c("Btrecase", "GEU-CHRIS"))


table <- tab2bplot(dt=btrecase.sig, var= c("Signif", "rej.level"), colors=c(None="#999999", Btrecase="yellow3", `GEU-CHRIS`="#0072B2", Both= "#D55E00"))

tables.btrecase.sig <- lapply(unique(table$rej.level), function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables.btrecase.sig, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=unique(table$rej.level), grob=gl)

lab <- paste('r = \u00B1',unique(table$rej.level))
names(lab) <- unique(table$rej.level)

##' # Btrecase (with ASE only) vs GEU_CHRIS (5% FDR)

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=btrecase.sig[Signif != "None",] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("Beta",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trecase",
                   yl="eQTL effect GEU-CRHIS",
                   col=c("Btrecase", "GEU-CHRIS"),
                   title=paste0("Btrecase based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
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


##' # Summary of associations relative to GEU-CHRIS

## Trec: Significant and total associations by rejection level

rbindlist(lapply(tables.trec.sig, format.tab , "GEU-CHRIS"))


## DESeq:

rbindlist(lapply(tables.dseq.sig, format.tab, "GEU-CHRIS"))


## Btrecase:

rbindlist(lapply(tables.btrecase.sig, format.tab, "GEU-CHRIS"))


############################################################################################################
##' # Compare Bayesian trec, DEseq and Btrecase with 373 EUR GEUVADIS samples significant associations at 5% FDR  assuming missing associations in GEUVADIS signif dataset are not significant

rej <- c(0,r)

trec.geu <- merge(post.dt, sig.geu, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
trec.geu[is.na(null), null:="yes"]
trec.geu <- add.signif(trec.geu, x1="null.rej", x2="null", col=c("trec", "GEU-EUR"))


table <- tab2bplot(dt=trec.geu[rej.level %in% rej ,], var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", `GEU-EUR`="#0072B2", Both= "#D55E00"))

tables.trec.geu <- lapply(rej, function(i) {
    table[rej.level==i,]
})


gl <- lapply(tables.trec.geu, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=sort(rej), grob=gl)

lab <- paste('r = \u00B1',rej)
names(lab) <- rej


#+ fig.width= 12, fig.height=21
btrecase.plot(dt=trec.geu[Signif != "None" & rej.level %in% rej,] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("rvalue",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trec",
                   yl="eQTL effect GEUVADIS",
                   col=c("trec", "GEU-EUR"),
                   title=paste0("Trec estimates based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
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

## DEseq

dseq.geu <- merge(dseq.l, sig.geu, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
dseq.geu[is.na(null), null:="yes"]
dseq.geu <-  add.signif(dseq.geu, x1="null.fdr", x2="null", col=c("DEseq","GEU-EUR") )

tab <- tab2bplot(dseq.geu, var= c("Signif", "Fdr"), colors=c(None="#999999",DEseq="yellow3", Both="#D55E00", `GEU-EUR`="#0072B2"))

tables.dseq.geu <- lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tab[Fdr==i,]
})


gl <- lapply(tables.dseq.geu, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(Fdr=c(0.1, 0.05, 0.01, 0.001), grob=gl)


lab <- paste('FDR(%) =', c(0.1, 0.05, 0.01, 0.001)*100)
names(lab) <- c(0.1, 0.05, 0.01, 0.001)

##' # DEseq with different FDR relative to GEU-EUR
#+ fig.width= 12, fig.height=21
btrecase.plot(dt=dseq.geu[Signif != "None" ,] , x1=c(rep("log2FoldChange",3),"null.fdr") ,
              x2=c(rep("rvalue",3), "null"),
              xl="eQTL effect DEseq",
              s=nrow(dseq.geu[Signif != "None" ,]),
              yl="eQTL effect GEU-EUR",
              col=c("DEseq", "GEU-EUR"),
              title="DEseq2 vs GEU-EUR (5% FDR)",
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) +  facet_grid(Fdr~., labeller=labeller(Fdr=lab))+
    theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14)) +
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)


## Btrecase

btrecase.geu <- merge(post.btrecase95.1$gt, sig.geu, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id"), all.x=T)
btrecase.geu[is.na(null), null:="yes"]
btrecase.geu <- add.signif(btrecase.geu, x1="null.rej", x2="null", col=c("Btrecase", "GEU-EUR"))


table <- tab2bplot(dt=btrecase.geu, var= c("Signif", "rej.level"), colors=c(None="#999999", Btrecase="yellow3", `GEU-EUR`="#0072B2", Both= "#D55E00"))

tables.btrecase.geu <- lapply(unique(table$rej.level), function(i) {
    table[rej.level==i,]
})

gl <- lapply(tables.btrecase.geu, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(rej.level=unique(table$rej.level), grob=gl)

lab <- paste('r = \u00B1',unique(table$rej.level))
names(lab) <- unique(table$rej.level)

##' # Btrecase (with ASE only) vs GEU_EUR (5% FDR)

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=btrecase.geu[Signif != "None",] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("rvalue",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trecase",
                   yl="eQTL effect GEU-EUR",
                   col=c("Btrecase", "GEU-EUR"),
                   title=paste0("Btrecase based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
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


##' # Summary of associations relative to GEU-EUR

## Trec:

rbindlist(lapply(tables.trec.geu, format.tab , "GEU-EUR"))


## DESeq:

rbindlist(lapply(tables.dseq.geu, format.tab, "GEU-EUR"))


## Btrecase:

rbindlist(lapply(tables.btrecase.geu, format.tab, "GEU-EUR"))


#########################################################################################################

## Get summary of associations using as gold standard significant in gtex or geu, but tested in both?

gtex.geu.eur <- merge(ebv, sig.geu, by.x=c("Gene_id", "SNP"), by.y=c("Gene_id", "id"), suffixes=c(".ebv", ".geu.eur"), all=T)

gtex.geu <- merge(gtex.geu.eur , sig.qtl, by.x=c("Gene_id", "SNP"), by.y=c("Gene_id", "id"), all=T)

## Code "null" to "no" if a signif association in any dataset

gtex.geu[ (null.ebv =="no" | null.geu.eur == "no"), null := "no"]

gtex.geu[is.na(null), null:="yes"]


## Merge gtex.geu with trec, btrecase btrecase rna and deseq

trec.gtex.geu <- merge(post.dt, gtex.geu, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))

trec.gtex.geu<- add.signif(trec.gtex.geu, x1="null.rej", x2="null", col=c("trec", "Gtex_GEU"))

table <- tab2bplot(dt=trec.gtex.geu[rej.level %in% rej ,], var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", `Gtex_GEU`="#0072B2", Both= "#D55E00"))

tables.trec.gtex.geu <- lapply(rej, function(i) {
    table[rej.level==i,]
})

## Btrecase GT

btrecase.gtex.geu <- merge(post.btrecase95.1$gt, gtex.geu, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
btrecase.gtex.geu<- add.signif(btrecase.gtex.geu, x1="null.rej", x2="null", col=c("Btrecase", "Gtex_GEU"))

table <- tab2bplot(dt=btrecase.gtex.geu, var= c("Signif", "rej.level"), colors=c(None="#999999", Btrecase="yellow3", `Gtex_GEU`="#0072B2", Both= "#D55E00"))

tables.btrecase.gtex.geu <- lapply(unique(table$rej.level), function(i) {
    table[rej.level==i,]
})

## DEseq

dseq.gtex.geu <- merge(dseq.l, gtex.geu, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
dseq.gtex.geu <-  add.signif(dseq.gtex.geu, x1="null.fdr", x2="null", col=c("DEseq","Gtex_GEU") )

tab <- tab2bplot(dseq.gtex.geu, var= c("Signif", "Fdr"), colors=c(None="#999999",DEseq="yellow3", Both="#D55E00", `Gtex_GEU`="#0072B2"))

tables.dseq.gtex.geu <- lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tab[Fdr==i,]
})


## Btrecase RNA

rna.gtex.geu <- merge(post.btrecase95$rna, gtex.geu, by.x=c("Gene_id", "tag.gt"), by.y=c("Gene_id", "SNP"))
rna.gtex.geu<- add.signif(rna.gtex.geu, x1="null.rej", x2="null", col=c("Hidden-GT", "Gtex_GEU"))

table <- tab2bplot(dt=rna.gtex.geu, var= c("Signif", "rej.level"), colors=c(None="#999999", `Hidden-GT`="yellow3", `Gtex_GEU`="#0072B2", Both= "#D55E00"))

tables.rna.gtex.geu <- lapply(unique(table$rej.level), function(i) {
    table[rej.level==i,]
})
##' # Summary of associations relative to Gtex_GEU

## Trec:

rbindlist(lapply(tables.trec.gtex.geu, format.tab , "Gtex_GEU"))


## DESeq:

rbindlist(lapply(tables.dseq.gtex.geu, format.tab, "Gtex_GEU"))


## Btrecase:

rbindlist(lapply(tables.btrecase.gtex.geu, format.tab, "Gtex_GEU"))

## RNA:

rbindlist(lapply(tables.rna.gtex.geu, format.tab, "Gtex_GEU"))


##########################################################################################################################
##' # Look at TDR by maf,  distance to TSS, and gene expression


## add EAF to trec and dseq

post.dt <- merge(post.dt, unique(post.btrecase95.1$gt[,.(Gene_id,tag,tag.EAF)]), by=c("Gene_id", "tag"))
trec.gtex.geu <- merge(trec.gtex.geu, unique(post.btrecase95.1$gt[,.(Gene_id,tag,tag.EAF)]), by=c("Gene_id", "tag"))
dseq.gtex.geu <- merge(dseq.gtex.geu, unique(post.btrecase95.1$gt[,.(Gene_id,tag,tag.EAF)]), by=c("Gene_id", "tag"))

## create categorical variables with EAF and gene distance: EAF low if <=0.25 | >=0.75, high otherwise
## gene distance: far >= 200KB, medium 200-100, close <100


post.dt <- rec.cat(dt=post.dt, col="tag.EAF", newcol="EAF", cats=c(0, 0.25,0.75,1), labs=c("low", "high", "low"))
trec.gtex.geu <- rec.cat(trec.gtex.geu, col="tag.EAF", newcol="EAF", cats=c(0, 0.25,0.75,1), labs=c("low", "high", "low"))
dseq.gtex.geu <- rec.cat(dt=dseq.gtex.geu, col="tag.EAF", newcol="EAF", cats=c(0, 0.25,0.75,1), labs=c("low", "high", "low"))
btrecase.gtex.geu <- rec.cat(btrecase.gtex.geu, col="tag.EAF", newcol="EAF", cats=c(0, 0.25,0.75,1), labs=c("low", "high", "low"))
rna.gtex.geu <- rec.cat(rna.gtex.geu, col="tag.EAF.gt", newcol="EAF.gt", cats=c(0, 0.25,0.75,1), labs=c("low", "high", "low"))

trec.gtex.geu <- rec.cat(dt=trec.gtex.geu, col="gene.dist", fun="abs", newcol="Distance2gene", cats=c(-1,10^3,5*10^5), labs=c("close",  "far"))
dseq.gtex.geu <- rec.cat(dt=dseq.gtex.geu, col="gene.dist", fun="abs", newcol="Distance2gene", cats=c(-1,10^3,5*10^5), labs=c("close", "far"))

btrecase.gtex.geu <- rec.cat(dt=btrecase.gtex.geu, col="gene.dist", fun="abs", newcol="Distance2gene", cats=c(-1,10^3,5*10^5), labs=c("close",  "far"))

rna.gtex.geu <- rec.cat(rna.gtex.geu, col="gene.dist.gt", fun="abs", newcol="Distance2gene", cats=c(-1,10^3,5*10^5), labs=c("close", "far"))

##' # Tables by EAF, low <=0.25, high >0.25

## Trec:
TDR.var(dt=trec.gtex.geu, var= c("Signif", "rej.level", "EAF"), cols= c( "rej.level", "EAF"),  gold= "Gtex_GEU")

## Dseq:

TDR.var(dt=dseq.gtex.geu, var= c("Signif", "Fdr", "EAF"), cols= c( "Fdr", "EAF"),  gold= "Gtex_GEU")

## Btrecase GT:

TDR.var(dt=btrecase.gtex.geu, var= c("Signif", "rej.level", "EAF"), cols= c( "rej.level", "EAF"),  gold= "Gtex_GEU")

  
## Btrecase Hidden-GT:

TDR.var(dt=rna.gtex.geu, var= c("Signif", "rej.level", "EAF.gt"), cols= c( "rej.level", "EAF.gt"),  gold= "Gtex_GEU")


##' # Tables by Gene distance, close <= 100,000, far >100,000 up to 5*10^5


##' Trec:
TDR.var(dt=trec.gtex.geu, var= c("Signif", "rej.level", "Distance2gene"), cols= c( "rej.level", "Distance2gene"),  gold= "Gtex_GEU")

##' Dseq:

TDR.var(dt=dseq.gtex.geu, var= c("Signif", "Fdr", "Distance2gene"), cols= c( "Fdr", "Distance2gene"),  gold= "Gtex_GEU")

##' Btrecase GT:

TDR.var(dt=btrecase.gtex.geu, var= c("Signif", "rej.level", "Distance2gene"), cols= c( "rej.level", "Distance2gene"),  gold= "Gtex_GEU")

  
##' Btrecase Hidden-GT:

TDR.var(dt=rna.gtex.geu, var= c("Signif", "rej.level", "Distance2gene"), cols= c( "rej.level", "Distance2gene"),  gold= "Gtex_GEU")

##' # TDR by Gene expression

counts <- fread(counts.f)

## get median expression per gene

med.g <- data.table(Gene_id = counts$gene_id, median.g = apply(counts[,-1], 1, median))
med.g[, Gene_expression:='low'][median.g >= median(median.g), Gene_expression:='high']

## add median expression to datasets:

trec.gtex.geu <- merge(trec.gtex.geu, med.g[,.(Gene_id,Gene_expression)], by=c("Gene_id"))
dseq.gtex.geu <- merge(dseq.gtex.geu, med.g[,.(Gene_id,Gene_expression)], by=c("Gene_id"))
btrecase.gtex.geu <- merge(btrecase.gtex.geu, med.g[,.(Gene_id,Gene_expression)], by=c("Gene_id"))
rna.gtex.geu  <- merge(rna.gtex.geu,  med.g[,.(Gene_id,Gene_expression)], by=c("Gene_id"))

##' # Tables by Gene expression, low < median, high >= median


##' Trec:
TDR.var(dt=trec.gtex.geu, var= c("Signif", "rej.level", "Gene_expression"), cols= c( "rej.level", "Gene_expression"),  gold= "Gtex_GEU")[order(Gene_expression,rej.level),]

##' Dseq:

TDR.var(dt=dseq.gtex.geu, var= c("Signif", "Fdr", "Gene_expression"), cols= c( "Fdr", "Gene_expression"),  gold= "Gtex_GEU")[order(Gene_expression,-Fdr),]

##' Btrecase GT:

TDR.var(dt=btrecase.gtex.geu, var= c("Signif", "rej.level", "Gene_expression"), cols= c( "rej.level", "Gene_expression"),  gold= "Gtex_GEU")[order(Gene_expression,rej.level),]

  
##' Btrecase Hidden-GT:

TDR.var(dt=rna.gtex.geu, var= c("Signif", "rej.level", "Gene_expression"), cols= c( "rej.level", "Gene_expression"),  gold= "Gtex_GEU")[order(Gene_expression,rej.level),]

##' Select some cut-offs and plot the data, high maf, close to gene and high expression

## recode an effect size for geu.gtex. when is.na(Gtex) use GEU-EUR, otherwise GEU-CHRIS

trec.gtex.geu[,merged.estimate:=slope][is.na(slope), merged.estimate:=Beta][is.na(slope) & is.na(Beta), merged.estimate:=rvalue]

dseq.gtex.geu[,merged.estimate:=slope][is.na(slope), merged.estimate:=Beta][is.na(slope) & is.na(Beta), merged.estimate:=rvalue]

btrecase.gtex.geu[,merged.estimate:=slope][is.na(slope), merged.estimate:=Beta][is.na(slope) & is.na(Beta), merged.estimate:=rvalue]

rna.gtex.geu[,merged.estimate:=slope][is.na(slope), merged.estimate:=Beta][is.na(slope) & is.na(Beta), merged.estimate:=rvalue]



##' Trec:
TDR.var(dt=trec.gtex.geu[Signif != "None" & EAF=="high" & Distance2gene=="close" & Gene_expression =="high",], var= c("Signif", "rej.level"), cols= c( "rej.level"),  gold= "Gtex_GEU")

##' Dseq:

TDR.var(dt=dseq.gtex.geu[Signif != "None" & EAF=="high" & Distance2gene=="close" & Gene_expression =="high",], var= c("Signif", "Fdr"), cols= c( "Fdr"),  gold= "Gtex_GEU")

##' Btrecase GT:

TDR.var(dt=btrecase.gtex.geu[Signif != "None" & EAF=="high" & Distance2gene=="close" & Gene_expression =="high",], var= c("Signif", "rej.level"), cols= c( "rej.level"),  gold= "Gtex_GEU")

  
##' Btrecase Hidden-GT:

TDR.var(dt=rna.gtex.geu[Signif != "None" & EAF.gt=="high" & Distance2gene=="close" & Gene_expression =="high",], var= c("Signif", "rej.level"), cols= c( "rej.level"),  gold= "Gtex_GEU")

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=trec.gtex.geu[Signif != "None" & EAF=="high" & Distance2gene=="close" & Gene_expression =="high",] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("merged.estimate",3), "null"),
                   s=50000,
                   xl="eQTL effect Trec",
                   yl="eQTL effect Gtex-GEU",
                   col=c("Trec", "Gtex-GEU"),
                   title=paste0("Trec based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
              ) + facet_grid(rej.level~.)#, labeller=labeller(rej.level=lab))


#+ fig.width= 12, fig.height=21
btrecase.plot(dt=dseq.gtex.geu[Signif != "None" & EAF=="high" & Distance2gene=="close" & Gene_expression =="high",] , x1=c(rep("log2FoldChange",3),"null.fdr") ,
                   x2=c(rep("merged.estimate",3), "null"),
                   s=50000,
                   xl="eQTL effect DEseq",
                   yl="eQTL effect Gtex-GEU",
                   col=c("DEseq", "Gtex-GEU"),
                   title=paste0("DEseq by FDR"),
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
              ) + facet_grid(Fdr~.)#, labeller=labeller(rej.level=lab))

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=btrecase.gtex.geu[Signif != "None" & EAF=="high" & Distance2gene=="close" & Gene_expression =="high",] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                   x2=c(rep("merged.estimate",3), "null"),
                   s=50000,
                   xl="eQTL effect Bayesian trecase",
                   yl="eQTL effect Gtex-GEU",
                   col=c("Btrecase", "Gtex-GEU"),
                   title=paste0("Btrecase based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
                   ) + facet_grid(rej.level~.)#, labeller=labeller(rej.level=lab))

#+ fig.width= 12, fig.height=21
btrecase.plot(dt=rna.gtex.geu[Signif != "None" & EAF.gt=="high" & Distance2gene=="close" & Gene_expression =="high",
                              -which(names(rna.gtex.geu) =="op.dir"), with=F] ,
              x1=c(rep("log2_aFC_mean.rna",3),"null.rej") ,
              x2=c(rep("merged.estimate",3), "null"),
              s=50000,
              xl="eQTL effect Btrecase hidden-GT",
              yl="eQTL effect Gtex-GEU",
              col=c("Hidden-GT", "Gtex-GEU"),
              title=paste0("Hidden-GT based on ", 100*post.level, "% posterior out of the indicated\n rejection zone"),
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) + facet_grid(rej.level~.)#, labeller=labeller(rej.level=lab))


##' # Test prior fit to Gtex-EBV data

## Randomly select 10^6 rows. Select rows every 150 aiming to get independent SNPs. File ordered by gene and SNP

## ebv <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz | awk 'NR == 1 || NR % 150 == 0' | head -n 1000000")
ebv <- fread(cmd=paste("zcat", snakemake@input[['gtex']], "| awk 'NR == 1 || NR % 150 == 0' | head -n 1000000" ))

## Fit a mixture of Gaussians

mixmdl = normalmixEM(ebv$slope,lambda=c(0.05,0.95),mu=c(0,0), mean=c(0,0))
x <- rmvnormmix(n=nrow(ebv), lambda=c(0.3, 0.7), mu=c(0,0), sigma=c(0.14, 0.46))

dt <- data.table(ebv=ebv$slope)


## save for Chris: Could you send me a .csv file with one column of the gtex effect estimate and one column the standard error of the estimate, for a large sample of independent snps?

ebv <- ebv[, .(slope, slope_se)]

##write.table(ebv, file="/home/ev250/newshare/elena/trecase/objects/gtex_ebv_indpendent_snps.csv", sep=",", row.names=F)
