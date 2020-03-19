#' ---
#' title: Set threshold for making significant calls with Bayesian model using Trec with mixed Gaussian prior
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
## library(biomaRt)
library(mixtools)


source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")
source('/home/ev250/Cincinatti/Functions/various.R')

##' # Compare associations using  Bayesian trec with old and new priors

## trec <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.stan.summary.txt")

trec <- comb.files(path=snakemake@params[['trec']], pattern="trec.stan.summary.txt")

## Add EAF

## le.file <-  "/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz"
le.file <- snakemake@input[['lefile']]

EAF <-  snp.eaf(le.file, unique(trec$tag))
setnames(EAF, "eaf", "tag.EAF")

trec <- merge(trec, EAF, by.x="tag", by.y="snp")


## new priors: mix1, mix2, mix3

## trec.m <- lapply(paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrec/MixedPrior", c("","2", "3")), function(i) comb.files(path=i,pattern=".stan.summary.txt"))

trec.m <- lapply(c(snakemake@params[['trec_mix']], snakemake@params[['trec_mix2']], snakemake@params[['trec_mix3']]) , function(i) comb.files(i, pattern=".stan.summary.txt"))

trec.m  <- lapply(trec.m, add.null)

names(trec.m) <- paste0("mix", 1:3)


## merge

trec.comp <- lapply(trec.m, function(i) merge(trec, i, by=c("Gene_id", "tag"), suffixes=c(".norm", ".mix")))

trec.comp <- lapply(trec.comp, function(i) add.signif(i, "log2_aFC_null", "null.95", c("Normal", "Mix")))


cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%")

##' Normal prior: N(0,0.3)

##' Mix-1: 70% N(0,0.02) + 30% N(0,0.025)

##' Mix-2: 90% N(0,0.02) + 10% N(0,0.025)

##' Mix-3: 11% N(0,0.036) + 46% N(0,0.01) + 43% N(0,0.07)
l1 <- lapply(seq_along(trec.comp), function (i) btrecase.plot(dt=trec.comp[[i]][Signif != "None" ,] ,
                                                        x1=c(paste0(cols, ".norm") ,"log2_aFC_null" ), 
                                                        x2=c(paste0(cols,  ".mix"), "null.95" ),
                                                        #s=nrow(trec.comp[[i]][Signif != "None" ,]),
                                                        xl="eQTL effect normal prior",
                                                        yl="eQTL effect mix prior",
                                                        col=c("Normal", "Mix"),
                                                        title=paste0("Trec with normal vs\n Gaussian mix-",i ," prior"),
                                                        title.size=16,
                                                        axis.title = 14,
                                                        axis.text=12,
                                                        legend.title=14,
                                                        legend.text=12,
                                                        legend.symbol=5,
                                                        point.size=3
                                                        ))
#+ fig.width= 8, fig.height=18
plot_grid(plotlist=l1, nrow=3)
lapply(trec.comp, function(i) i[,.N,Signif])

###############################################################################################################################################
##' # Compare trec model run with different priors to DEseq and Gtex-ebv

## append normal prior to trec.m

trec.all <- c(list(normal=trec), trec.m)

## add gene distance to trec

## gene.coord <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt")

gene.coord <- fread(snakemake@input[['geneStEnd']]) 

## select genes in chrom 22

gt22 <- gene.coord[chrom==22,]

## add tag distance to gene (closest to start or end)

trec.all <- lapply(trec.all, function(i) gene.d(i, gt22[, .(gene_id, start,end,chrom)]))


## Bayesian trec.all: using normal approximation calculate proportion of
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

post.dt <- lapply(trec.all, function(i) {
    dt <- rej.recode(a=rej, b=i,c=post.level)
    setkey(dt, Gene_id,tag)
    })


## freq: add fdr at 0.1, 0.01, 0.001
## post.per: 0.5-0.95, 0.99

fdr <- c(0.001, 0.01, 0.1)


####################################################################################################################################################################
## Compare with DEseq2 NB model (run by Wei-Yu)

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

## Add gene distance

dseq <- gene.d(dseq, gt22[, .(gene_id, start,end,chrom)], snp="tag")

## Add EAF.tag

dseq <- merge(dseq, EAF, by.x="tag", by.y="snp")


###############  Compare trec.all  with Dseq2 NB model ###############

## merge post.l and dseq

btrec.dsq <- lapply(post.dt, function(i) {
    dt <- merge(i,dseq, by=c("Gene_id", "tag"))
    setkey(dt, p.adj)
    dt <- add.signif(dt, x1="null.rej", x2="null.fdr5", col=c("trec","dseq") )
    })


##' #Plot thresholds for calling significant associations in Bayesian negative binomial  and compare with NB-Dseq2  model at 5%FDR
lab <- paste('r = \u00B1',rej)
names(lab) <- rej

#+ fig.width= 8, fig.height=21
l2 <- mapply(function(k,z){
    table <- k[,.N,.(Signif, rej.level)]
    setnames(table, "N", "SNPs")
    table[Signif=="None", color:="#999999"][Signif=="trec", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="dseq", color:="#0072B2"]
    table[, Signif:=factor(Signif, levels=c("None","trec", "dseq","Both"))]
    setkey(table,Signif)
    tables <- lapply(rej, function(i) {
        table[rej.level==i,]
    })
    gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))
    dt.tab <- data.table(rej.level=rej, grob=gl)
    p <- btrecase.plot(dt=k[Signif != "None",] , x1=c(rep(cols[1],3),"null.rej") ,
              x2=c(rep("log2_aFC",3), "null.fdr5"),
              #s=nrow(k[Signif != "None" ,]),
              xl="eQTL effect Trec",
              yl="eQTL effect Dseq2",
              col=c("trec", "dseq"),
              title=paste0("DEseq 5% FDR and ", post.level*100 ,"% PIP  by rejection zone\n Trec with ",z, " prior" ),
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
    
    print(p)
    return(tables)
    }, k=btrec.dsq, z=names(btrec.dsq), SIMPLIFY=F)


##' # Extend analysis using a range of FDR and a range of %of posterior and calculate Jaccard index (intersection over union)

for (i in fdr){
    dseq[,eval(paste0("null.fdr",i*100)):= "yes"][p.adj<=i, paste0("null.fdr",i*100):="no"]
}

btrec.dsq2 <- lapply(post.dt, function(i) {
    dt <- merge(i,dseq, by=c("Gene_id", "tag"))
    setkey(dt, p.adj)
    })

## Prepare tables for Jaccard index (Both/Total excluding None)
post.per <- c(seq(0.5,0.9, 0.1), 0.95, 0.99)

JI <- lapply(btrec.dsq2, function(k) {
    rbindlist(lapply(post.per, function(i) {
    rbindlist(lapply(c(fdr, 0.05), function(j){
        tmp <- k[,.N , .(get(paste0("null.fdr", j*100)), rej.level, post.out >=i)]
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
    })

temp <- mapply(function(i,j){
    ##+ fig.width= 25, fig.height=10   
    p <- ggplot(i, aes(fdr*100, JI, color=as.factor(post.per))) +
        geom_line() + geom_point() +
        facet_grid(.~rej.level) +
        ggtitle(paste0("Jaccard index by FDR and PIP(%)\nTrec with ", j)) + xlab("FDR (%)") +
        theme(axis.text.x = element_text(angle = 45)) +
        scale_x_log10(breaks =c(0.1,  1,5,10)) + labs(color = "Post (%)")
    print(p)
}, i=JI, j=names(JI) )

  
maxJI <- lapply(JI, function(i) i[JI==max(JI),])


btrec.dsq2 <- lapply(btrec.dsq2, function(i) add.signif(i, x1="null.rej", x2="null.fdr10", col=c("Btrec","Dseq2") ))

sig.cols <- c("None"="#999999","Btrec"="yellow3", "Both"="#D55E00", "Dseq2"="#0072B2")

##+ fig.width= 4, fig.height=5
temp <- mapply(function(a,b,d) {
    
    tab <- tab2bplot(a[rej.level==d$rej.level,], colors=sig.cols)
   
   
    p <- btrecase.plot(dt=a[rej.level==d$rej.level& Signif !="None" & !(post.out<=d$post.per & eval(paste0("null.fdr",d$fdr)) =="yes")  ,] ,
                  x1=c(rep(cols[1],3),"null.rej") ,
                  x2=c(rep("log2_aFC",3), "null.fdr10"),
                  xl="eQTL effect Btrec",
                  yl="eQTL effect Dseq2",
                  col=c("Btrec", "Dseq2"),
                  title=paste0("DEseq ",100*d$fdr, "% FDR and\n Trec ", d$post.per*100,"% PIP for\n rejection zone = ", d$rej.level,"\n Trec prior is ", b ),
                  ) +
        annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=10, ymax=15) +
        annotation_custom(grobTree(textGrob(paste("Jaccard index =", round(d$JI,2)), x=0.2, y=0.7, gp=gpar(fontsize=10))))
    print(p)
    }, a=btrec.dsq2, b=names(btrec.dsq2), d=maxJI, SIMPLIFY=F)
    

#########################################################################################################################################
##' # Compare associations with 95% posterior probability for different rejection levels with Gtex EBV

## Look at Gtex ebv for chromosome 22, second  

## ebv <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz | awk -F'\t' 'NR == 1 || $2~/^22_/' ", header=T, sep="\t")
ebv <- fread(cmd=paste("zcat", snakemake@input[['gtex']], "| awk -F'\t' 'NR==1 || $2~/^22_/' "))

## ebv.sig <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz")

ebv.sig <- fread(cmd=paste("zcat", snakemake@input[['sigGtex']]))

ebv[, Gene_id:=gsub("\\..*","",gene_id)]
ebv <- ebv[Gene_id %in% unique(trec$Gene_id),]
ebv[, SNP:=gsub("^22_|_b37$", "", variant_id)][, SNP:=gsub("_", ":",SNP)]

ebv.sig <- ebv.sig[variant_id %in% ebv$variant_id,][, null:="no"]
ebv <- merge(ebv,ebv.sig[,.(gene_id, variant_id, null)], by=c("gene_id", "variant_id"), all.x=T)
ebv[is.na(null), null:="yes"]



##' # Compare Bayesian trec and DEseq to Gtex EBV cell lines

## merge post.dt with ebv 

trec.ebv <- lapply(post.dt, function(i) {
    dt <- merge(i, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
    dt <- add.signif(dt, x1="null.rej", x2="null", col=c("trec", "Gtex-ebv"))
    })
##' Bayesian trec to Gtex EBV cell lines
#+ fig.width= 8, fig.height=21
tables.trec <- mapply(function(a, b) {
    table <- tab2bplot(dt=a[rej.level %in% rej,], var= c("Signif", "rej.level"), colors=c(None="#999999", trec="yellow3", `Gtex-ebv`="#0072B2", Both= "#D55E00"))

    tab.trec <- lapply(rej, function(i) {
        table[rej.level==i,]
    })

    gl <- lapply(tab.trec, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

    dt.tab <- data.table(rej.level=rej, grob=gl)
    



    p <- btrecase.plot(dt=a[Signif != "None" & rej.level %in% rej,] , x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                  x2=c(rep("slope",3), "null"),
                  #s=50000,
                  xl="eQTL effect Bayesian trec",
                  yl="eQTL effect Gtex",
                  col=c("trec", "Gtex-ebv"),
                  title=paste0("Trec using ", 100*post.level, "% PIP rejection zone\n with prior ", b),
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
    print(p)
    return(tab.trec)
    } ,
    a=trec.ebv,
    b=names(trec.ebv),
    SIMPLIFY=F)
    
## DEseq at various FDR to Gtex EBV cell lines at 5% FDR

dseq.l <- rbindlist(lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tmp <- dseq[, null.fdr:="yes"][p.adj<=i, null.fdr:="no"]
    tmp[, Fdr:= i]
    tmp <- tmp[, .(Gene_id, tag, tag.EAF, log2FoldChange,log2_aFC, p.adj, null.fdr, Fdr, gene.dist)]
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


##################################################################################################################################
##' Summary of associations relative to Gtex-EBV


## Trec: Significant and total associations by rejection level

lapply(tables.trec, function(i) rbindlist(lapply(i, format.tab , "Gtex-ebv")))


## DESeq:  Significant and total associations by FDR

rbindlist(lapply(tables.dseq, format.tab, "Gtex-ebv" ))

####################################################################################################################################
##' Look at TDR by maf,  distance to TSS, and gene expression



## create categorical variables with EAF and gene distance: EAF low if <=maf | >=maf, high otherwise
## gene distance: far >= d, close <d

maf <- c(0.1, 0.25)
d <- 100000

trec.ebv <- lapply(trec.ebv, function(i) rec.cat(i, col="tag.EAF", newcol="EAF",
                                                 cats=sort(c(0, maf, 1-maf,1)), labs=c("low", "med","high","med", "low")))
dseq.ebv <- rec.cat(dt=dseq.ebv, col="tag.EAF", newcol="EAF",
                    cats=sort(c(0, maf, 1-maf,1)), labs=c("low", "med","high","med", "low"))

trec.ebv <-lapply(trec.ebv, function(i)  rec.cat(dt=i, col="gene.dist", fun="abs", newcol="Distance2gene", cats=c(-1,d,5*10^5), labs=c("close",  "far")))
dseq.ebv <- rec.cat(dt=dseq.ebv, col="gene.dist", fun="abs", newcol="Distance2gene", cats=c(-1,d,5*10^5), labs=c("close", "far"))

##' Tables by EAF 

## Trec:

lapply(trec.ebv, function(i) setorder(TDR.var(dt=i, var= c("Signif", "rej.level", "EAF"), cols= c( "rej.level", "EAF"),  gold= "Gtex-ebv"),
         EAF, rej.level)[,])

## Dseq:

setorder(TDR.var(dt=dseq.ebv, var= c("Signif", "Fdr", "EAF"), cols= c( "Fdr", "EAF"),  gold= "Gtex-ebv"),
         EAF, -Fdr)[,]

##'  Tables by Gene distance, close <= 100,000, far >100,000 up to 5*10^5

## Trec:
lapply(trec.ebv, function(i) TDR.var(dt=i, var= c("Signif", "rej.level", "Distance2gene"), cols= c( "rej.level", "Distance2gene"),  gold= "Gtex-ebv"))

## Dseq:

TDR.var(dt=dseq.ebv, var= c("Signif", "Fdr", "Distance2gene"), cols= c( "Fdr", "Distance2gene"),  gold= "Gtex-ebv")


##' # TDR by Gene expression
counts.f <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'
counts <- fread(counts.f)

## get median expression per gene

med.g <- data.table(Gene_id = counts$gene_id, median.g = apply(counts[,-1], 1, median))
med.g[, Gene_expression:='low'][median.g >= median(median.g), Gene_expression:='high']

## add median expression to datasets:

trec.ebv <- lapply(trec.ebv, function(i) merge(i, med.g[,.(Gene_id,Gene_expression)], by=c("Gene_id")))
dseq.ebv <- merge(dseq.ebv, med.g[,.(Gene_id,Gene_expression)], by=c("Gene_id"))

##' Tables by Gene expression, low < median, high >= median


## Trec:
lapply(trec.ebv, function(i) TDR.var(dt=i, var= c("Signif", "rej.level", "Gene_expression"), cols= c( "rej.level", "Gene_expression"),
        gold= "Gtex-ebv")[order(Gene_expression,rej.level),])

## Dseq:

TDR.var(dt=dseq.ebv, var= c("Signif", "Fdr", "Gene_expression"), cols= c( "Fdr", "Gene_expression"),  gold= "Gtex-ebv")[order(Gene_expression,-Fdr),]

##################################################################################################
##' Comparing Btrecase run with normal or mix of 2 Gaussian normals (Mix1) as priors

## btrecase.n <- comb.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/GT', "^refbias\\.ENSG[0-9]+.*stan.summary.txt")

## btrecase <- mapply(function(i,j)  {
##     dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
##     dt <- add.null(dt)
##     dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
##     return(dt)
## },
## i=c('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/GT',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/MixPrior/GT'),
## j=c("refbias", "rbias"),
## SIMPLIFY=F)

btrecase <- mapply(function(i,j)  {
    dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
    dt <- add.null(dt)
    dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
    return(dt)
    },
    i=c(snakemake@params[['gt_btrecase_normal']], snakemake@params[['gt_btrecase_m1']]),
    j=c("refbias", "rbias"),
    SIMPLIFY=F)
                   
names(btrecase) <- c("normal", "mix")

###################################################
##' Compare Btrecase estimates with different priors

btrecase.comp <- Reduce(function(a,b) {
    dt <- merge(a,b, by=c("Gene_id", "tag", "tag.EAF", "gene.dist"), suffixes=paste0(".", names(btrecase)))
    s <- sapply(names(btrecase), function(i) paste0("null.95.", i))
    dt <- add.signif(dt, s[1],s[2], names(btrecase))
    return(dt)
    },
    btrecase)

cols.b <- c(cols, "null.95")
tab <- tab2bplot(btrecase.comp, colors=c("None"="#999999","normal"="yellow3", "Both"="#D55E00", "mix"="#0072B2"))
btrecase.plot(dt=btrecase.comp[Signif != "None" ,] ,
              x1=paste0(cols.b, ".normal"), 
              x2=paste0(cols.b,  ".mix"),
              ##s=nrow(trec.comp[[i]][Signif != "None" ,]),
              xl="eQTL effect Btrecase normal prior",
              yl="eQTL effect Btrecase mix normal prior",
              col=c("normal", "mix"),
              title=paste0("Btrecase with normal or Gaussian mix prior"),
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=1.5, ymax=4)



################################################################
##' Compare effect size of associations run with btrecase  with or without ASE info

btrecase.comp.l <- reshape(btrecase.comp[,.(Gene_id, tag, log2_aFC_mean.normal, log2_aFC_mean.mix, model.normal, model.mix)],
                           direction="long",
                           varying=list(c( "log2_aFC_mean.normal",  "log2_aFC_mean.mix"), c("model.normal" ,"model.mix" )),
                           v.names=c("log2_aFC_mean", "model"),
                           times=c("normal", "mix"),
                           timevar="prior"
                           )
model <- ggplot(btrecase.comp.l, aes(x=log2_aFC_mean, color=model)) +
    geom_density() +
    facet_grid(prior ~.) +
    ggtitle("Effect size for associations with or without ASE info\n by model and prior")


prior <- ggplot(btrecase.comp.l, aes(x=log2_aFC_mean, color=prior)) +
    geom_density() +
    facet_grid(model ~.) +
    ggtitle("Effect size for associations with or without ASE info\n by model and prior")
#+ fig.width= 6, fig.height=10
plot_grid(model, prior, nrow=2)



############################################################################################
##' Compare trec and btrecase using same prior distribution but only association with ASE info

l3 <- mapply(function(a,b,d,e,f,cols) {
    dt <- merge(a,b, by=c("Gene_id", "tag"), suffixes=c(".trec",".trecase"))
    dt <- add.signif(dt, d, f, c("trec", "trecase"))
    tab <- tab2bplot(dt[model.trecase == "trec-ase",], colors=c("None"="#999999","trec"="yellow3", "Both"="#D55E00", "trecase"="#0072B2"))
    p <- btrecase.plot(dt=dt[Signif != "None" & model.trecase=="trec-ase",] ,
              x1=c(paste0(cols, ".trec") ,d ), 
              x2=c(paste0(cols,  ".trecase"), f ),
              ##s=nrow(trec.comp[[i]][Signif != "None" ,]),
              xl=paste0("eQTL effect Trec ", e),
              yl=paste0("eQTL effect Trecase ", e),
              col=c("trec", "trecase"),
              title=paste0("Trec vs Trecase with ", e),
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-2, xmax=-1, ymin=0.5, ymax=1.5)
    print(p)
    return(dt)
},
a=trec.all[1:2], b=btrecase, d=c("log2_aFC_null",  "null.95.trec"), e=c("normal prior", "mix1 prior"), f=c("null.95", "null.95.trecase"),
MoreArgs=list(cols =c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%") ), SIMPLIFY=F)




btrecase <- lapply(btrecase, function(i) gene.d(i, gt22[, .(gene_id, start,end,chrom)]))

btrecase.rej=seq(0,0.2,0.1)
post.btrecase <- lapply(btrecase, function(i) {
    dt <- rej.recode(a=btrecase.rej, b=i,c=post.level)
    setkey(dt, Gene_id,tag)
    })

###################################################################################
##' Assess normal approximation by comparing significant calls for rection level =0

## trec

post.dt$normal[rej.level==0,.N,.(log2_aFC_null,null.rej)]

lapply(post.dt[paste0("mix",1:3)], function(i) i[rej.level==0, .N, .(null.95, null.rej)])

## trecase

lapply(post.btrecase, function(i) i[model=="trec-ase" & rej.level==0,.N,.(null.95,null.rej)])


##' Normal approximation overestimates the number of significant associations by two-fold

####################################################################################################################
##' Compare trecase with trec using rejection regions of 0 and 0.2 for trec and mix1 prior for both trec and trecase

##+ fig.width= 8, fig.height=21
trec.trecase <- lapply(c(0, 0.2), function(i) {
    dt <- merge(post.dt$mix1[rej.level==i,], post.btrecase$mix, by=c("Gene_id", "tag", "tag.EAF", "gene.dist"), suffixes=c(".trec", ".trecase"))
    dt <- add.signif(dt, "null.rej.trec" , "null.rej.trecase" , col=c("trec","trecase") )

    ## returns all associations but plot only with ASE
    tmp <- dt[model.trecase=="trec-ase",]
    table <- tab2bplot(dt=tmp, var= c("Signif", "rej.level.trecase"), colors=c(None="#999999", trec="yellow3", trecase="#0072B2", Both= "#D55E00"))
    u <- unique(tmp$rej.level.trecase)
    tables <- lapply(u, function(i) {
        table[rej.level.trecase==i,]
    })
    
    gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

    tmp.tab <- data.table(rej.level.trecase=u, grob=gl)

    lab <- paste('r = \u00B1',u)
    names(lab) <-u
    p <- btrecase.plot(dt=tmp[Signif != "None" ,] , x1=c(rep("log2_aFC_mean.trec",3),"null.rej.trec") ,
                       x2=c(rep("log2_aFC_mean.trecase",3), "null.rej.trecase"),
                       xl="eQTL effect Bayesian trec",
                                        #s=nrow(trec.trecase),
                       yl="eQTL effect Btrecase",
                       col=c("trec", "trecase"),
                       title=paste0("eQTL estimates with Trec ", post.level*100, "% PIP out of ",i ,"\n and Trecase by rejection zone for associations with ASE"),
                       title.size=16,
                       axis.title = 14,
                       axis.text=12,
                       legend.title=14,
                       legend.text=12,
                       legend.symbol=5,
                       point.size=3
                       ) + facet_grid(rej.level.trecase~., labeller=labeller(rej.level.trecase=lab))+
        theme(strip.background=element_rect(fill="white"),strip.text.y = element_text(size = 14))+
        geom_custom(data=tmp.tab, aes(grob=grob), x = 0.15, y = 0.75)
    
    print(p)

    ## returns all asscociations
    return(dt)
    })

names(trec.trecase) <- paste("Trec rej level:", c(0,0.2))

##' Number of associations by model
lapply(trec.trecase,
       function(i) rbindlist(lapply(c(0,0.2),
                                    function(j) i[rej.level.trecase == j,
                                                  .N,.(Signif, rej.level.trecase, model.trecase)][order(rej.level.trecase,model.trecase),])))

######################################################################################################
##' Compare trecase with deseq2

btrecase.dsq <- lapply(post.btrecase, function(i){
    dt <- merge(i,dseq, by=c("Gene_id", "tag", "tag.EAF","gene.dist"))
    setkey(dt, p.adj)
    dt <- add.signif(dt, x1="null.rej", x2="null.fdr5", col=c("trecase","dseq") )
    })

#+ fig.width= 8, fig.height=18
l4 <-  mapply(function(y,z,cols){
    k <- y[model=="trec-ase",]
    table <- k[,.N,.(Signif, rej.level)]
    setnames(table, "N", "SNPs")
    table[Signif=="None", color:="#999999"][Signif=="trecase", color:="yellow3"][Signif=="Both", color:="#D55E00"][Signif=="dseq", color:="#0072B2"]
    table[, Signif:=factor(Signif, levels=c("None","trecase", "dseq","Both"))]
    setkey(table,Signif)
    rej <- unique(k$rej.level)
    lab <- paste('r = \u00B1',rej)
    names(lab) <-rej
    tables <- lapply(rej, function(i) {
        table[rej.level==i,]
    })
    
    gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 16,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-1.5, xmax=-0.5, ymin=1.5, ymax=2.5)))
    dt.tab <- data.table(rej.level=rej, grob=gl)
    p <- btrecase.plot(dt=k[Signif != "None",] ,
                       x1=c(rep(cols[1],3),"null.rej") ,
                       x2=c(rep("log2_aFC",3), "null.fdr5"),
                                        #s=nrow(k[Signif != "None" ,]),
                       xl="eQTL effect Bayesian NB",
                       yl="eQTL effect Dseq2",
                       col=c("trecase", "dseq"),
                       title=paste0("DEseq 5% FDR and ", post.level*100 ,"% PIP  by rejection zone\n Trecase with ",z, " prior" ),
                       title.size=16,
                       axis.title = 14,
                       axis.text=12,
                       legend.title=14,
                       legend.text=12,
                       legend.symbol=5,
                       point.size=3
                       ) + facet_grid(rej.level~., labeller=labeller(rej.level=lab))+
        theme(strip.background=element_rect(fill="white"), strip.text.y = element_text(size = 14))+
        geom_custom(data=dt.tab, aes(grob=grob), x = 0.2, y = 0.8)
                                        
    print(p)
    return(tables)
}, y=btrecase.dsq, z=names(btrecase.dsq), SIMPLIFY=F,
MoreArgs=list(cols =c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%") ))

##############################################################################
##' Compare rasqual with deseq2

## rasq <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/output", "ENSG[0-9]+.*txt", full.names=T)
rasq <- list.files(snakemake@params[['rasqual']], "ENSG[0-9]+.*txt", full.names=T)
rasq.header <- paste0(snakemake@params[['rasqual']], "/rasqual.header.txt")

rasqual <- rbindlist(lapply(rasq, format_rasqual, top.hits="no", header=rasq.header))

## rasqual didnt run some genes "Estimated computational time is too long...", remove from output

rasqual <- rasqual[rs_id != "SKIPPED",]

rasqual[, p_adjust:= p.adjust(p, method="BH")]

## transform Fold change to log2aFC
rasqual[, log2_aFC:=log2(Fold_change)]

## add various fdr cut-offs to rasqual

for (i in c(10^(seq(-3,-1,1)), 0.05)){
    rasqual[,eval(paste0("null.fdr",i*100)):= "yes"][p_adjust<=i, paste0("null.fdr",i*100):="no"]
}


## merge with deseq2

rasq.dseq <- merge(dseq, rasqual, by.x=c("Gene_id", "tag"), by.y=c("gene_id", "rs_id"), suffixes=c(".dseq", ".rasqual"))


rasq.dseq.p<- lapply(c(10^(seq(-3,-1,1)), 0.05), function(i){
    dt <- copy(rasq.dseq)
    null.rasq <- paste0("null.fdr",i*100,".rasqual")
    dt <- add.signif(dt, null.rasq, "null.fdr5.dseq", col=c("rasqual","deseq"))
    table  <- tab2bplot(dt, var="Signif", colors=c(None="#999999", rasqual="yellow3", deseq="#0072B2", Both= "#D55E00"))
   
    p <- btrecase.plot(dt=dt[Signif != "None",] , x1=c(rep("log2_aFC.rasqual",3), null.rasq) ,
              x2=c(rep("log2_aFC.dseq",3), "null.fdr5.dseq"),
              #s=nrow(k[Signif != "None" ,]),
              xl="eQTL effect Rasqual",
              yl="eQTL effect Dseq2",
              col=c("rasqual", "deseq"),
              title=paste0("DEseq 5% FDR and\n ", i*100 ,"% FDR in Rasqual" ),
              title.size=16,
              axis.title = 14,
              axis.text=12,
              legend.title=14,
              legend.text=12,
              legend.symbol=5,
              point.size=3
              ) +
        annotation_custom(tableGrob(table[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=0.5, xmax=2, ymin=-10, ymax=-2) 
        
    return(p)
})
#+ fig.width= 10, fig.height=21
plot_grid(plotlist=rasq.dseq.p, nrow=4)

#####################################################################################
##' Rasqual vs Btrec-mix1

#+ fig.width= 10, fig.height=21
rasq.trec.p <- mapply(rasq.bay.p,
                      a=list(rasqual),
                      b=post.dt['mix1'],
                      MoreArgs=list(rej=0.2, suf=c("Rasqual", "Trec"), null=c("null.Fdr", "null.rej")),
                      SIMPLIFY=F)


rbindlist(lapply(rasq.trec.p, function(i) rbindlist(lapply(unique(i[['Fdr_per']]) , function(j) format.tab(i[Fdr_per == j, ], "Rasqual")))))
#######################################
##' Rasqual vs Btrecase

#+ fig.width= 10, fig.height=21
rasq.btrecase.p <- mapply(rasq.bay.p,
                      a=list(rasqual),
                      b=post.btrecase['mix'],
                      MoreArgs=list(rej=c(0.2), suf=c("Rasqual", "Trecase"), null=c("null.Fdr", "null.rej")),
                      SIMPLIFY=F)
rbindlist(lapply(rasq.btrecase.p, function(i) rbindlist(lapply(unique(i[['Fdr_per']]) , function(j) format.tab(i[Fdr_per == j, ], "Rasqual")))))

########################################
##' Rasqual vs Gtex

rasq.l <- reshape(rasqual[,.(gene_id, rs_id, log2_aFC, null.fdr0.1, null.fdr1, null.fdr5, null.fdr10)],
                    direction="long",
                    varying=list( c("null.fdr0.1", "null.fdr1", "null.fdr5", "null.fdr10")),
                    v.names="null.Fdr",
                    times=as.character(c(0.1, 1, 5, 10)),
                    timevar="Fdr_per")
rasq.ebv <- merge(rasq.l, ebv, by.x=c("gene_id", "rs_id"), by.y=c("Gene_id", "SNP"))
rasq.ebv <- add.signif(rasq.ebv, x1="null.Fdr", x2="null", col=c("Rasqual", "Gtex-ebv"))
setnames(rasq.ebv, "gene_id", "Gene_id")
tab <- tab2bplot(rasq.ebv, var= c("Signif", "Fdr_per"), colors=c("None"="#999999","Rasqual"="yellow3", "Both"="#D55E00", "Gtex-ebv"="#0072B2"))

tables.ras <- lapply(unique(tab$Fdr_per), function(i) {
    tab[Fdr_per==i,]
})


gl <- lapply(tables.ras, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(Fdr_per=unique(tab$Fdr_per), grob=gl)


lab <- paste('FDR(%) =', unique(tab$Fdr_per))
names(lab) <- unique(tab$Fdr_per)

##' # Rasqual with different FDR relative to Gtex-EBV
#+ fig.width= 12, fig.height=21
btrecase.plot(dt=rasq.ebv[Signif != "None" ,] , x1=c(rep("log2_aFC",3),"null.Fdr") ,
                   x2=c(rep("slope",3), "null"),
                   xl="eQTL effect Rasqual",
                   yl="eQTL effect Gtex",
                   col=c("Rasqual", "Gtex-ebv"),
                   title="Rasqual vs Gtex-EBV (5% FDR)",
                   title.size=16,
                   axis.title = 14,
                   axis.text=12,
                   legend.title=14,
                   legend.text=12,
                   legend.symbol=5,
                   point.size=3
              ) +  facet_grid(Fdr_per~., labeller=labeller(Fdr_per=lab))+
    theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14)) +
    geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)


rbindlist(lapply(tables.ras, format.tab, "Gtex-ebv" ))

##' Rasqual vs Gtex-ebv for the same associations tested by trec
rasq.ebv.sub <- merge(trec[,.(Gene_id, tag, tag.EAF)], rasq.ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "rs_id"))
tab <- tab2bplot(rasq.ebv.sub, var= c("Signif", "Fdr_per"), colors=c("None"="#999999","Rasqual"="yellow3", "Both"="#D55E00", "Gtex-ebv"="#0072B2"))

tables.ras.sub <- lapply(unique(tab$Fdr_per), function(i) {
    tab[Fdr_per==i,]
})

rbindlist(lapply(tables.ras.sub, format.tab, "Gtex-ebv" ))


############################################
##' Check prior fit to Gtex-EBV data

## Randomly select 10^6 rows. Select rows every 150 aiming to get
## independent SNPs. File ordered by gene and SNP

## ebv2 <- fread(cmd="zcat /mrc-bsu/scratch/ev250/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.allpairs.txt.gz | awk 'NR == 1 || NR % 150 == 0' | head -n 1000000")
ebv2 <- fread(cmd=paste("zcat", snakemake@input[['gtex']],
                        "| awk 'NR == 1 || NR % 150 == 0' | head -n 1000000" ))


##' Mixture of 2 Gaussians: mix1

mixmdl2 = normalmixEM(ebv2$slope, lambda=c(0.05,0.5),mu=c(0,0), mean=c(0,0))

summary(mixmdl2)

qq.mix(data=ebv2$slope, normals=list(mix=c(0.1, 0.9), sd=mixmdl2$sigma, mu=mixmdl2$mu))

##' Mixture 2

qq.mix(data=ebv2$slope, normals=list(mix=mixmdl2$lambda, sd=mixmdl2$sigma, mu=mixmdl2$mu))


##' Mixture of 3 Gaussians: mix3

mixmdl = normalmixEM(ebv2$slope,k=3, lambda=c(0.05,0.5,0.45),mu=c(0,0,0), mean=c(0,0,0))

summary(mixmdl)

qq.mix(data=ebv2$slope, normals=list(mix=mixmdl$lambda, sd=mixmdl$sigma, mu=mixmdl$mu))



##' Mixture of 2 Gaussians for spike-slab (Chis)

qq.mix(data=ebv2$slope, normals=list(mix=c(0.96, 0.44), sd=c(0,sqrt(0.0835)), mu=c(0,0)))

##' Mixture of 3 Gaussians for spike-slab (Chis)

qq.mix(data=ebv2$slope, normals=list(mix=c(0.83, 0.15, 0.017), sd=c(0,sqrt(0.00984), sqrt(0.1592)), mu=c(0,0,0)))
 

##' Compare prior with sd 0 and prior with sd 0.001

p1 <- list(mix=c(0.83, 0.15, 0.017), sd=c(0,sqrt(0.00984), sqrt(0.1592)), mu=c(0,0,0))
p2 <- list(mix=c(0.83, 0.15, 0.017), sd=c(0.001,sqrt(0.00984), sqrt(0.1592)), mu=c(0,0,0))


qq.p(p1,p2,n=nrow(ebv2))

qq.p(list(mix=c(0.96, 0.04), sd=c(0, sqrt(0.0835)), mu=c(0,0)),
     list(mix=c(0.96, 0.04), sd=c(0.001,sqrt(0.0835)), mu=c(0,0)),
     n=nrow(ebv2))
