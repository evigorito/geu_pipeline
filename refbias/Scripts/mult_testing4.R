#' ---
#' title: Choose prior and set threshold for making significant calls with Bayesian model
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


#################################################################################
##' # Compare associations using Trec and Btrecase  with old and new priors #####
#################################################################################

##########################
## Open trec and format ##
##########################
## trec <- comb.files(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/GT',pattern="trec.stan.summary.txt")

trec <- comb.files(path=snakemake@params[['trec']], pattern="trec.stan.summary.txt")

## Add EAF

## le.file <-  "/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz"
le.file <- snakemake@input[['lefile']]

EAF <-  snp.eaf(le.file, unique(trec$tag))
setnames(EAF, "eaf", "tag.EAF")

trec <- merge(trec, EAF, by.x="tag", by.y="snp")


## new priors

trec.m <- lapply(paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrec/" , c("SpikeMixV3_2", "SpikeMixV3_3")),
                 function(i) comb.files(path=i,pattern=".stan.summary.txt"))
names(trec.m) <- c("SpikeMixV3_2", "SpikeMixV3_3")

trec.m <- lapply(snakemake@params[['trec_other']] , function(i) comb.files(i, pattern=".stan.summary.txt"))

trec.m  <- lapply(trec.m, add.null)

names(trec.m) <- basename(snakemake@params[['trec_other']])

## add gene distance to trec

## gene.coord <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt")

gene.coord <- fread(snakemake@input[['geneStEnd']]) 

## select genes in chrom 22

gt22 <- gene.coord[chrom==22,]

## add tag distance to gene (closest to start or end)

trec.m <- lapply(trec.m, function(i) gene.d(i, gt22[, .(gene_id, start,end,chrom)]))

## remove bad runs

trec.m <- lapply(trec.m , function(i) i[Rhat <= 1.01,])

## merge after removing bad runs

trec.comp <- lapply(trec.m, function(i) merge(trec, i, by=c("Gene_id", "tag"), suffixes=c(".norm", ".mix")))

trec.comp <- lapply(trec.comp, function(i) add.signif(i, "log2_aFC_null", "null.95", c("Normal", "Mix")))


cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%")


#############################
## Open Btrecase and format ##
##############################

## btrecase <- mapply(function(i,j)  {
##     dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
##     dt <- add.null(dt)
##     dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
##     return(dt)
## },
## i=c('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/GT',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/GT',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_3/GT'
##     ),
## j=c("refbias", rep("rbias",2)),
## SIMPLIFY=F)

## names(btrecase) <- basename(dirname(c('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/GT',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/GT',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_3/GT'
##     )))

btrecase <- mapply(function(i,j)  {
    dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
    dt <- add.null(dt)
    dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
    return(dt)
    },
    i=c(snakemake@params[['gt_btrecase_normal']], snakemake@params[['gt_btrecase_mix']]),
    j=c("refbias", rep("rbias",2)),
    SIMPLIFY=F)


names(btrecase) <- basename(dirname(c(snakemake@params[['gt_btrecase_normal']], snakemake@params[['gt_btrecase_mix']])))

btrecase.comp <- lapply(btrecase[2:length(btrecase)],function(i) {
    dt <- merge(btrecase[[1]], i, by=c("Gene_id", "tag","tag.EAF", "gene.dist"), suffixes=c(".norm", ".mix"))
    
    dt <- add.signif(dt, "null.95.norm", "null.95.mix", c("Normal", "Mix"))
    return(dt)})

cols.b <- c(cols, "null.95")




##########################################################
##'  Trec
l1 <- lapply(seq_along(trec.comp), function (i) btrecase.plot(dt=trec.comp[[i]][Signif != "None" ,] ,
                                                        x1=c(paste0(cols, ".norm") ,"log2_aFC_null" ), 
                                                        x2=c(paste0(cols,  ".mix"), "null.95" ),
                                                        #s=nrow(trec.comp[[i]][Signif != "None" ,]),
                                                        xl="eQTL effect normal prior",
                                                        yl="eQTL effect mix prior",
                                                        col=c("Normal", "Mix"),
                                                        title=paste0("Trec with normal vs\n Gaussian ",names(trec.m)[i] ," prior"),
                                                        title.size=16,
                                                        axis.title = 14,
                                                        axis.text=12,
                                                        legend.title=14,
                                                        legend.text=12,
                                                        legend.symbol=5,
                                                        point.size=3
                                                        ))
##+ fig.width= 8.3, fig.height=4.78
plot_grid(plotlist=l1, ncol=length(trec.m))
lapply(trec.comp, function(i) i[,.N,Signif])

#################################################################
##' Trecase

btrecase.tab <- lapply(names(btrecase.comp), function(i) {
    tab <- tab2bplot(btrecase.comp[[i]], colors=c("None"="#999999","Normal"="yellow3", "Both"="#D55E00", "Mix"="#0072B2"))
    p <- btrecase.plot(dt=btrecase.comp[[i]][Signif != "None" ,] ,
                       x1=paste0(cols.b, ".norm"), 
                       x2=paste0(cols.b,  ".mix"),
                       ##s=nrow(trec.comp[[i]][Signif != "None" ,]),
                       xl="eQTL effect Btrecase normal prior",
                       yl="eQTL effect Btrecase mix normal prior",
                       col=c("Normal", "Mix"),
                       title=paste0("Btrecase with normal or Gaussian mix prior ", i),
                       title.size=16,
                       axis.title = 14,
                       axis.text=12,
                       legend.title=14,
                       legend.text=12,
                       legend.symbol=5,
                       point.size=3
                       ) +
        annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-3, xmax=-1.5, ymin=1.5, ymax=3)
    print(p)
    return(tab)
})



## ################################################################
##' Compare effect size of associations run with btrecase  with or without ASE info

btrecase.com.l <- lapply(names(btrecase.comp), function(i) {
    dt <-  reshape(btrecase.comp[[i]][,.(Gene_id, tag, log2_aFC_mean.norm, log2_aFC_mean.mix, model.norm, model.mix)],
                           direction="long",
                           varying=list(c( "log2_aFC_mean.norm",  "log2_aFC_mean.mix"), c("model.norm" ,"model.mix" )),
                           v.names=c("log2_aFC_mean", "model"),
                           times=c("normal", "mix"),
                           timevar="prior"
                           )
    model <- ggplot(dt, aes(x=log2_aFC_mean, color=model)) +
        geom_density() +
        facet_grid(prior ~., scales="free") +
        ggtitle(paste0("Effect size for associations with or without ASE info\n by model and prior ", i))

    

    prior <- ggplot(dt, aes(x=log2_aFC_mean, color=prior)) +
        geom_density() +
        facet_grid(model ~., scale="free") +
        ggtitle(paste0("Effect size for associations with or without ASE info\n by model and prior ", i))
                                        
    print(plot_grid(model, prior, ncol=1))
    })


###########################################################################################################################################

##' Calculate PIP by rejection zone in trec.m and add relevant columns

## Bayesian trec : using normal approximation calculate proportion of
## posterior out of rejection zone. If the mean of the posterior is
## within the rejection zone (-r,r) I set the posterior out of the
## rejection zone as 0% as I dont want to call any of those
## associations significant. If the rejection zone is narrow I could
## have a high % of the posterior out of the zone. The I coded a
## variable, null.rej "yes" if the % of the posterior out of the
## rejection zone is below a threshold I define (post.level) and "no"
## otherwise.



rej<-c(0,0.06, 0.08) ## this is for "SpikeMixV3_2", "SpikeMixV3_3", based on 95% and 99% of prior within +/- rej



## Trec
post.dt <- lapply(trec.m, function(i){
    rbindlist(mapply(function(x,y,z) {
        dt <- rej.recode(a=x,b=y,c=z)
        dt[, post.level:=z]
        setkey(dt, Gene_id, tag)
    }, x=list(c(0,0.06), c(0, 0.08)), z=list(0.95, 0.99), MoreArgs=list(y=i), SIMPLIFY=F))
    })

## Trecase
post.btrecase <-lapply(btrecase[2:3], function(i){
    rbindlist(mapply(function(x,y,z) {
        dt <- rej.recode(a=x,b=y,c=z)
        dt[, post.level:=z]
        setkey(dt, Gene_id, tag)
    }, x=list(c(0,0.06), c(0, 0.08)), z=list(0.95, 0.99), MoreArgs=list(y=i), SIMPLIFY=F))
    })




######################################################################################################################################
##' # Use Gtex-ebv as gold standard and compare to Trec, 
##' # DEseq2, Rasqual and Trecase. Choose FDR in Dseq that gives a similar number of
##' # associations as Gtex-ebv
########################################################################################################################################

##' Open and format datasets: Gtex-ebv, DEseq and Rasqual

########################################
## Look at Gtex ebv for chromosome 22 ##
########################################

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

###########
## DEseq ##
##########
## Open and format DEseq2 output (run by Wei-Yu)

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

## DEseq at various FDR to Gtex EBV cell lines at 5% FDR

dseq.l <- rbindlist(lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tmp <- dseq[, null.fdr:="yes"][p.adj<=i, null.fdr:="no"]
    tmp[, Fdr:= i]
    tmp <- tmp[, .(Gene_id, tag, tag.EAF, log2FoldChange,log2_aFC, p.adj, null.fdr, Fdr, gene.dist)]
    return(tmp)
    }))

#############
## Rasqual ##
#############

## rasq <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/output", "ENSG[0-9]+.*txt", full.names=T)
## rasq.header <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/output", "rasqual.header.txt", full.names=T)


rasq <- list.files(snakemake@params[['rasqual']], "ENSG[0-9]+.*txt", full.names=T)
rasq.header <- paste0(snakemake@params[['rasqual']], "/rasqual.header.txt")

rasqual <- rbindlist(lapply(rasq, format_rasqual, top.hits="no", header=rasq.header))

## rasqual didnt run some genes "Estimated computational time is too long...", remove from output

rasqual <- rasqual[rs_id != "SKIPPED",]

## select associations run in Trec, no tagging was done in rasqual

rasqual <-  merge(trec[,.(Gene_id, tag, tag.EAF)], rasqual, by.x=c("Gene_id", "tag"), by.y=c("gene_id", "rs_id"))

rasqual[, p_adjust:= p.adjust(p, method="BH")]

## transform Fold change to log2aFC
rasqual[, log2_aFC:=log2(Fold_change)]

## add various fdr cut-offs to rasqual

for (i in c(10^(seq(-3,-1,1)), 0.05)){
    rasqual[,eval(paste0("null.fdr",i*100)):= "yes"][p_adjust<=i, paste0("null.fdr",i*100):="no"]
}

###################################################
##' Compare Bayesian trec to Gtex EBV cell lines
###################################################



## merge post.dt with ebv but on the same associations as rasqual

trec.ebv <- lapply(post.dt, function(i) {
    dt <- merge(i, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
    dt <- merge(dt, rasqual[, .(Gene_id, tag, Chrom)], by=c("Gene_id", "tag"))
    dt <- add.signif(dt, x1="null.rej", x2="null", col=c("trec", "Gtex-ebv"))
})


##+ fig.width= 9.54, fig.height=7.15

tables.trec <- plot_tab(a=trec.ebv, fac=list(rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))), colors=c(None="#999999", trec="yellow3", `Gtex-ebv`="#0072B2", Both= "#D55E00"),
                        title = "Gtex-ebv at 5%FDR vs Trec by PIP and rejection level\n using prior ",
                        xpos=0.25, ypos=0.7,
                        x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                        x2=c(rep("slope",3), "null"),
                                        #s=50000,
                        xl="eQTL effect Bayesian trec",
                        yl="eQTL effect Gtex",
                        col=c("trec", "Gtex-ebv"),
                        title.size=16,
                        axis.title = 14,
                        axis.text=12,
                        legend.title=14,
                        legend.text=12,
                        legend.symbol=5,
                        point.size=3
                        )

lapply(tables.trec, function(i) rbindlist(lapply(i, format.tab , "Gtex-ebv")))

#########################################
##'  Compare DEseq to Gtex EBV 
#########################################

dseq.ebv <- merge(dseq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))

## make it comparable with rasqual

dseq.ebv <- merge(dseq.ebv, rasqual[, .(Gene_id, tag, Chrom)], by=c("Gene_id", "tag"))

dseq.ebv <-  add.signif(dseq.ebv, x1="null.fdr", x2="null", col=c("DEseq","Gtex-ebv") )

tab <- tab2bplot(dseq.ebv, var= c("Signif", "Fdr"), colors=c("None"="#999999","DEseq"="yellow3", "Both"="#D55E00", "Gtex-ebv"="#0072B2"))

tables.dseq <- lapply(c(0.1, 0.05, 0.01, 0.001), function(i) {
    tab[Fdr==i,]
})


gl <- lapply(tables.dseq, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(Fdr=c(0.1, 0.05, 0.01, 0.001), grob=gl)


lab <- paste('FDR(%) =', c(0.1, 0.05, 0.01, 0.001)*100)
names(lab) <- c(0.1, 0.05, 0.01, 0.001)

##'  DEseq with different FDR relative to Gtex-EBV
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


## DESeq:  Significant and total associations by FDR

rbindlist(lapply(tables.dseq, format.tab, "Gtex-ebv" ))


##  I carry on with DEseq2 model using fdr 0.01 and 0.001 which are the best matches to Gtex-ebv


###################
##' Rasqual vs Gtex
####################
rasq.l <- reshape(rasqual[,.(Gene_id, tag, log2_aFC, null.fdr0.1, null.fdr1, null.fdr5, null.fdr10)],
                    direction="long",
                    varying=list( c("null.fdr0.1", "null.fdr1", "null.fdr5", "null.fdr10")),
                    v.names="null.Fdr",
                    times=as.character(c(0.1, 1, 5, 10)),
                    timevar="Fdr_per")
rasq.ebv <- merge(rasq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
rasq.ebv <- add.signif(rasq.ebv, x1="null.Fdr", x2="null", col=c("Rasqual", "Gtex-ebv"))

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

##' Carry on with 0.01 and 0.001 FDR

###################
##' Trecase vs Gtex
###################


trecase.ebv <- lapply(post.btrecase, function(i){
    dt <- merge(i, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
    dt <- merge(dt, rasqual[, .(Gene_id, tag, Chrom)], by=c("Gene_id", "tag"))
    dt <- add.signif(dt, x1="null.rej", x2="null", col=c("trecase", "Gtex-ebv"))
})


##+ fig.width= 9.78, fig.height=7.97

tables.trecase <- plot_tab(a=trecase.ebv, fac=list(rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))), colors=c(None="#999999", trecase="yellow3", `Gtex-ebv`="#0072B2", Both= "#D55E00"),
                        title = "Gtex-ebv at 5%FDR vs Trecase by PIP and rejection level\n using prior ",
                        xpos=0.25, ypos=0.75,
                        x1=c(rep("log2_aFC_mean",3),"null.rej") ,
                        x2=c(rep("slope",3), "null"),
                                        #s=50000,
                        xl="eQTL effect trecase",
                        yl="eQTL effect Gtex",
                        col=c("trecase", "Gtex-ebv"),
                        title.size=16,
                        axis.title = 14,
                        axis.text=12,
                        legend.title=14,
                        legend.text=12,
                        legend.symbol=5,
                        point.size=3
                        )

lapply(tables.trecase, function(i) rbindlist(lapply(i, format.tab , "Gtex-ebv")))


##########################################################################
##' # Compare trec, trecase and rasqual  with Dseq2  model ###############
######################################################################### 

## merge post.dt and dseq
fdr <- c(0.001, 0.01)
btrec.dsq <- lapply(post.dt, function(i) {
    rbindlist(lapply(fdr, function(j) {
        dt <- merge(i,dseq.l[Fdr %in% j,], by=c("Gene_id", "tag", "tag.EAF", "gene.dist"))
        setkey(dt, p.adj)
        dt <- add.signif(dt, x1="null.rej", x2="null.fdr", col=c("trec","dseq") )
    }))
    })

###################################
##' Trec vs Dseq2  model at 0.1% FDR
###################################

lab <- paste('r = \u00B1',rej)
names(lab) <- rej

#+ fig.width= 7.4, fig.height=6.7

l2 <- lapply(1:2, function(i) plot_tab(a=btrec.dsq[i], fac=list(Fdr=rep(0.001,4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
               colors=c(None="#999999", trec="yellow3", dseq="#0072B2", Both= "#D55E00"),
               title = "DEseq at 0.1 %FDR vs Trec by PIP\n and rejection level\n using prior ",
               var=setNames(0.001, "Fdr"),
               xpos=0.25, ypos=0.75,
               facet.fac=c("rej.level", "post.level"),
               x1=c(rep(cols[1],3),"null.rej"),
               x2=c(rep("log2_aFC",3), "null.fdr"),
               xl="eQTL effect Trec",
               yl="eQTL effect Dseq2",
               col=c("trec", "dseq"),
               title.size=12,
               axis.title = 10,
               axis.text=10,
               legend.title=12,
               legend.text=10,
               legend.symbol=3,
               tab.size=8,
               point.size=1.5))


##' Comparing Trec with DEseq by prior and Fdr

lapply(fdr, function(j) lapply(only_tab(a=btrec.dsq, fac=list(Fdr=rep(j, 4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2)))),
       function(i) rbindlist(lapply(i, format.tab,"dseq"))))

## merge post.btrecase and dseq

btrecase.dsq <- lapply(post.btrecase, function(i) {
    rbindlist(lapply(fdr, function(j) {
        dt <- merge(i,dseq.l[Fdr %in% j,], by=c("Gene_id", "tag", "tag.EAF", "gene.dist"))
        setkey(dt, p.adj)
        dt <- add.signif(dt, x1="null.rej", x2="null.fdr", col=c("trecase","dseq") )
    }))
    })

#######################################
##' Trecase vs Dseq2  model at 0.1% FDR
#######################################

#+ fig.width= 7.4, fig.height=6.7

btrecase.dsq.p <- lapply(1:2, function(i) plot_tab(a=btrecase.dsq[i], fac=list(Fdr=rep(0.001,4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
                                                   colors=c(None="#999999", trecase="yellow3", dseq="#0072B2", Both= "#D55E00"),
                                                   var=setNames(0.001, "Fdr"),
                                                   title = "DEseq at 0.1 %FDR vs Trecase by PIP\n and rejection level\n using prior ",
                                                   xpos=0.25, ypos=0.75,
                                                   facet.fac=c("rej.level", "post.level"),
                                                   x1=c(rep(cols[1],3),"null.rej"),
                                                   x2=c(rep("log2_aFC",3), "null.fdr"),
                                                   xl="eQTL effect Trec",
                                                   yl="eQTL effect Dseq2",
                                                   col=c("trecase", "dseq"),
                                                   title.size=16,
                                                   axis.title = 14,
                                                   axis.text=12,
                                                   legend.title=14,
                                                   legend.text=12,
                                                   legend.symbol=5,
                                                   tab.size=12,
                                                   point.size=3))

##' Comparing Trecase with DEseq by prior and Fdr

lapply(fdr, function(j) lapply(only_tab(a=btrecase.dsq, fac=list(Fdr=rep(j, 4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2)))),
       function(i) rbindlist(lapply(i, format.tab,"dseq"))))

#################################################
##' Compare rasqual with deseq2 
#################################################


## merge with deseq2

rasq.dseq <- merge(dseq.l[Fdr==fdr[1],], rasq.l, by=c("Gene_id", "tag"),  suffixes=c(".dseq", ".rasqual"))

rasq.dseq <- rbindlist(lapply(fdr, function(i) {
    dt <- merge(dseq.l[Fdr == i,], rasq.l[Fdr_per == i*100,], by=c("Gene_id", "tag"),
                suffixes=c(".dseq", ".rasqual"))
    dt <- add.signif(dt, "null.fdr", "null.Fdr", col=c("deseq", "rasqual"))
    return(dt)
}
))

rasq.dseq.tab <- tab2bplot(rasq.dseq, var=c("Signif", "Fdr"), colors=c(None="#999999", deseq="yellow3", rasqual="#0072B2", Both= "#D55E00"))

tables.rasq.dseq <- lapply(fdr, function(i) {
    rasq.dseq.tab[Fdr==i,]
})


gl <- lapply(tables.rasq.dseq, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

dt.tab <- data.table(Fdr=fdr, grob=gl)


lab <- paste('FDR =', fdr)
names(lab) <- fdr

#############################
##'  DEseq vs rasqual by FDR
#############################

#+ fig.width= 6, fig.height=8
btrecase.plot(dt=rasq.dseq[Signif != "None" ,] , x1=c(rep("log2_aFC.dseq",3),"null.fdr") ,
                   x2=c(rep("log2_aFC.rasqual",3), "null.Fdr"),
                   xl="eQTL effect DEseq",
                   yl="eQTL effect Rasqual",
                   col=c("deseq", "rasqual"),
                   title="DEseq2 vs Rasqual by FDR",
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


rbindlist(lapply(tables.rasq.dseq, format.tab,"deseq"))


##############################################################################################
##' Compare trec and btrecase using same prior distribution but only association with ASE info
##############################################################################################

##' Compare trecase with trec using rejection regions and post.level

##+ fig.width= 8.8, fig.height=8.3

trec.trecase.ase <- mapply(function(a,b) {
    dt <- merge(a,b[model=="trec-ase",], by=c("Gene_id", "tag", "tag.EAF", "gene.dist", "post.level", "rej.level"), suffixes=c(".trec", ".trecase"))
    dt <- add.signif(dt, "null.rej.trec" , "null.rej.trecase" , col=c("trec","trecase") )
    }, a=post.dt, b=post.btrecase, SIMPLIFY=F)


tables.trec.trecase <-  plot_tab(trec.trecase.ase,
                                 fac=list(rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
                                 colors=c(None="#999999", trec="yellow3", trecase="#0072B2", Both= "#D55E00"),
                                 title = "Trec vs trecase by PIP and rejection level\n using prior ",
                                 xpos=0.25, ypos=0.75,
                                 x1=c(rep("log2_aFC_mean.trec",3),"null.rej.trec") ,
                                 x2=c(rep("log2_aFC_mean.trecase",3), "null.rej.trecase"),
                                        #s=50000,
                                 xl="eQTL effect trec",
                                 yl="eQTL effect trecase",
                                 col=c("trec", "trecase"),
                                 title.size=16,
                                 axis.title = 14,
                                 axis.text=12,
                                 legend.title=14,
                                 legend.text=12,
                                 legend.symbol=5,
                                 point.size=3
                                 )


##' Number of associations by model

lapply(tables.trec.trecase, function(i) rbindlist(lapply(i, format.tab, "trec" )))


################################
##' Rasqual vs Btrec
################################
#+ fig.width= 7.4, fig.height=6.7

rasq.trec <- lapply(post.dt, function(i) {
    rbindlist(lapply(fdr, function(j) {
        dt <- merge(i,rasq.l[Fdr_per %in% as.character(j*100),], by=c("Gene_id", "tag"))
        dt <- add.signif(dt, x1="null.rej", x2="null.Fdr", col=c("trec","rasqual") )
    }))
})


rasq.trec.p <- lapply(1:2, function(i) plot_tab(a=rasq.trec[i], fac=list(Fdr_per=rep(0.1,4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
                                                colors=c(None="#999999", trec="yellow3", rasqual="#0072B2", Both= "#D55E00"),
                                                var=setNames(0.1, "Fdr_per"),
                                                title = "Rasqual at 0.1 %FDR vs Trec by PIP\n and rejection level\n using prior ",
                                                xpos=0.25, ypos=0.75,
                                                facet.fac=c("rej.level", "post.level"),
                                                x1=c(rep(cols[1],3),"null.rej"),
                                                x2=c(rep("log2_aFC",3), "null.Fdr"),
                                                xl="eQTL effect Trec",
                                                yl="eQTL effect Rasqual",
                                                col=c("trec", "rasqual"),
                                                title.size=16,
                                                axis.title = 14,
                                                axis.text=12,
                                                legend.title=14,
                                                legend.text=12,
                                                legend.symbol=5,
                                                tab.size=10,
                                                point.size=3))

##' Comparing Rasqual to Trec by prior and Fdr

lapply(fdr*100, function(j) lapply(only_tab(a=rasq.trec, fac=list(Fdr_per=rep(j, 4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2)))),
       function(i) rbindlist(lapply(i, format.tab,"rasqual"))))

#####################################
##' Rasqual vs Btrecase ASE only
#####################################


#+ fig.width= 7.4, fig.height=6.7

fdr=c(0.001, 0.01, 0.05, 0.1)

rasq.trecase <- lapply(post.btrecase, function(i) {
    rbindlist(lapply(fdr, function(j) {
        dt <- merge(i[model=="trec-ase",] ,rasq.l[Fdr_per %in% as.character(j*100),], by=c("Gene_id", "tag"))      
        dt <- add.signif(dt, x1="null.rej", x2="null.Fdr", col=c("trecase","rasqual") )
    }))
})


rasq.trecase.p <-lapply(1:2, function(i)  plot_tab(a=rasq.trecase[i], fac=list(Fdr_per=rep(0.1,4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
                                                   colors=c(None="#999999", trecase="yellow3", rasqual="#0072B2", Both= "#D55E00"),
                                                   var=setNames(0.1, "Fdr_per"),
                                                   title = "Rasqual at 0.1 %FDR vs Trecase by PIP\n and rejection level\n using prior ",
                                                   xpos=0.25, ypos=0.75,
                                                   facet.fac=c("rej.level", "post.level"),
                                                   x1=c(rep(cols[1],3),"null.rej"),
                                                   x2=c(rep("log2_aFC",3), "null.Fdr"),
                                                   xl="eQTL effect Trecase",
                                                   yl="eQTL effect Rasqual",
                                                   col=c("trecase", "rasqual"),
                                                   title.size=16,
                                                   axis.title = 14,
                                                   axis.text=12,
                                                   legend.title=14,
                                                   legend.text=12,
                                                   legend.symbol=5,
                                                   point.size=3))

##' Comparing Rasqual to Trecase by prior and Fdr

lapply(fdr*100, function(j) lapply(only_tab(a=rasq.trecase, fac=list(Fdr_per=rep(j, 4), rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2)))),
       function(i) rbindlist(lapply(i, format.tab,"rasqual"))))


#######################################################################################
##' Hidden-GT vs Trec and Btrecase
#######################################################################################


## read file with old output with tags matching GT and hidden-GT
## match.tags <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/results/refbias.gt.rna.p054.txt")

match.tags <- fread(snakemake@input[['gt_rna_sum']])

match.tags <- match.tags[,.(Gene_id, tag.rna,tag.gt, op.dir)]

## ## change sign of effect size of op.dir=="yes"
## nogt <- mapply(function(i,j)  {
##     dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
##     dt <- add.null(dt)
##     dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
##     dt <- merge(match.tags, dt, by.x=c("Gene_id","tag.rna"), by.y=c("Gene_id", "tag"))
##     ##dt <- dt[!Gene_id %in% c( "ENSG00000093072", "ENSG00000100364", "ENSG00000241973"),]
##     dt[op.dir=="yes", log2_aFC_mean:= -log2_aFC_mean]
##     return(dt)
## },
## i=c('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/RNA',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_3/RNA'
##     ),
## j=c( rep("rbias",2)),
## SIMPLIFY=F)

## names(nogt) <- basename(dirname(c('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/RNA',
##     '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_3/RNA'
##     )))

nogt <- mapply(function(i,j)  {
    dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
    dt <- add.null(dt)
    dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
    dt <- merge(match.tags, dt, by.x=c("Gene_id","tag.rna"), by.y=c("Gene_id", "tag"))   
    dt[op.dir=="yes", log2_aFC_mean:= -log2_aFC_mean]
    return(dt)
    },
    i=snakemake@params[['nogt_mix']],
    j=c(rep("rbias",2)),
    SIMPLIFY=F)


names(nogt) <- basename(dirname(snakemake@params[['nogt_mix']]))

post.nogt <-lapply(nogt, function(i){
    rbindlist(mapply(function(x,y,z) {
        dt <- rej.recode(a=x,b=y,c=z)
        dt[, post.level:=z]
        setkey(dt, Gene_id, tag.gt)
    }, x=list(c(0,0.06), c(0, 0.08)), z=list(0.95, 0.99), MoreArgs=list(y=i), SIMPLIFY=F))
    })

###################################
##' Compare observed vs hidden-GT
##################################


trecase.nogt <- mapply(function(a,b) {
    dt <- merge(a ,b, by.x=c("Gene_id", "tag", "post.level", "rej.level"),
                by.y=c("Gene_id", "tag.gt", "post.level", "rej.level"),
                suffixes=c(".trecase", ".noGT"))
    dt <- add.signif(dt, "null.rej.trecase" , "null.rej.noGT" , col=c("obs_GT","hidden_GT") )
    }, a=post.btrecase, b=post.nogt, SIMPLIFY=F)

#+ fig.width= 8.8, fig.height=8.3

tables.trecase.nogt <-  plot_tab(trecase.nogt,
                                 fac=list(rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
                                 colors=c(None="#999999", obs_GT="yellow3", hidden_GT="#0072B2", Both= "#D55E00"),
                                 title = "Observed vs hidden GT by PIP and rejection level\n using prior ",
                                 xpos=0.25, ypos=0.75,
                                 x1=c(paste0(cols, ".trecase"),"null.rej.trecase") ,
                                 x2=c(paste0(cols,".noGT"), "null.rej.noGT"),
                                        #s=50000,
                                 xl="eQTL effect obs-GT",
                                 yl="eQTL effect hidden-GT",
                                 col=c("obs_GT","hidden_GT"),
                                 title.size=16,
                                 axis.title = 14,
                                 axis.text=12,
                                 legend.title=14,
                                 legend.text=12,
                                 legend.symbol=5,
                                 point.size=1.5
                                 )


##' Number of associations by model

lapply(tables.trecase.nogt, function(i) rbindlist(lapply(i, format.tab, "obs_GT" )))


######################
##' Trec vs hidden-GT#
######################

trec.nogt <-  mapply(function(a,b) {
    dt <- merge(a,b, by.x=c("Gene_id", "tag", "post.level", "rej.level"),
                by.y=c("Gene_id", "tag.gt", "post.level", "rej.level"),
                suffixes=c(".trec", ".noGT"))
    dt <- add.signif(dt, "null.rej.trec" , "null.rej.noGT" , col=c("obs_GT","hidden_GT") )
    }, a=post.dt, b=post.nogt, SIMPLIFY=F)

#+ fig.width= 8.8, fig.height=8.3

tables.trec.nogt <-  plot_tab(trec.nogt,
                                 fac=list(rej.level=c(0, 0.06, 0, 0.08), post.level=c(rep(0.95,2), rep(0.99,2))),
                                 colors=c(None="#999999", obs_GT="yellow3", hidden_GT="#0072B2", Both= "#D55E00"),
                                 title = "Trec vs hidden GT by PIP and rejection level\n using prior ",
                                 xpos=0.25, ypos=0.75,
                                 x1=c(paste0(cols, ".trec"),"null.rej.trec") ,
                                 x2=c(paste0(cols,".noGT"), "null.rej.noGT"),
                                        #s=50000,
                                 xl="eQTL effect obs-GT",
                                 yl="eQTL effect hidden-GT",
                                 col=c("obs_GT","hidden_GT"),
                                 title.size=16,
                                 axis.title = 14,
                                 axis.text=12,
                                 legend.title=14,
                                 legend.text=12,
                                 legend.symbol=5,
                                 point.size=3
                                 )


##' Number of associations by model

lapply(tables.trec.nogt, function(i) rbindlist(lapply(i, format.tab, "obs_GT" )))


###################################################################################################################
##' Repeating analysis with trec and trecase with and without genotypes using posterior instead of normal approx.
###################################################################################################################


######################
##' Obs vs hidden-GT #
######################

## select spikeSlab priors
btrecase.m <- btrecase[2:3]




p <- lapply(c("null.95", "null.99"), function(i) {
    mapply(merge.plot,
           a=btrecase.m,
           prior=names(btrecase.m),
           b=nogt,
           MoreArgs=list(null=i,byx=c("Gene_id", "tag") , byy=c("Gene_id", "tag.gt"), suffixes=c(".gt", ".nogt"),
                         colsig=c("obs_GT", "hidden_GT"), siglevel=gsub("null.", "",i), xl="eQTL effect observed-GT",
                         yl="eQTL effect hidden-GT", title=paste0("PIP = ", gsub("null.","", i), "%"),
                         size.tab=8, xmin=-1,  xmax=-.5, ymin=.3, ymax=.5),
           SIMPLIFY=F )})

#+ fig.width= 7, fig.height=5.5
plot_grid(plotlist=c(p[[1]], p[[2]]), nrow=2)
