#' ---
#' title: Look at associations with SpikeMixV3_2 prior
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



source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")
source('/home/ev250/Cincinatti/Functions/various.R')
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_assoc.R")

###################################
##' Open and format trecase with GT
####################################

## gene.coord <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt")

gene.coord <- fread(snakemake@input[['geneStEnd']]) 

## select genes in chrom 22

gt22 <- gene.coord[chrom==22,]

## trecase <- comb.files('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/GT', paste0("^", "rbias","\\.ENSG[0-9]+.*stan.summary.txt"))

trecase <- comb.files(snakemake@params[['gt_btrecase_mix']], paste0("^", "rbias","\\.ENSG[0-9]+.*stan.summary.txt"))


## add tag distance to gene (closest to start or end)
trecase <- gene.d(trecase, gt22[, .(gene_id, start,end,chrom)])

trecase <- add.null(trecase)

## remove bad runs

trecase <- trecase[Rhat < 1.1,]

###########################################
##' Open and format trecase with hidden-GT
###########################################

## read file with old output with tags matching GT and hidden-GT
## match.tags <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SingleCounting/results/refbias.gt.rna.p054.txt")

match.tags <- fread(snakemake@input[['gt_rna_sum']])

match.tags <- match.tags[,.(Gene_id, tag.rna,tag.gt, op.dir, tag.EAF.gt, tag.EAF.rna)]


merge.for <- function(i,j){
    dt <- comb.files(i, paste0("^", j,"\\.ENSG[0-9]+.*stan.summary.txt"))
    dt <- add.null(dt)
    ##remove bad runs
    dt <- dt[Rhat <1.1,]
    dt <- gene.d(dt, gt22[, .(gene_id, start,end,chrom)])
    return(dt)
    }

match.nogt <- function(i,j)  {
    dt <- merge.for(i,j)
    dt <- merge(match.tags, dt, by.x=c("Gene_id","tag.rna","tag.EAF.rna" ), by.y=c("Gene_id", "tag", "tag.EAF"))   
    dt[op.dir=="yes", log2_aFC_mean:= -log2_aFC_mean]
    return(dt)
}


## nogt <- match.nogt('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/RNA', 'rbias')

nogt <- match.nogt(snakemake@params[['nogt_mix']], 'rbias')

## all.nogt <- merge.for('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/RNA', 'rbias')

all.nogt <- merge.for(snakemake@params[['nogt_mix']], 'rbias')

## work with 100,000 kb

d <- 10^5

#############################
##' All associations by PIP #
#############################


## compare gt with no gt with null.99 and 95 within d

plots <- lapply(c("null.95", "null.99"), function(i) {
    merge.plot(trecase[abs(gene.dist)<=d,], nogt, null=i,
                                                       byx=c("Gene_id", "tag", "tag.EAF") ,
                                                       byy=c("Gene_id", "tag.gt", "tag.EAF.gt"),
                                                       suffixes=c(".gt", ".nogt"),
                                                       colsig=c("obs_GT", "hidden_GT"),
                                                       siglevel=gsub("null.", "",i),
                                                       xl="eQTL effect observed-GT",
                                                       yl="eQTL effect hidden-GT",
                                                       title=paste0("PIP = ", gsub("null.","", i), "%"),
               size.tab=8, xmin=-1,  xmax=-.5, ymin=.2, ymax=.4)
    } )
#+ fig.width= 7, fig.height=5.5
plot_grid(plotlist=plots, nrow=2)

###################################
##' Associations by gene distance #
###################################

l <- list(GT=trecase ,RNA=all.nogt)

tabs <- lapply(names(l), function(i) {
    lapply(c(95, 99), function(j){
        n <- paste0("null.",j)
        print(paste("Total associations with ", i ,"PIP = ",j))
        print(l[[i]][,.N, get(n)])
        print(paste("Associations within",d/1000, "KB from gene with", i, "and PIP =",j))
        print(l[[i]][abs(gene.dist)<=d,.N, get(n)])
        return(l[[i]][abs(gene.dist)<=d,.N, get(n)])
    })})




##################################################################################
## Look at individual genes with gt and rna only ##
###################################################################################



## to use previous functions need to add ".ngt" data, I will just duplicate rna for simplicity

rna.ngt <- merge(nogt,  nogt, by=c("Gene_id", "tag.rna", "tag.gt", "op.dir", "tag.EAF.rna","tag.EAF.gt", "gene.dist"), suffixes=c(".ngt", ".rna"))

rna.ngt[,tag:=tag.rna]

rna.ngt100 <- rna.ngt[abs(gene.dist)<=d,]


## merge trecase with nogt by level of significance

gt.rna <-  merge(trecase, nogt, by.x=c("Gene_id", "tag", "tag.EAF"), by.y=c("Gene_id", "tag.gt", "tag.EAF.gt" ), suffixes=c(".gt", ".rna"))
    
gt.rna.100 <- gt.rna[abs(gene.dist.gt) <= d,]
gt.rna.100 <- add.name(gt.rna.100)


all.100 <- merge(gt.rna.100, rna.ngt[,c("Gene_id", "tag.rna", grep(".ngt$", names(rna.ngt), value=T)), with=F], by.x=c("Gene_id", "tag.rna"), by.y=c("Gene_id", "tag.rna"))

## add col tag.ngt

all.100[ , tag.ngt:=tag.rna][, tag.gt:=tag]

plots.95 <- lapply(unique(all.100$Gene_id), function(i) gene.plot2b(x=trecase[abs(gene.dist) <= d,],
                   
                   null.x="null.95",
                   y=rna.ngt100,
                  
                   null.y="null.95",
                   
                   z=all.100,
                   gene=i))

plots.99 <- lapply(unique(all.100$Gene_id), function(i) gene.plot2b(x=trecase[abs(gene.dist) <= d,],
                   
                   null.x="null.99",
                   y=rna.ngt100,
                  
                   null.y="null.99",
                   
                   z=all.100,
                   gene=i))


## plot good examples

unique(all.100[Gene_id ==unique(all.100$Gene_id)[11], Gene_name])

#+ fig.width= 7, fig.height=5.5
plots.99[[11]]

unique(all.100[Gene_id ==unique(all.100$Gene_id)[76], Gene_name])
#+ fig.width= 7, fig.height=5.5
plots.99[[76]]

####################################################################################################################
## External validity/calibration: GTEX, MuTHER and GEUVADIS 
####################################################################################################################

## Use Gtex-ebv as gold standard and compare Deseq, Rasqual and Btrecae-GT

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
ebv <- ebv[Gene_id %in% unique(trecase$Gene_id),]
ebv[, SNP:=gsub("^22_|_b37$", "", variant_id)][, SNP:=gsub("_", ":",SNP)]

ebv.sig <- ebv.sig[variant_id %in% ebv$variant_id,][, null:="no"]
ebv <- merge(ebv,ebv.sig[,.(gene_id, variant_id, null)], by=c("gene_id", "variant_id"), all.x=T)
ebv[is.na(null), null:="yes"]

## add FDR by R for comparison

ebv[,p.adj:= p.adjust(pval_nominal,method = "BH")]
ebv[,null.fdr5:="yes"][p.adj<=0.05, null.fdr5:="no"]

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

## dseq <- merge(dseq, EAF, by.x="tag", by.y="snp")


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

## select associations run in dseq, as no tagging was done in rasqual

rasqual <-  merge(dseq[,.(Gene_id, tag)], rasqual, by.x=c("Gene_id", "tag"), by.y=c("gene_id", "rs_id"))

rasqual[, p_adjust:= p.adjust(p, method="BH")]

## transform Fold change to log2aFC
rasqual[, log2_aFC:=log2(Fold_change)]

## add various fdr cut-offs to rasqual

for (i in sort(c(10^(seq(-4,-1,1)), 0.05, 0.5))){
    rasqual[,eval(paste0("null.fdr",i*100)):= "yes"][p_adjust<=i, paste0("null.fdr",i*100):="no"]
}

#############################################################
## Select same associations in rasqual, btrecase-GT and dseq

assoc <- Reduce(fintersect, list(rasqual[,.(Gene_id,tag)], dseq[,.(Gene_id, tag)], trecase[, .(Gene_id,tag)]))

common.all <- lapply(list(Rasqual=rasqual, Btrecase=trecase, Deseq=dseq), function(i) merge(i, assoc, by=names(assoc)) )

## Convert rasqual and dseq to long format by FDR and btrecase using null using common associations only to ease comparison with Gtex-ebv

## make rasqual long format
rasq.l <- reshape(common.all$Rasqual[,c("Gene_id", "tag","log2_aFC", grep("null.fdr",names(common.all$Rasqual), value=T)), with=F],
                  direction="long",
                  varying=list( grep("null.fdr",names(common.all$Rasqual), value=T)),
                  v.names="null.Fdr",
                  times= as.numeric(gsub("null.fdr", "", grep("null.fdr",names(common.all$Rasqual), value=T)))/100,
                  timevar="Fdr")


## DEseq at various FDR to Gtex EBV cell lines at 5% FDR

dseq.l <- rbindlist(lapply(c(0.25, 0.1, 0.05, 0.01, 0.001), function(i) {
    tmp <- common.all$Deseq[, null.fdr:="yes"][p.adj<=i, null.fdr:="no"]
    tmp[, Fdr:= i]
    tmp <- tmp[, .(Gene_id, tag, log2FoldChange,log2_aFC, p.adj, null.fdr, Fdr, gene.dist)]
    return(tmp)
    }))

btrecase.gt.l <- rbindlist(lapply(c(95, 99), function(i) {
    dt <- rej.recode(a=0,b=common.all$Btrecase, c=i/100)
    dt[,PEP:=1-post.out]
    null <- paste0( "null.", i)
    tmp <- dt[, c("Gene_id", "tag", "log2_aFC_mean", "PEP", "null.rej", null), with=F]
    setnames(tmp, eval(null), "null")
    tmp[, PIP:= i/100]
    
}))

## extend btrecase.gt.l to PIP=0.8 and 0.9, null.rej corresponds to normal approx., copy to null to merge with btrecase.gt.l

post.dt <- rbindlist(lapply(c(0.8, 0.9), function(i) {
    dt <- rej.recode(a=0,b=common.all$Btrecase, c=i)
    dt[,PEP:=1-post.out]
    dt[,PIP:=i]
    dt[, null:=null.rej]
    dt <- dt[,.(Gene_id, tag, log2_aFC_mean,PEP,null.rej,null, PIP)]
}))

btrecase.gt.l <- rbind(btrecase.gt.l, post.dt)   

### Merge with gtex-ebv

rasq.ebv <- merge(rasq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
rasq.ebv <- add.signif(rasq.ebv, x1="null.Fdr", x2="null", col=c("Rasqual", "Gtex-ebv"))

dseq.ebv <- merge(dseq.l, ebv, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"))
dseq.ebv <-  add.signif(dseq.ebv, x1="null.fdr", x2="null", col=c("DEseq","Gtex-ebv") )

btrecase.gt.ebv <- merge(btrecase.gt.l, ebv,  by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "SNP"), suffixes=c(".btrecaseGT", ".ebv"))
btrecase.gt.ebv <- add.signif(btrecase.gt.ebv, x1="null.btrecaseGT", x2="null.ebv", col=c("Btrecase", "Gtex-ebv"))

## calculate FDR for trecase (GT) using normal approximation

btrecase.fdr <- rbindlist(lapply(c(0.99, 0.95, .9,.8), function(j){
    dt <- rej.recode(0,trecase,c=j)
    dt[,PEP:=1-post.out]
    rej <- nrow(dt[null.rej=="no",])
    false.pos <- dt[null.rej=="no",sum(PEP)]
    return(data.table(post.level=j, total.rej=rej, total.false.pos=false.pos, FDR=false.pos/rej))
    }))

btrecase.fdr


## calculate FDR using only associations shared with RASQUAL, Deseq2 and Gtex, very similar to FDR with all associations

fdr.shared <- rbindlist(lapply(unique(btrecase.gt.l$PIP), function(i) {
    dt <- btrecase.gt.l[PIP==i & null.rej=="no",]
    tmp <- data.table(post.level=i, total.rej=nrow(dt), total.fp=dt[,sum(PEP)], FDR=dt[,sum(PEP)]/nrow(dt))
    return(tmp)
}))

fdr.shared <- fdr.shared[order(FDR),]

###################################################
#### Make tables/plot  by number of associations
###################################################

all.ebv <- list(Rasqual=rasq.ebv, DEseq=dseq.ebv,Btrecase=btrecase.gt.ebv)

tabs <- mapply(function(a,b,c) {
    dt <- tab2bplot(dt=a,var=c("Signif", b))
    tables <- lapply(sort(as.numeric(unique(dt[[b]]))), function(i) dt[get(b) ==i,])
    rbindlist(lapply(tables, function(i) format.tab(i, c)))

},
a=all.ebv,
b=list("Fdr", "Fdr",  "PIP"),
MoreArgs=list(c="Gtex-ebv"), SIMPLIFY=F)

tabs.dt <- rbindlist(tabs, idcol="Method", fill=T)
tabs.dt[, Significant.level:=as.numeric(Fdr)][!is.na(PIP),Significant.level:=as.numeric(PIP)]
tabs.dt[, Discoveries:=Rasqual][!is.na(DEseq), Discoveries:=DEseq][!is.na(Btrecase), Discoveries:=Btrecase]

scaleFUN <- function(x) sprintf("%.2f", x)  ## for making consistent scales across plots

lab <- c("BaseQTL", "DESeq2", "RASQUAL")

## geom_path to connect line by FDR (order in dataset)
p.a <- ggplot(tabs.dt, aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method, label=Discoveries)) +
    geom_path(linetype="dotted") +
    geom_point() +
    geom_text_repel() +
    xlab("Positive predicted value") +
    ylab("Sensitivity") +
    ggtitle('"True" associations: 112') +
##     labs(title="External validity  by significant level (FDR or CI)",
##          subtitle="True discoveries in Gtex-EBV: 112",
##          caption="The number of total discoveries in GEUVADIS by each method is indicated\n
## FDR: Deseq: 0.25-5e-3, Rasqual: 0.5-1e-4, Btrecase CI: 0.8-0.99") +
##     theme(
##         plot.title = element_text(hjust = 0.5),
##         plot.subtitle = element_text(hjust = 0.5),
##         plot.caption = element_text(hjust = 0, lineheight = 0.5)            
    ##     ) +
    ## change label in legend
    scale_colour_manual(values=c("#009E73", "#CC79A7", "grey42"), labels=lab) +
    scale_shape_discrete(labels=lab) +
    scale_y_continuous(labels=scaleFUN) 
    
    

###############################
#### Make tables/plot by genes
###############################


sig.genes <- mapply(function(a,b,c,d, f, j, k){   
    rbindlist(lapply(sort(as.numeric(unique(a[[b]]))), function(i) {
         r <- unique(a[get(b) == i & get(j) =="no", Gene_id])
         e <- unique(a[get(b) == i & get(k) =="no", Gene_id])
         dt <- data.table(Signif=c(d, f, c), sig.level=i, SNPs=c(length(setdiff(r,e)), length(intersect(r,e)), length(setdiff(e,r))))
         return(format.tab(dt,c))
         }))},
    a=all.ebv,
    b=list("Fdr", "Fdr",  "PIP"),
    j=list("null.Fdr", "null.fdr", "null.btrecaseGT"),
    k=list("null", "null", "null.ebv"),
    d=names(all.ebv),
    MoreArgs=list(c="Gtex-ebv", f="Both"), SIMPLIFY=F)

tabs.g <- rbindlist(sig.genes, fill=T, idcol="Method")

tabs.g[, Discoveries:=Rasqual][!is.na(DEseq), Discoveries:=DEseq][!is.na(Btrecase), Discoveries:=Btrecase]

setorder(tabs.g, Method,  sig.level)


## geom_path to connect line by FDR (order in dataset)
p.g <- ggplot(tabs.g, aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method, label=Discoveries)) +
    geom_path(linetype="dotted") +
    geom_point() +
    geom_text_repel() +
    xlab("Positive predicted value") +
    ylab("Sensitivity") +
    ggtitle('"True" eGenes: 26') +
    scale_colour_manual(values=c("#009E73", "#CC79A7", "grey42"), labels=lab) +
    scale_shape_discrete(labels=lab) +
    scale_y_continuous(labels=scaleFUN) 
    
## get legend to share between the 2 plots
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

myleg <- g_legend(p.a)

p.ebv <- plot_grid(p.a + theme(legend.position="none"),
          NULL,
          p.g + theme(legend.position="non"),
          myleg, ncol=4,
          rel_widths=c(4.5,0.3, 4.5,1.5),
          labels=c('a', '', 'b', ''),
          label_x=c(0.1,0, 0.12,0)
          )

##########################################################################################
## Comparing associations to GEUVADIS bigger sample size (only sig associations available)
########################################################################################


## geu.eur <- fread("/mrc-bsu/scratch/ev250/EGEUV1/array_express_eqtl/EUR373.gene.cis.FDR5.all.rs137.txt")
geu.eur <- fread(snakemake@input[['geu_eur']])
geu.eur <- geu.eur[CHR_SNP==22 & CHR_GENE==22,][,Gene_id:=gsub("\\..*","", GENE_ID)]

## get ref and alt allele for each SNP_ID

snpmart <- useEnsembl(biomart="snp", dataset="hsapiens_snp", GRCh="37")

snp <- data.table(getBM(attributes=c('refsnp_id', 'chrom_start', 'allele', 'allele_1', 'minor_allele'),
             filters = 'snp_filter',
             values = unique(geu.eur[['SNP_ID']]),
             mart=snpmart))


## allele_1 is ancestral allele, allele is a list of all alleles

## GEUVADIS analysis from Chris
## select  chrom 22 and format compatible with gt
## sig.qtl <- fread('/home/ev250/rds/rds-cew54-wallace-share/Data/GEUVADIS/Geuvadis_sig_eqtl')
sig.qtl <- fread(snakemake@input[['geu_chris']])
sig.qtl <- sig.qtl[CHROM==22,][,Gene_id:=gsub("\\..*","", Gene)]
sig.qtl[, id := paste(POS,REF,ALT, sep=":")][, "null":="no"]
sig.qtl[, aFC:=2*Beta]


## Assess whether significant associations are in GEUVADIS bigger dataset with same direction of effects
geu.btrecase <-   merge(trecase[null.99 =="no" & abs(gene.dist) <= d,],
          sig.qtl,
          by.x=c("Gene_id", "tag"),
          by.y=c("Gene_id", "id"))


common.assoc <- geu.btrecase[sign(log2_aFC_mean) == sign(Beta),.N,]
total.assoc <- trecase[null.99 =="no"  & abs(gene.dist) <= d,.N,]

## 90% of assoc in GT are in GEU 

mapply(function(a,b) round(a*100/b),
       a=common.assoc,
       b=total.assoc)

## By genes

common.genes <- length(unique(geu.btrecase$Gene_id))
total.genes <-  length(unique(trecase[null.99 =="no" & abs(gene.dist) <= d ,Gene_id]))

mapply(function(a,b) round(a*100/b),
       a=common.genes,
       b=total.genes)



##########################################################################
## Compare all methods to GEU bigger dataset, selecting same associations
##########################################################################


## start from rasq.l, dseq.l and btrecase.gt.l (same gene-SNP associations)

common.all.l <- list(Rasqual=rasq.l, DEseq=dseq.l, Btrecase=btrecase.gt.l)


## get number of sig associations by FDR
tot.sig.assoc <- rbindlist(mapply(function(a,b,c) {
    dt <- a[get(c) == "no",.N, get(b)]
    setnames(dt, "get", b)
},
                        a=common.all.l,
                        b=c("Fdr", "Fdr", "PIP"),
                        c=c("null.Fdr", "null.fdr", "null"),
                        SIMPLIFY=F), idcol="Method", fill=T)
                        
## merge with sig.qtl (bigger sample size)
all.geu <-  lapply(common.all.l , function(i) merge(i, sig.qtl, by.x=c("Gene_id", "tag"), by.y=c("Gene_id", "id")))

## Count associations in same direction only, only interested in sig associations

sig.geu <- rbindlist(mapply(function(a, b, c, d) {
    dt <- a[sign(get(b)) == sign(get(c)), .N, d ]
    setnames(dt, names(dt)[2], "null")
    dt <- dt[null == "no",]
    return(dt)
},
a=all.geu,
b=c("log2_aFC", "log2_aFC", "log2_aFC_mean"),
d=list(c( "Fdr", "null.Fdr"),
       c( "Fdr", "null.fdr"),
       c( "PIP", "null.x")),
MoreArgs=list(c="Beta"),
SIMPLIFY=F), fill=T, idcol="Method")


## Merge total with sig.geu and format, make plot

comp.geu.assoc <- help.geu(tot.sig.assoc, sig.geu,  btrecase.fdr, N=all.geu$DEseq[,.N, Fdr][1,N] )

plots.a <- ggplot(comp.geu.assoc, aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method, label=Discoveries)) +
    geom_path(linetype="dotted") +
    geom_point() +
    geom_text_repel() +
    xlab("Positive predicted value") +
    ylab("Sensitivity") +
    ggtitle(paste0('"True" associations: ', unique(comp.geu.assoc$sig.geu))) +
    ## change label in legend
    scale_colour_manual(values=c("#009E73", "#CC79A7", "grey42"), labels=lab) +
    scale_shape_discrete(labels=lab) +
    scale_y_continuous(labels=scaleFUN)  +
    xlim(0.15,1)
    


## By genes
tot.sig.genes <- rbindlist(mapply(function(a,b,c) {
    dt <- a[get(c) == "no", unique(Gene_id), get(b)][,.N, get]
    setnames(dt, "get", b)
},
a=common.all.l,
b=c("Fdr", "Fdr", "PIP"),
c=c("null.Fdr", "null.fdr", "null"), SIMPLIFY=F),
idcol="Method", fill=T)

sig.geu.g <- rbindlist(mapply(function(a, b, c, d) {
    dt <- a[sign(get(b)) == sign(get(c)), unique(Gene_id), d ][ ,.N,d]
    setnames(dt, names(dt)[2], "null")
    dt <- dt[null == "no",]
    return(dt)
},
a=all.geu,
b=c("log2_aFC", "log2_aFC", "log2_aFC_mean"),
d=list(c( "Fdr", "null.Fdr"),
       c( "Fdr", "null.fdr"),
       c( "PIP", "null.x")),
MoreArgs=list(c="Beta"),
SIMPLIFY=F), fill=T, idcol="Method")


comp.geu.g <- help.geu(tot.sig.genes, sig.geu.g, btrecase.fdr, N=all.geu$DEseq[,unique(Gene_id), Fdr][,.N,Fdr][1,N] )

plots.g <- ggplot(comp.geu.g, aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method, label=Discoveries)) +
    geom_path(linetype="dotted") +
    geom_point() +
    geom_text_repel() +
    xlab("Positive predicted value") +
    ylab("Sensitivity") +
    ggtitle(paste0('"True" eGenes: ', unique(comp.geu.g$sig.geu))) +
    ## change label in legend
    scale_colour_manual(values=c("#009E73", "#CC79A7", "grey42"), labels=lab) +
    scale_shape_discrete(labels=lab) +
    scale_y_continuous(labels=scaleFUN) +
    xlim(0.15, 1)
    

## Combine p.a, p.g, plots.a and plots.g

ex.all <- plot_grid(plots.a + theme(legend.position="none"),
          NULL,
          plots.g + theme(legend.position="non"),
          myleg,
          p.a + theme(legend.position="none"),
          NULL,
          p.g + theme(legend.position="non"),
          myleg, ncol=4,
          rel_widths=c(4.5,0.3, 4.5,1.5),
          labels=c('a', '', 'b', '', 'c', '', 'd', ''),
          label_x=rep(c(0.1,0, 0.12,0),2)
          )

##########################################################################
## Compare all methods to MuTHER, selecting same associations
##########################################################################


## muther <- fread('/mrc-bsu/scratch/ev250/MuTHER/MuTHER_cis_results_chr22_B37.txt')
muther <- fread(snakemake@input[['muther37']])

## remove Gene NA

muther <- muther[!is.na(Gene),]

## I have rs_id and A1 allele, need ref/alt to match with my tags. To reduce biomart time I select SNPs with same positions as my tags (assoc). I add POS, ref and alt cols to assoc.

assoc[, POS:=as.numeric(gsub(":.*","",tag))][, ref:=sub(".*[0-9]:([A-Z]):[A-Z]", "\\1", tag)][, alt:=sub(".*:", "", tag)]
length(unique(assoc$tag))

## select unique tags

utags <- unique(assoc[, .(tag, POS,ref,alt)])

## select SNPs in muther in same positions as utags

muther <- muther[ POS %in% utags$POS,]


## get ref and alt allele for SNPs.m

snp.m <- data.table(getBM(attributes=c('refsnp_id', 'chrom_start', 'allele', 'allele_1', 'minor_allele'),
             filters = 'snp_filter',
             values =  unique(muther$SNP),
             mart=snpmart))

## some rsids map to more than 1 location, select those in same position as muther

snp.m <- merge(snp.m , unique(muther[,.(POS, SNP, A1)]), by.x=c("chrom_start", "refsnp_id"), by.y=c("POS", "SNP"))

## add snp.m to muther

muther <- merge(muther, snp.m, by.x=c("POS", "SNP", "A1"), by.y=c("chrom_start", "refsnp_id", "A1"))


## to select same Gene-SNP associations I need gene name to match to assoc

assoc <- add.name(assoc)

## select same gene-snp associations in muther and assoc, snp by position

muther <- merge(muther, assoc, by.x=c("Gene", "POS"), by.y=c("Gene_name", "POS"))


## Not further follow up as I can work with Chris data for 462 ind from GEUVADIS
