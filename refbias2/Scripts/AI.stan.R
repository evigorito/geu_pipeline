#' ---
#' title: Compare stan output with new AI estimates
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule external_validity from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias2/Snakefile

library(data.table)
library(grid)
library(gridExtra)
source('/home/ev250/Cincinatti/Functions/various.R')
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

## dirs with stan summaries
oldir <-  "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/GT"
newdir <-  "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"
oldir <- snakemake@params[['old_dir']]
newdir <- snakemake@params[['new_dir']]

## AI estimates
oAI <- snakemake@input[['oldAI']]
nAI <- snakemake@input[['newAI']]

## Get stan summaries
maxRhat <- 1.1

dirs <- setNames(c(oldir, newdir), c("old", "new"))
baseq <- lapply(dirs, function(i) {
    dt <- comb.files(path=i,
                      pattern=paste0(paste0("^",
                                            rbias,
                                            "\\.ENSG[0-9]+.*stan.summary.txt")))
    dt[Rhat<maxRhat,]
})
names(baseq) <- names(dirs)

## merge and compare
s <- c(".old", ".new")
x <- paste0("null.99", s)
base.w <- Reduce(function(a,b) merge(a,b, by=c("Gene_id", "tag"), suffixes=s), baseq)
base.w <- add.signif(base.w, x[1], x[2], gsub("\\.", "",s))

cols <- c(None="#999999", old="yellow3", new="#0072B2", Both= "#D55E00")
tab <- tab2bplot(dt=base.w, colors=cols)
btrecase.plot(dt=base.w[Signif != "None",],
                        x1=paste0(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', "null.99"),s[1]),
                        x2= paste0(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', 'null.99'),s[2]),
                        xl='eQTL-effect (old)', yl='eQTL-effect (new)',
                        col=gsub("\\.", "",s),axis.title=12, axis.text=10,
                        legend.title=12, legend.text=10,
                        legend.symbol=4, point.size=3 ,
                        title="Old vs new estimates", title.size=12) +
    annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=-1.1, xmax=-0.3, ymin=0.35, ymax=0.45)


## Compare both datasets with GEUADIS big data analysis

## sig.qtl <- fread('/home/ev250/rds/rds-cew54-wallace-share/Data/GEUVADIS/Geuvadis_sig_eqtl')
sig.qtl <- fread(snakemake@input[['geu_chris']])
sig.qtl <- sig.qtl[CHROM==22,][,Gene_id:=gsub("\\..*","", Gene)]
sig.qtl[, id := paste(POS,REF,ALT, sep=":")][, "null":="no"]
sig.qtl[, aFC:=2*Beta]

baseqtl.fdr <- lapply(baseq, function(j){
    rbindlist(lapply(c(99,95,90,85, 80)/100, function(i){
    dt <- rej.recode(a=0,b=j, c=i)
    dt[,PEP:=1-post.out]
    ## get total rejections and false positives to calc Fdr
    rej <- nrow(dt[null.rej=="no",])
    false.pos <- dt[null.rej=="no",sum(PEP)]   
    dt[,Fdr:=false.pos/rej]
    dt[, null.Fdr:=null.rej]   
    dt <- dt[,.(Gene_id, tag, log2_aFC_mean, log2_aFC_sd, null.Fdr, Fdr)]
    }))
})

tot.sig.a <- rbindlist(lapply(baseqtl.fdr,
                              function(a) {
                                  dt <- a[null.Fdr == "no",.N, Fdr]
                                  ## check if any Fdr has 0 total associations
                                  w <- which(! unique(a$Fdr) %in% dt$Fdr)
                                  if(length(w)) dt <- rbind(data.table(Fdr=unique(a$Fdr)[w], N=0), dt)
                                  return(dt)
                                  }
                             ),
                       idcol="Method", fill=T)

all.geu <-  lapply(baseqtl.fdr , function(i) merge(i, sig.qtl,
                                                    by.x=c("Gene_id", "tag"),
                                                    by.y=c("Gene_id", "id")))

sig.count <-  rbindlist(lapply(all.geu, function(i) {
    dt <- i[sign(log2_aFC_mean) == sign(Beta), .N, c( "Fdr", "null.Fdr")]
    dt <- dt[null.Fdr == "no",]
    
    return(dt)
}) , fill=T, idcol="Method")

## Merge total with sig.geu and format, make plot

comp.geu.assoc <- merge(sig.count[, .(Method, Fdr, N)],
                        tot.sig.a, by=c("Method","Fdr"),
                        suffixes=c(".sig", ".total"), all.y=T)

## add total sig assoc in geu that have been tested by methods
comp.geu.assoc[, sig.geu:=unique(unlist(lapply(all.geu, function(i) i[, .N,Fdr]$N)))]

## calculate Power and PPV, PPV gives NA when 0/0, change to 0
comp.geu.assoc[, Sensitivity:=N.sig/sig.geu][, PPV:=N.sig/N.total][is.na(PPV), PPV:=0]

fdr.plot(comp.geu.assoc,
          x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="associations", col="sig.geu",
         path="no",
         label="no")
