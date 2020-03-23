#' ---
#' title: Benchmarking eQTL methods v2
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule external_validity from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/refbias/Snakefile

library(data.table)
source('/home/ev250/Cincinatti/Functions/various.R')
source("/home/ev250/Bayesian_inf/trecase/Functions/out.trecase.many.genes.R")

##' Open methods results: Dseq2, RASQUAL, lm and BaseQTL

## params to get relevant output: BaseQTL reference panel bias correctiona and observed genotypes
rbias <- snakemake@params[['rbias']]
## rbias="rbias"


## FDR range for frequentist tests

r.fdr <- sort(c(10^(seq(-7,-1,1)), 0.05, 0.015, 0.15, 0.5, 0.8, 1))

## add BH corrected pvalue and significance by fdr cut-off
## @param dt data table with output
## @param pv name of column with raw p-value
## @param r.fdr vector with a range of fdr cut-offs
 
cut.fdr <- function(dt, pv="pvalue", r.fdr){
    dt[, p_adjust:=p.adjust(get(pv) ,method = "BH")]
    for (i in r.fdr){
        dt[,eval(paste0("null.fdr",i*100)):= "yes"][p_adjust<=i, paste0("null.fdr",i*100):="no"]
    }
    return(dt)
}


###########
## Deseq2 #
##########

## Open and format DEseq2 output (run by Wei-Yu), various shrinkage options

## files

dseq.f <- snakemake@input[['dseq']]
## dseq.f <- list.files(path="/mrc-bsu/scratch/wyl37/ElenaData/RunNBmodelshrinkage", pattern="RunNBmodelshrinkagebatch", full.names=T)

dseq <- rbindlist(lapply(dseq.f,  fread))

## format dseq

dseq.fmt <- function(i) {    
    i[, SNP:=NULL]
    ## remove pval NA
    i <- i[!is.na(pvalue),]
    ## add log2_aFC to dseq
    i[, log2_aFC_mean:=log2FoldChange*2]
    i[, log2_aFC_se_mean:= 2*lfcSE]
    i[, pvalue:=2*pnorm(log2FoldChange/lfcSE, lower.tail=FALSE)]  ## temporary until Wei-Yu fixes pvalues
    ## add BH correction to dseq
    i <- cut.fdr(dt=i, r.fdr=r.fdr)
}

dseq.l <- lapply(unique(dseq$shrink.method), function(i) dseq.fmt(dseq[shrink.method==i,]))
names(dseq.l) <- unique(dseq$shrink.method)

#############
## RASQUAL ##
#############

## rasqual files and header

rasq86 <- list.files(snakemake@params[['rasqual_86']], "ENSG[0-9]+.*txt", full.names=T)
rasqPer <- list.files(snakemake@params[['rasqual_per']], "ENSG[0-9]+.*txt", full.names=T)
rasq.size.dir  <- snakemake@params[['rasqual_size']]
rasq.size <- lapply(rasq.size.dir, function(i) list.files(i, "ENSG[0-9]+.*txt", full.names=T))

rasq.header <- snakemake@input[['rasqual_header']]


## rasqPer <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/outASE_per/cis1_10_5/","ENSG[0-9]+.*txt", full.names=T)
## rasq86 <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/outASE/cis1_10_5/","ENSG[0-9]+.*txt", full.names=T)
## rasq.header <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/output/rasqual.header.txt"

## rasq.size.dir  <- sapply(c(5,25), function(i) paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/out",
## i, "ase/cis1_10_5"))
## rasq.size <- lapply(rasq.size.dir, function(i) list.files(i, "ENSG[0-9]+.*txt", full.names=T))

## read and format rasqual output, adjust pvalue by BH and add significant at different Fdr cut-offs

## format rasqual output
## @param x character vector with rasqual output by gene
## @param y data table with Gene_id and tags to select from rasqual
## @param fdr whether to adjust fdr using BH
 
rasq <- function(x,y, fdr="BH"){ ## format rasqual input

    rasqual <- rbindlist(lapply(x, format_rasqual, top.hits="no", header=rasq.header))
    

    ## rasqual didnt run some genes "Estimated computational time is too long...",
    ## remove from output
    ##print(nrow(rasqual[rs_id == "SKIPPED",]))
    rasqual <- rasqual[rs_id != "SKIPPED",]
    ## make sure rs_id is POS:REF:ALT
    if(all(unique(rasqual$rs_id) == ".")) rasqual[, rs_id := paste(SNP_pos, Ref, Alt, sep=":")]
    ## select same tags run in dseq, all elements of dseq have same assoc
    rasqual <-  merge(unique(y[,.(Gene_id, tag)]), rasqual, by.x=c("Gene_id", "tag"),
                      by.y=c("gene_id", "rs_id"))
    if(fdr=="BH"){
        rasqual <- cut.fdr(dt=rasqual, pv="p", r.fdr=r.fdr)
    } 
    rasqual[, log2_aFC_mean:=log2(Fold_change)]
    return(rasqual)
}

rasq1 <- rasq(rasq86, dseq.l[[1]])
rasq2 <- rasq(rasqPer, dseq.l[[1]])

rasq_comp <- merge(rasq1, rasq2, by=c("Gene_id", "tag"), suffixes= c(".ok", ".perm"))

ggplot(rasq_comp, aes(log2_aFC_mean.ok, log2_aFC_mean.perm)) + geom_point()

r.comp.l <- rbindlist(list(Real=rasq1, Permutated=rasq2), idcol="Test")

ggplot(r.comp.l, aes(p)) +
    geom_histogram(color="black",fill="white") +
    facet_grid(Test~.) +
    ggtitle("RASQUAL with input permutation") +
    xlab("pvalue")



## format rasqual runs
rasqual <- lapply(c(list(rasq86), rasq.size), rasq, y=dseq.l[[1]])

## name by number of inds
names(rasqual) <- c("rasq86",
                    paste0("rasq",
                           sapply(rasq.size.dir,
                                  function(i) gsub(".*out([0-9]{1,2}).*",
                                                   "\\1", i))))
#### look at FDR using permutation p-value
## Plot FDR=#per-tests | pvalue-perm < alpha / #tests |pvalue<alpha

## x data table with rasqual test and permutation test
## y range of alpha to calculate fdr based on permutations
## z suffixes for tests without and with permutated inputs
fdr.per <- function(x,y,z=c(".ok", ".perm")){
    tmp <- lapply(y, function(a) {
        per <- x[get(paste0("p", z[2])) <a,.N]
        ok <- x[get(paste0("p", z[1])) <a,.N]
        if(ok==0) return(0)
        return(per/ok)
    })
    return(unlist(tmp))
}

alpha <- exp(seq(from=log(10^-15), to=log(0.01), length.out=500))
fdr.a <- fdr.per(rasq_comp, alpha)   
plot(x=alpha, y=fdr.a)

dt <- data.table(fdr=fdr.a, a=alpha)

## Select "a" to get fdr approx c( 0.001, 0.010, 0.015, 0.050)

f.approx <- c( 0.001, 0.010, 0.015, 0.050, 0.1,0.15)
## get closest alpha to f.approx

a.approx <- rbindlist(lapply(f.approx, function(i) {
    dt[, paste0("fdr.", i) := abs(fdr-i)]
    tmp <- dt[get(paste0("fdr.", i)) == min(get(paste0("fdr.", i))), c("fdr", paste0("fdr.", i), "a"), with=F]
    tmp[, f.approx:=i]
    tmp <- tmp[, paste0("fdr.", i):=NULL]
    ## when no false positives (no assoc in perm below p value cut-off) select the highest a (more assoc in rasqual)
    if( unique(tmp$fdr) == 0 ){
        tmp <- tmp[a==max(a),]
    }
    
    return(tmp)
}))


## add permutated corrected pvalue and significance by fdr cut-off
## @param dt data table with output not adjusted
## @param pv name of column with raw p-value
## @param dt2 data table with fdr calcualted by permutations, cols fdr, a, f.approx
 
per.fdr <- function(dt, pv="pvalue", dt2){
    for (i in 1:nrow(dt2)){
        fdr <- dt2[i, f.approx]
        dt[,eval(paste0("null.fdr",fdr*100)):= "yes"][get(pv) < dt2[i,a], paste0("null.fdr",fdr*100):="no"]
    }
    return(dt)
}

rasq.perm <- per.fdr(dt=rasq(rasq86, dseq.l[[1]], fdr="no"), pv="p", dt2=a.approx)

## add to rasqual
rasqual[['rasqPerm']] <- rasq.perm


#########
## lm ###
#########

## lm dirs
lm86.dir <- snakemake@params[['lm_86']]
lm.size <- snakemake@params[['lm_size']]
## lm86.dir <- c("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/lm/log_counts",
##               "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/lm/Inds86_gclibsize")

## lm.size <- c("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/lm/Inds5",
##              "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/lm/Inds25")

## get files for each run
lm.f <- lapply(c(lm86.dir, lm.size), function(i) list.files(i,
                                                            "ENSG[0-9]+.*summary.txt",
                                                            full.names=T))
names(lm.f) <- paste0("lm.", sapply(c(lm86.dir, lm.size), basename))

## get lm output, add BH p_adjusted and Fdr cut-off

lm <- lapply(lm.f, function(i) {
    dt <-rbindlist(lapply(i,  fread))
    ## if standar error isn't in output, add it
    if(!"log2.aFC_se" %in% names(dt)){
        dt[, log2_aFC_se_mean:= abs(log2.aFC_97.5 - log2.aFC_2.5)/3.92]
    } else{ ## change name for consistency with other datasets
        setnames(dt, "log2.aFC_se", "log2_aFC_se_mean")
    }
    ## change SNP for tag and "log2.aFC_mean" for "log2_aFC_mean"
    ## for consistency with other datasets
    if(any(names(dt) == "SNP")) setnames(dt, "SNP", "tag")
    setnames(dt, "log2.aFC_mean", "log2_aFC_mean")
    ## adjust pvalue
    dt <- cut.fdr(dt, pv="log2.pvalue" , r.fdr)
    return(dt)
})



################
## BaseQTL #####
################

maxRhat <- 1.1

## btrecase_dir  <- paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/", c("GT", "GT25"))
btrecase_dir <- snakemake@params[['btrecase']]

## open files and apply QC
baseqtl <-lapply(btrecase_dir, function(i){
    dt <- comb.files(path=i,
                      pattern=paste0(paste0("^",
                                            rbias,
                                            "\\.ENSG[0-9]+.*stan.summary.txt")))
    dt <- dt[Rhat < maxRhat,]
    ## change SNP for tag for consistency with other datasets
    if(any(names(dt) == "SNP")) setnames(dt, "SNP", "tag")
    return(dt)
    })
names(baseqtl) <- unlist(lapply(basename(btrecase_dir), function(i) gsub("GT", "BaseQTL",i)))

##base2 <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT/rbias.ENSG00000159958.stan.summary.txt")

## base3 <- merge(baseqtl, base2, by=c("Gene_id", "tag"))

## for various rejection rules calculate PEP of opposite sign as
## mean log2aFC, use normal approximation of posterior and calculate Fdr
## Calculate fdr for stan summary
## @param baseqtl data table with stan sum 
base.fdr <- function(baseqtl) {
    rbindlist(lapply(c(99,95,90,85, 80)/100, function(i){
    dt <- rej.recode(a=0,b=baseqtl, c=i)
    dt[,PEP:=1-post.out]
    ## get total rejections and false positives to calc Fdr
    rej <- nrow(dt[null.rej=="no",])
    false.pos <- dt[null.rej=="no",sum(PEP)]   
    dt[,Fdr:=false.pos/rej]
    dt[, null.Fdr:=null.rej]   
    dt <- dt[,.(Gene_id, tag, log2_aFC_mean, log2_aFC_sd, null.Fdr, Fdr)]
    }))
}

baseqtl.fdr <- lapply(baseqtl, base.fdr)

##################
## save objects
###################

saveRDS(lm, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/lm.rds")
saveRDS(dseq.l, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/dseq.l.rds")
saveRDS(rasqual, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/rasqual.rds")
saveRDS(baseqtl,  "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/baseqtl.GT.gc.rds")
saveRDS(baseqtl.fdr,  "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/baseqtl.fdr.GT.gc.rds")

###############################################################
## Select same associations in rasqual, baseqtl, dseq and lm
################################################################

## common associations
assoc <- Reduce(function(a,b) {
    fintersect(a[,.(Gene_id,tag)], b[,.(Gene_id,tag)])
}, c(dseq.l, rasqual, lm, baseqtl))


common.all <- lapply(c(dseq.l, rasqual, lm, baseqtl),
                     function(i) merge(i, assoc, by=names(assoc)) )

####################################################
## Convert to long format by FDR to ease plotting ##
####################################################

## format long
## @param dt data table to make long
## @param n name of dt cols to keep in output
long <- function(dt, n){
    tmp <- dt[,n, with=F]
    tmp2 <- reshape(tmp,
            direction="long",
            varying=list( grep("null.fdr",names(tmp), value=T)),
            v.names="null.Fdr",
            times= as.numeric(gsub("null.fdr", "", grep("null.fdr",names(tmp), value=T)))/100,
            timevar="Fdr")
    tmp2[, id:=NULL]
    
}


dseq2.l <- lapply(common.all[names(dseq.l)],
                  function(i) long(i, c("Gene_id", "tag","log2_aFC_mean","p_adjust","log2_aFC_se_mean",
                                  grep("null.fdr",names(i), value=T))))

lm.l <- lapply(common.all[names(lm)],
               function(i) long(i,
                                c("Gene_id", "tag","log2_aFC_mean","p_adjust","log2_aFC_se_mean",
                                  grep("null.fdr",names(i), value=T))))

rasq.l <- mapply(function(i,j) {
    tmp <- c("Gene_id", "tag","log2_aFC_mean", j, grep("null.fdr",names(i), value=T)) ## names in data table
    long(i, tmp)
},
                i= common.all[names(rasqual)],
                j=c(rep("p_adjust",3), "p"),
SIMPLIFY=F)


## restrict to common assoc
baseqtl.l <- lapply(baseqtl.fdr, function(i) merge(i, assoc, by=names(assoc)))

###########################################################
## Comparing associations to GEUVADIS bigger sample size ##
###########################################################


## GEUVADIS analysis from Chris
## select  chrom 22 and format compatible with methods outputs
## sig.qtl <- fread('/home/ev250/rds/rds-cew54-wallace-share/Data/GEUVADIS/Geuvadis_sig_eqtl')
sig.qtl <- fread(snakemake@input[['geu_chris']])
sig.qtl <- sig.qtl[CHROM==22,][,Gene_id:=gsub("\\..*","", Gene)]
sig.qtl[, id := paste(POS,REF,ALT, sep=":")][, "null":="no"]
sig.qtl[, aFC:=2*Beta]


## gather datasets with common assoc in long form
common.all.l <- c(dseq2.l, rasq.l,  lm.l, baseqtl.l)

## save
saveRDS(common.all.l, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/common.all.l.rds")


tot.sig.a <- rbindlist(lapply(common.all.l,
                              function(a) {
                                  dt <- a[null.Fdr == "no",.N, Fdr]
                                  ## check if any Fdr has 0 total associations
                                  w <- which(! unique(a$Fdr) %in% dt$Fdr)
                                  if(length(w)) dt <- rbind(data.table(Fdr=unique(a$Fdr)[w], N=0), dt)
                                  return(dt)
                                  }
                             ),
                       idcol="Method", fill=T)

## merge with sig.qtl (bigger sample size)
all.geu <-  lapply(common.all.l , function(i) merge(i, sig.qtl,
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

## add Fdr.approx col to plot baseqtl in the same scale as freq methods
setkey(comp.geu.assoc, Method, Fdr)
comp.geu.assoc[, Fdr.approx:=Fdr][Method %in% c("BaseQTL","BaseQTL25"),
                                  Fdr.approx:=c( 0.001, 0.010, 0.050, 0.1, 0.15)][, Fdr.approx:=as.factor(Fdr.approx)]

## Save

saveRDS(comp.geu.assoc, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/comp.geu.assoc.rds")


## compare all methods with 86 inds with Fdr 0.001-0.01 (as before), exclude deseq2
assoc <- fdr.plot(dt=comp.geu.assoc[Fdr.approx %in% c( 0.001, 0.010, 0.050, 0.1, 0.15) & 
                                    Method %in% c( "BaseQTL","lm.Inds86_gclibsize", "rasqPerm"), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="associations", col="sig.geu",
         path="no",
         label="yes")


assoc

fdr.plot(dt=comp.geu.assoc[Method %in% c("BaseQTL","lm.Inds86_gclibsize", "rasq86"), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="associations", col="sig.geu",
         path="no",
         label="no")

## Compare rasq BH with rasq permutation adjusted
fdr.plot(dt=comp.geu.assoc[Method %in% c("lm.Inds86_gclibsize", "rasq86", "rasq.perm"), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="associations", col="sig.geu",
         path="no",
         label="no")

## rasq perm reduces discoveries for similar FDR
comp.geu.assoc[Method %in% c( "rasq86", "rasq.perm"),1:7, with=F ]

## compare methods with 25 individuals
fdr.plot(dt=comp.geu.assoc[Fdr.approx %in% c(1e-05, 1e-04 ,0.001, 0.010, 0.050, 0.1, 0.15) &
                           Method %in% c(names(lm)[4], "rasq25",  "BaseQTL25"), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="associations", col="sig.geu",
         path="no",
         label="no")


## compare lm with different covariates



## compare dseq

## fdr.plot(dt=comp.geu.assoc[Fdr.approx %in% c(0.001, 0.010, 0.050, 0.10) &
##                                        Method %in% c(names(dseq.l), "lm.Inds86_gclibsize" ), ],
##                   x="PPV",
##                   y="Sensitivity",
##                   lab.col="N.total",
##          type="associations", col="sig.geu",
##          path="no")

## look at egenes

tot.sig.g <- rbindlist(lapply(common.all.l,
                              function(a) {
                                  dt <- a[null.Fdr == "no", unique(Gene_id), Fdr][,.N, Fdr]
                                  ## check if any Fdr has 0 total associations
                                  w <- which(! unique(a$Fdr) %in% dt$Fdr)
                                  if(length(w)) dt <- rbind(data.table(Fdr=unique(a$Fdr)[w], N=0), dt)
                                  return(dt)
                              }),
                       idcol="Method", fill=T)

sig.geu.gene <- rbindlist(lapply(all.geu,
                                 function(a) {
                                     dt <- a[sign(log2_aFC_mean) == sign(Beta), unique(Gene_id), c("Fdr", "null.Fdr")][,.N, .(Fdr, null.Fdr)]
                                     dt <- dt[null.Fdr == "no",]
                                     return(dt)
                                 }),
                          fill=T, idcol="Method")

comp.geu.g <- merge(sig.geu.gene, tot.sig.g,  by=c("Method","Fdr"),
                    suffixes=c(".sig", ".total"), all.y=T)


## add total sig assoc in geu that have been tested by methods
comp.geu.g[, sig.geu:= unique(unlist(lapply(all.geu, function(i) i[,unique(Gene_id), Fdr][,.N,Fdr][1,N]))) ]

## calculate Power and PPV, PPV gives NA when 0/0, change to 0
comp.geu.g[, Sensitivity:=N.sig/sig.geu][, PPV:=N.sig/N.total][is.na(PPV), PPV:=0]

## add Fdr.approx col to plot baseqtl in the same scale as freq methods
setkey(comp.geu.g, Method, Fdr)


comp.geu.g[, Fdr.approx:=Fdr][Method %in% c("BaseQTL", "BaseQTL25"),
                                  Fdr.approx:=c( 0.001, 0.010, 0.050, 0.10, 0.15)][, Fdr.approx:=as.factor(Fdr.approx)]


## change label for method
comp.geu.g[Method =="GT", Method:="BaseQTL"][Method=="GT25", Method:="BaseQTL25"]

##save
saveRDS(comp.geu.g, "/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/comp.geu.g.rds")

## compare all methods with 86 inds
fdr.plot(dt=comp.geu.g[  Method %in%  c( "BaseQTL","lm.Inds86_gclibsize", "rasq86"), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
                  type="eGenes", col="sig.geu",
                  path="no")



genes <- fdr.plot(dt=comp.geu.g[ Fdr.approx %in% c( 0.001, 0.010, 0.050, 0.10, 0.15) & Method %in%  c("BaseQTL","lm.Inds86_gclibsize", "rasqPerm"), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="eGenes", col="sig.geu",
         path="no")

plot <- plot_grid(assoc, genes, nrow=1)
ggsave("/mrc-bsu/scratch/ev250/EGEUV1/quant/outputPaper/assoc_genes.png", plot,width=9.2, height=3.5)

## dif number of inds
print(names(lm)[4])
fdr.plot(dt=comp.geu.g[ Fdr.approx %in% c(1e-07,1e-06, 1e-05, 1e-04, 0.001, 0.010, 0.050, 0.10, 0.15,.5) &
                        Method %in%  c(names(rasqual)[c(3)],"BaseQTL25",  names(lm)[c(4)]), ],
                  x="PPV",
                  y="Sensitivity",
                  lab.col="N.total",
         type="eGenes", col="sig.geu",
         path="no")




#########################
## eGENES as in Rasqual##
#########################


## get genes by fdr and add col wether the gene is true positive or false positive
## @param l list in long format by Fdr with methods
## @param assoc data table with common associations across methods to define tested genes
## @param egenes data table with genes in gold standard

egenes.f <- function(l,assoc, egenes){
    ## get true positives: Genes in assoc and in egenes
    tp <- unique(assoc$Gene_id)[unique(assoc$Gene_id) %in% unique(egenes$Gene_id)]
    ## get true negative: Genes in assoc and not in egenes
    tn <- unique(assoc$Gene_id)[!unique(assoc$Gene_id) %in% unique(egenes$Gene_id)]

    rbindlist(mapply(function(a,b,c) {
        ## select method positives
        u <- unique(a[get(c)== "no" , c("Gene_id", b), with=F])
        ## get positive genes in gold standard and add aux col
        dt <- merge(u, data.table(Gene_id = unique(egenes$Gene_id),
                                  col=1:length(unique(egenes$Gene_id))), by="Gene_id", all.x=T)
        ## na values for col are false positives
        TP=dt[!is.na(col), .N, get(b)]
        names(TP) <- c(b, "TP")
        FP <- dt[is.na(col), .N, get(b)]
        names(FP) <- c(b, "FP")
        tmp <- merge(TP, FP, by=b, all=T)
        ## if tmp has NA convert to 0
        tmp[is.na(tmp)] <- 0
        ## add number of gold standard positives and negatives
        ## to calculate sensitivity and specificity
        tmp[,GEUpos:=length(tp)][,GEUneg:=length(tn)]
        tmp[, Sensitivity:=100*TP/GEUpos][, FPR:=100*FP/GEUneg]

        return(tmp)  
        
    },
    a=l,
    MoreArgs=list(b="Fdr",
                  c="null.Fdr"),
    SIMPLIFY=F),
    idcol="Method", fill=T)
}



## compare all methods with 86 inds (as before)
tpfp.geu <- egenes.f(common.all.l[c(#"DESeq2",
                                    "BaseQTL","lm.Inds86_gclibsize", "rasq86")], assoc, sig.qtl)
roc1 <- ggplot(tpfp.geu, aes(FPR, Sensitivity, color=Method))+ geom_line() + geom_point()

roc1

## compare all methods with 86 inds (as before)
tpfp.geu <- egenes.f(common.all.l[c(#"DESeq2",
                                    "BaseQTL","lm.Inds86_gclibsize", "rasq.perm")], assoc, sig.qtl)
roc.perm <- ggplot(tpfp.geu, aes(FPR, Sensitivity, color=Method))+ geom_line() + geom_point()

roc.perm


## compare lm 86 inds

lm.geu <- egenes.f(common.all.l[names(lm)[1:2]], assoc, sig.qtl)
roc2 <-  ggplot(lm.geu, aes(FPR, Sensitivity, color=Method))+ geom_line() + geom_point()
roc2


## compare all methods with 86 inds (as before)
tpfp.geu <- egenes.f(common.all.l[c(#"DESeq2",
                                    "BaseQTL25","lm.Inds25", "rasq25")], assoc, sig.qtl)
roc.25 <- ggplot(tpfp.geu, aes(FPR, Sensitivity, color=Method))+ geom_line() + geom_point()

roc.25

## compare lm/rasq by inds

## rasq.lm <- egenes.f(common.all.l[c(names(lm)[2:4], names(rasqual))], assoc, sig.qtl)
## roc3 <- mapply(function(a, b) ggplot(rasq.lm[Method==a | Method==b,],
##                                      aes(FPR, Sensitivity, color=Method)) +
##                               geom_point()+ geom_line(),
##                a=names(lm)[2:4],
##                b=names(rasqual),
##                SIMPLIFY=F)

#+ fig.width=5.5, fig.height=5.8
#plot_grid(plotlist=roc3, ncol=1)


## Add adjustement method for each method and save

## long.all <-rbindlist( mapply(function(a,b) {
##     a[, adjustment:=b]
##     ## merge log2_aFC_se_mean with log2_aFC_sd
##     if(any(names(a) == "log2_aFC_se_mean")){
##         setnames(a, "log2_aFC_se_mean", "log2_aFC_se/sd")
##     } else {
##         setnames(a, "log2_aFC_sd", "log2_aFC_se/sd", skip_absent=T)
##     }
    
## },
##                    a=common.all.l,
##                    b=c("internal_normalization",
##                        rep("offset_lib.size/GC content",3),
##                        "covariate_log(libsize)",
##                        rep("covariate_lib.szie/GC content",3),
##                        "covariate_log(libsize)"),
##                    SIMPLIFY=F),
##                    idcol="Method", fill=T)

## write.table(long.all, "~/rds/rds-cew54-wallace-share/Projects/trecase/data/method.comparison.txt",
##             row.names=F)
