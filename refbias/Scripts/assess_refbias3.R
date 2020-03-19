#' ---
#' title: Geuvardis eQTL reference bias correction2 and corrected likelihood
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
source("/home/ev250/Bayesian_inf/trecase/Functions/various.R")

## get wildcards to identify files 
source <- snakemake@params[['source']]
rbias <- snakemake@params[['rbias']]

btrec_dir <- snakemake@params[['btrec_dir']]


btrec <- lapply(source, function(i) {
    tmp <- lapply(rbias, function(j) {
        ## For GT both trec and trec-ase were run, select only entries for trec-ase when both are available
            
        tmp <- comb.files(path=paste0(btrec_dir,"/", i), pattern=paste0("^",j, "\\.ENSG[0-9]+.*stan.summary.txt"))
        ## remove bad runs as recommended in stan
        tmp <- tmp[Rhat < 1.01,]
        ## add 95% CI based null column (null.95)
        tmp <- add.null(tmp)
                   })
    names(tmp) <- rbias
    return(rbindlist(tmp, idcol="rbias", fill=T))})
names(btrec) <- source

## Get tags run with GT to combine with RNA

gt.tags <- comb.files(path=paste0(btrec_dir,"/", source[1]), pattern=paste0("^",rbias[1], "\\.ENSG[0-9]+.*tags.lookup.txt"))
gt.tags.run <- merge(gt.tags, btrec[[ source[1]]][rbias=="refbias",.(Gene_id,tag)], by= c("Gene_id","tag"))

## get start end fo genes
## gene.coord <- fread("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/gene_inputs/gene_info.txt")

gene.coord <- fread(snakemake@input[['geneStEnd']]) 

## select genes in chrom 22

gt22 <- gene.coord[chrom==22,]

## add tag distance to gene (closest to start or end)


btrec <- lapply(btrec, function(i) gene.d(dt1=i, dt2=gt22[, .(gene_id, start,end,chrom)]))

############################################################################################################################
######## Compare gt vs rna within 100kB window #############################################################################

## merge gt and rna considering different tags, effect size in RNA already corrected by op.dir, only use trec-ase to compare same model

gt.rna <- comp.ngt(path.s=paste0(btrec_dir, "/",  source[2]), #'/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/refbias2/rna99',
                   pattern.s = paste0("^",rbias[2], "\\.ENSG[0-9]+.*stan.summary.txt"),
                   pattern.t= paste0("^",rbias[2], "\\.ENSG[0-9]+.*tags.lookup.txt"), #"eqtl.tags.lookup.txt",
                   lm.sum=btrec[["GT"]][rbias == "refbias" & model=="trec-ase", ] ,
                   gt.sum= btrec[["GT"]][rbias == "refbias" & model == "trec-ase", ],
                   s=c(".rem", ".gt", ".rna"),
                   tags.m2=gt.tags.run)


## effect size in rna is based on gt tag maf (op.dir column)

## remove ".rem" col, artifact to run conp.ngt function

gt.rna[ , grep("rem$", names(gt.rna), value=T):=NULL]


## exclude inconsistencies due to changes in eaf being too close to 0.5 

gt.rna[(null.99.gt=='no' & null.99.rna == 'no' & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean.rna)), .(tag.gt, tag.EAF.gt, tag.rna, tag.EAF.rna, Gene_id)]

gt.rna<- gt.rna[!(null.99.gt=='no' & null.99.rna == 'no' & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean.rna)),]


## Select association within 100kB and plot

sig <- c(99,95)
sig <- matrix( c(rep(sig, each=2), rep(sig,each=1)), ncol=2)
d <- 10^5

## add null.95.rna col to gt.rna

gt.rna <- add.null(gt.rna, suffix=".rna")

## check and exclude inconsistencies when eaf in gt or rna is too close to 0.5 and are in opposite directions

gt.rna<- gt.rna[!(null.95.gt=='no' & null.95.rna == 'no' & sign(log2_aFC_mean.gt) != sign(log2_aFC_mean.rna)),]

## save gt.rna

write.table(gt.rna, snakemake@output[['gt_rna_sum']], row.names=F)

cols <- c(None="#999999", `obs-GT`="yellow3", `hidden-GT`="#0072B2", Both= "#D55E00")

rna.gt.Sig <- apply(sig,1, function(i) {
    
    gt.rna<- add.signif(gt.rna, x1=paste0("null.", i[2], ".gt"), x2=paste0("null.",i[1], ".rna"), col=c("obs-GT","hidden-GT"))
    
    gt.rna.tab <- tab2bplot(dt=gt.rna[ abs(gene.dist.gt) <= d ,], colors=cols)

    rna.gt <- btrecase.plot(dt=gt.rna[abs(gene.dist.gt) <= d & Signif != "None" ,],
                        x1=c('log2_aFC_mean.gt', 'log2_aFC_2.5%.gt','log2_aFC_97.5%.gt', paste0('null.',i[2],'.gt')),
                        x2= paste0(c('log2_aFC_mean', 'log2_aFC_0.5%','log2_aFC_99.5%', paste0('null.', i[1])), ".rna"),
                        xl='eQTL-effect (observed GT)', yl='eQTL-effect (hidden GT)',
                        col=c("obs-GT","hidden-GT"),axis.title=12, axis.text=10,
                        legend.title=12, legend.text=10,
                        legend.symbol=4, point.size=3 ,
                        title=paste0("Observed (",i[2] ," %CI) vs hidden genotypes (",i[1]," %CI)") , title.size=12) +
    annotation_custom(tableGrob(gt.rna.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(gt.rna.tab$color, rep("black", 4)))))), xmin=-1.3, xmax=-.5, ymin=0.5, ymax=1)

    return(rna.gt)
    })

#+ fig.width= 5, fig.height=9
plot_grid(plotlist=rna.gt.Sig, nrow=3)

############################################################################################################################
##################################################################################################
## Compare gt trec, trec-ase and rna in terms of associations tested and Signif by gene distance

## get trec-ase entries for  gt100 to format long to look at SNP distance to gene, all tested and signif


## only those associations with ASE info
## gt.rna.long <- rbindlist(lapply(btrec, function(i) i[rbias == "refbias" & model != "trec",]), idcol="source", fill=T)


## all associations
gt.rna.long <- rbindlist(lapply(btrec, function(i) i[rbias == "refbias" ,]), idcol="source", fill=T)


gt.rna.long <- gt.rna.long[abs(gene.dist) <= d,]

## Add Signif associations only for GT and RNA to make facet plot, indicate that in source col, use 99%CI 

gt.rna.sig <- rbindlist(lapply(unique(gt.rna.long$source), function(i) {
    temp <- gt.rna.long[source == i & null.99 =="no",]
    temp[, source:=paste0(source, ".Sig")]
    return(temp)
    }))


gt.rna.long <- rbind(gt.rna.long,gt.rna.sig)

## recode to have facet by Sig/All overlapping density for GT and RNA

gt.rna.long[, type:="All"][grep("Sig", source), type:="Signif"]
gt.rna.long[grep("GT", source), source:="obs-GT"][grep("RNA", source), source:="hidden-GT"]

## prepare tables to add to plot

table2d <- lapply(unique(gt.rna.long$type), function(i) {
    dt <- gt.rna.long[type==i,.N,source][order(1/N),]
    names(dt) <- c("", "SNPs")
    return(dt)
    
})


## change table2d[[2]][2,2] to make it compatible with gt.rna (tags issue)

#table2d[[2]][2,2] <- nrow(gt.rna[ Signif =="Both" | Signif == "hidden-GT" ,])



gl <- lapply(table2d, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c("yellow3", "#0072B2", rep("black", 2)))))))


dt <- data.table(type=unique(gt.rna.long$type), grob = gl )

#+ fig.width= 5, fig.height=4
ggplot(gt.rna.long, aes(abs(gene.dist)/1000, color=source, fill=source))  +      
    theme_bw() +
    xlab("Distance to gene (kB)") +
    ylab("Density") +
    theme(legend.title = element_blank()) +     
    guides(colour=guide_legend(override.aes = list(size = 3))) +
    geom_density(alpha = .1) +
    scale_color_manual(values=c("#0072B2","#F0E442")) +
    scale_fill_manual(values=c("#0072B2","#F0E442")) +
    facet_wrap(~type) +
    #facet_grid(cols= vars(type)) + #, "#D55E00"))
    theme(strip.text.x = element_text(size=12),
          strip.background=element_rect(fill="white")) + 
    geom_custom(data=dt, aes(grob=grob), x = 0.6, y = 0.9)
  


#########################################################################################################################
###########################  Compare obs-genotypes with ref bias correction to trec #########################################


trec <- comb.files(snakemake@params[['trec']], pattern="trec.stan.summary.txt")

trec.refbias <- merge(trec, btrec$GT[rbias=="refbias" & model =="trec-ase",], by=c("Gene_id", "tag"), suffixes=c(".trec",".ase"))
trec.refbias <- add.signif(trec.refbias, x2="null.95", x1= "log2_aFC_null" , col=c("trec", "trec-ase"))

trec.rb.tab <- tab2bplot(trec.refbias, colors=c(None="#999999", `trec-ase`="yellow3", trec="#0072B2", Both= "#D55E00"))

btrecase.plot(dt=trec.refbias[Signif != "None",] ,
              x2=c(paste0(c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%"), ".ase"), "null.95"),
              x1=c(paste0(c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%"), ".trec"),  "log2_aFC_null" ),
              yl="eQTL effect btrec-ase",
              xl="eQTL effect trec",
              col=c("trec-ase", "trec"),
              title="eQTL estimates"
              ) + 
    annotation_custom(tableGrob(trec.rb.tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(trec.rb.tab$color, rep("black", 4)))))), xmin=-2, xmax=-1, ymin=1, ymax=2.5)


                   
###################################################################################################################################
############################## Explore outliers 

outl <- gt.rna[Signif== "hidden-GT" & log2_aFC_mean.rna >0 & log2_aFC_mean.gt<0,][3,]

##            Gene_id      tag.rna       tag.gt  log2_aFC_mean.gt
## 1: ENSG00000128335 36616445:T:C 36608160:C:T    -0.1541765

gene <- outl$Gene_id
gt.tag <- outl$tag.gt
rna.tag <- outl$tag.rna

gt.bias <- readRDS(paste0(btrec_dir, "/GT/refbias.",gene,".GT.stan1.input.rds"))
gt.bias <- gt.bias[[gt.tag]]

#' ## Plot observed GT input
inp.qc.gt(gt.bias)



rna.bias <- readRDS(paste0(btrec_dir, "/RNA/refbias.",gene,".noGT.stan.input.rds"))
rna.bias <- rna.bias[[rna.tag]]

#' Plot hidden GT input
inp.qc.nogt(rna.bias)


## Looking at the inputs, p(Gi== het) for those individuals that are hets for the rSNP is about 0.98. p(Gi==0) for those inds that are Gi==0 (obs genotypes) is mostly 0.81 and p(Gi==het) is 0.17.

## get true hets based on obs genotypes

idx.hets <- which(abs(gt.bias$gm$g.ase)==1)

hets <- names(gt.bias$n)[idx.hets]

## select true hets from rna.bias and plot

rna.b.hets <- copy(rna.bias)
rna.b.hets$ase <- lapply(rna.bias$ase, function(i) i[hets])


## Plot with same NB side but only including "True" hets for ASE side

#' Plot hidden GT only including "True" hets for ASE side
inp.qc.nogt(rna.b.hets)


## select true homs from rna.bias and plot

rna.b.hom <- copy(rna.bias)
rna.b.hom$ase <- lapply(rna.b.hom$ase, function(i) {
    w <- which(names(i) %in% names(gt.bias$n)[-idx.hets])
    return(i[w])
    })

## Plot with same NB side but only including "True" homs for ASE side

#' Plot hidden-GT only including "True" homozygotes for ASE side
inp.qc.nogt(rna.b.hom)


## It appears that the swing towards the positive logit(ASE) sign when genotypes are hidden is coming from the ASE contribution of the "True" negative homozygous.




####################################################################################################################################
################################# reference panel bias ##########################################################################


## Compare eQTL estimates with and without ref panel bias correction

## Count significant associations with or without reference panel bias correction


null  <- "null.95"
    
r.bias.95 <- lapply(btrec, function(i) {
    u <- unique(i$rbias)
    names(u) <- c("Without", "With")
    x <- paste(null, u, sep=".")
    i <- i[model != "trec",]
    dt <- merge(i[rbias == u[1],], i[rbias ==u[2], ], by=c("Gene_id", "tag", "tag.EAF", "gene.dist"), suffixes= paste0(".",u))
    dt <- add.signif(dt, x1=x[1], x2=x[2], col=names(u))
    return(dt)
    })

cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%", null)
tab.colors <- c(None="#999999", Without="yellow3", With="#0072B2", Both= "#D55E00")

GT <- r.bias.95$GT
table1 <- tab2bplot(dt=GT, colors=tab.colors)

r.95 <- btrecase.plot(dt=GT[Signif != "None" , ] , x1=paste0(cols, ".norefbias"),
              x2=paste0(cols, ".refbias"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title=paste0("eQTL estimates with or without reference panel bias correction\n", gsub("null\\.", "", null),"% CI")
              ) +
    annotation_custom(tableGrob(table1[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table1$color, rep("black", 4)))))), xmin=-1.0, xmax=-.5, ymin=1, ymax=2)

null  <- "null.99"

r.bias.99 <- lapply(btrec, function(i) {
    u <- unique(i$rbias)
    names(u) <- c("Without", "With")
    x <- paste(null, u, sep=".")
    i <- i[model != "trec",]
    dt <- merge(i[rbias == u[1],], i[rbias ==u[2], ], by=c("Gene_id", "tag", "tag.EAF", "gene.dist"), suffixes= paste0(".",u))
    dt <- add.signif(dt, x1=x[1], x2=x[2], col=names(u))
    return(dt)
    })

cols <- c("log2_aFC_mean", "log2_aFC_0.5%", "log2_aFC_99.5%", null)


GT <- r.bias.99$GT

table1 <- tab2bplot(GT, tab.colors)

r.99 <- btrecase.plot(dt=GT[Signif != "None" , ] , x1=paste0(cols, ".norefbias"),
              x2=paste0(cols, ".refbias"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title=paste0("eQTL estimates with or without reference panel bias correction\n", gsub("null\\.", "", null),"% CI")
              ) +
    annotation_custom(tableGrob(table1[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table1$color, rep("black", 4)))))), xmin=-1.0, xmax=-.5, ymin=1, ymax=2)

#+ fig.width= 5, fig.height=7
plot_grid(r.95,r.99, nrow=2)

#################### Same plot for hidden-GT ############################

RNA <- r.bias.99$RNA

null <- 'null.99'

cols <- c("log2_aFC_mean", "log2_aFC_0.5%", "log2_aFC_99.5%", null)
table1 <- tab2bplot(RNA, tab.colors)

#+ fig.width= 5, fig.height=5
btrecase.plot(dt=RNA[Signif != "None"] , x1=paste0(cols, ".norefbias"),
              x2=paste0(cols, ".refbias"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates with or without\n reference panel bias correction\n hidden-GT"
              ) +
    annotation_custom(tableGrob(table1[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table1$color, rep("black", 4)))))), xmin=-.8, xmax=-.5, ymin=0.2, ymax=0.5)


##############################################################################################
##' # Look at effect size of cis-SNP stratified by being fSNP with genotypes

## start with GT, look for significant eQTL in close LD with a fSNP. Either the fSNP is eQTL or tagged at 0.9

## Get fSNPs and their target gene

## Get fsnps used in GT by gene

files <- list.files(path=paste0(btrec_dir, "/GT"),
                    pattern="^refbias\\.ENSG[0-9]+.GT.fsnps.with.counts.rds",
                    full.names=T)

u.fsnps <- lapply(files, function(i) unique(gsub("\\.[mn]$","",  readRDS(i))))
names(u.fsnps) <- sapply(files, function(i) regmatches(i, regexpr("ENSG[0-9]+", i)))

dt.fsnps <- data.table(Gene_id=unlist(lapply(seq_along(u.fsnps), function(i) rep(names(u.fsnps)[i], length(u.fsnps[[i]])))),
                       fSNP=unlist(u.fsnps))

## get tagging SNP for fsnps
tag.fsnp <- merge(dt.fsnps, gt.tags.run, by.x=c("Gene_id", "fSNP"), by.y=c("Gene_id", "SNP"))

## add EAF for SNP col to check direction of effects with tag
DT <- snp.eaf(file1='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz', snps = tag.fsnp$fSNP)
tag.fsnp <- merge(tag.fsnp, DT, by.x="fSNP", by.y="snp")
setnames(tag.fsnp, c("fSNP","eaf"), c("fSNP.id","EAF.fsnp"))


gt.l <- merge(btrec$GT[model == "trec-ase",], tag.fsnp, by=c("Gene_id", "tag"), all.x=T)
gt.l[, fSNP:="yes"][is.na(EAF.fsnp), fSNP:="no"]
gt.l[fSNP=="yes" , tag.fsnp.op:="no"]
gt.l[(EAF.fsnp < 0.5 & tag.EAF > 0.5) | (EAF.fsnp > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]

## create source column to use in plot to distinguish refbias-ASE, norefbias-ASE and trec, to be added below

gt.l[,source:="Ref. bias correction"][rbias=="norefbias", source:="Without correction"]

## select the same Gene-snp pairs in refbias and norefbias, use r.bias (95 or 99) to select common entries

gt.l <- merge(gt.l, r.bias.95$GT[ , .(Gene_id,tag,model.refbias)], by=c("Gene_id", "tag"))

gt.l[, model.refbias:=NULL]

## add trec

## get trec with same tags as gt.l run with ASE info
trec.sub <- merge(trec, unique(gt.l[,.(Gene_id, tag, tag.EAF)]), by=c("Gene_id", "tag"))
trec.sub[, source:="NegBinomial only"]

## trec.sub has same tags as gt.l, compare same associations
trec.sub <- merge(trec.sub, tag.fsnp, by=c("Gene_id", "tag"), all.x=T)

trec.sub[, fSNP:="yes"][is.na(EAF.fsnp), fSNP:="no"]
trec.sub[fSNP=="yes" , tag.fsnp.op:="no"][(EAF.fsnp < 0.5 & tag.EAF > 0.5) & (EAF.fsnp > 0.5 & tag.EAF < 0.5), tag.fsnp.op:="yes"]


## combine with gt.l
gt.trec <- rbind(gt.l, trec.sub, fill=T)
gt.trec[, Signif:="No significant"][log2_aFC_null=="no", Signif:="Significant"][null.95 == "no", Signif:="Significant"]


## Some fSNPs are in different direction as tag SNP, create new column to adjust eQTL estimates

gt.trec[, log2mean.eaf:=log2_aFC_mean][tag.fsnp.op=="yes", log2mean.eaf:= -log2_aFC_mean]
gt.trec[,source:=factor(source, levels=c( "Without correction" ,  "Ref. bias correction", "NegBinomial only"))]

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

## Check I am comparing same rSNPs in each condition

gt.trec[,.N,source]

#' # eQTL effects by fSNP
ref.p


##' # Distribution of eQTL effects by model

## Add column to indicate in all models whether the association is significant in Neg Binomial

gt.trec<- merge(gt.trec, unique(gt.trec[source == "NegBinomial only" ,.(Gene_id, tag, Signif)]), by =c("Gene_id", "tag"), suffixes=c("", ".NegBinom"))

gt.trec[, Signif.NegBinom:=paste(Signif.NegBinom, "in Neg Binom model")]


tabs <- tab2facet(dt=gt.trec, vars=c("Signif.NegBinom", "source"))

gl <- lapply(tabs, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 7,padding = unit(c(1.5, 1,1,1), "mm"), core=list(fg_params=list(col=rep(pal, 2))))))   

gt.trec <- gt.trec[order(match(source, levels(gt.trec[["source"]]))),]

## make sure names in dt are the same as in gt.trec (Signif and source)
gnb <- aux.tab2facet(gt.trec, c("Signif.NegBinom", "source"))
dt <- data.table(Signif.NegBinom=gnb$Var1, source=gnb$Var2, grob = gl )

ref.p2 <- ggplot(gt.trec, aes(log2_aFC_mean)) +
    xlab("eQTL-effect") +
    ylab("Density") +
    geom_density()+
    scale_color_manual(values=pal) +
    theme_bw() +
    geom_vline(xintercept= 0, linetype = "dashed",  color = "grey")+
    facet_grid(rows= vars(source), cols=vars(Signif.NegBinom), scales = "free_x")+
    theme(strip.background=element_rect(fill="white"))+
    geom_custom(data=dt, aes(grob=grob), x = 0.75, y = 0.75)


##' # Plot distribution of eQTL effects by significance in trec

ref.p2
##ref.norefbias 
##ref.refbias 

##' # Distribution of the difference between uncorrrected trecase and trec eQTL estimates and corrected/trec with genotypes


diff <- rbindlist(lapply(rbias, function(i) {
    dt <- merge(unique(gt.trec[rbias==i,.(Gene_id, tag, log2_aFC_mean, Signif.NegBinom)]), unique(gt.trec[model == "trec",.(Gene_id, tag, log2_aFC_mean)]), by =c("Gene_id", "tag"), suffixes = c(paste0(".", i), ".trec"))
    dt[ ,diff.with.trec := get(paste0("log2_aFC_mean.", i))- log2_aFC_mean.trec]
    dt[, grep("log2_aFC_mean.", names(dt), value=T):=NULL]
    dt[,RefBiasCorrection:=paste0("Obs-GT.", i )]
    return(dt)
}))

g <-  aux.tab2facet(diff, "Signif.NegBinom")

tabs <-  lapply(1:nrow(g), function(i) {
        
    dt <- diff[Signif.NegBinom==g[i,Var1]  ,.N,.(RefBiasCorrection, sign(diff.with.trec))]
        dt[, sign:=as.character(sign)][sign=="-1", sign:="<0"][sign=="1", sign:=">0"]
        setkey(dt, RefBiasCorrection, sign)
        setnames(dt, c("sign", "N"), c("Effect", "SNPs"))
        return(dt)
    })

gl <- lapply(tabs, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 7,padding = unit(c(1.5, 1,1,1), "mm"), core=list(fg_params=list(col=c( rep(pal, each=2), rep("black", 8)))))))   

##diff <- diff[order(match(source, levels(g$Var2))),]

## make sure names in dt are the same as in diff 
dt <- data.table(Signif.NegBinom=g$Var1,  grob = gl )


#+ fig.width= 6, fig.height=7
ggplot(diff, aes(diff.with.trec, color=RefBiasCorrection)) +
    xlab("eQTL-effect") +
    ylab("Density") +
    geom_density()+
    scale_color_manual(values=pal) +
    theme_bw() +
    geom_vline(xintercept= 0, linetype = "dashed",  color = "grey")+
    facet_grid(rows= vars(Signif.NegBinom), scales = "free_y")+
    theme(strip.background=element_rect(fill="white")) +
    ggtitle("Distribution of differences of uncorrected/corrected\n B-trecase relative to B-trec") #+
    #geom_custom(data=dt, aes(grob=grob), x = 0.75, y = 0.75)


#################### Same plot for hidden-GT ############################


##' # Distribution the difference between uncorrrected trecase and trec eQTL estimates and corrected/trec without genotypes

## For each tag in gt randomly select one equivalent rna tag
seed <- 1
gt.rna.one <- gt.rna[, .(tag.rna=sample(x=tag.rna, size=1)), .(Gene_id,tag.gt, tag.EAF.gt)]

## Select tag.rna as in gt.rna.one from RNA

rna.sub <- merge(gt.rna.one, RNA, by.y=c("Gene_id", "tag"), by.x=c("Gene_id", "tag.rna"))

## Add trec from gt.trec to have Signif.NegBinom col

trec.rna.sub <- merge(rna.sub, unique(gt.trec[model=="trec", .(Gene_id, tag,log2_aFC_mean, Signif.NegBinom )]) , by.x=c("Gene_id","tag.gt"), by.y=c("Gene_id", "tag"), suffixes=c(".trec", ""))


## calculate the difference between uncorrrected trecase and trec eQTL estimates and corrected/trec without genotypes

diff.rna <- rbindlist(lapply(rbias, function(i) {
    col <- paste0("log2_aFC_mean.", i)
    dt <- trec.rna.sub[, c("Gene_id", "tag.gt", "tag.rna", "tag.EAF.gt", "tag.EAF", col, "log2_aFC_mean", "Signif.NegBinom" ), with=F]
    ## correct direction of effects in rna if tags in op directions
    dt[(tag.EAF.gt <0.5 & tag.EAF >0.5) | (tag.EAF.gt >0.5 & tag.EAF <0.5), eval(col):= -get(col)]
    dt[, diff.with.trec:= get(col) - log2_aFC_mean ]
    dt[ , grep("log2_aFC_mean", names(dt), value=T):=NULL]
    dt[,RefBiasCorrection:= paste0("Hidden-GT.", i )]
    return(dt)
}))
       
## Combine diff with  diff.rna, selecting same tags.gt

diff.all <- rbind(merge( gt.rna.one, diff, by.x=c("Gene_id","tag.gt"), by.y=c("Gene_id", "tag")),
                  diff.rna, fill=T)

g <-  aux.tab2facet(diff.all, "Signif.NegBinom")
pal <- brewer.pal(n = 8, name = "Set1")[1:4]
tabs <-  lapply(1:nrow(g), function(i) {
        
    dt <- diff.all[Signif.NegBinom==g[i,Var1]  ,.N,.(RefBiasCorrection, sign(diff.with.trec))]
        dt[, sign:=as.character(sign)][sign=="-1", sign:="<0"][sign=="1", sign:=">0"]
        setkey(dt, RefBiasCorrection, sign)
        setnames(dt, c("sign", "N"), c("Effect", "SNPs"))
        return(dt)
    })

gl <- lapply(tabs, function(i) tableGrob(i, rows=NULL, theme=ttheme_minimal(base_size = 7,padding = unit(c(1.5, 1,1,1), "mm"), core=list(fg_params=list(col=c( rep(pal, each=2), rep("black", 2*nrow(i))))))))   

##diff <- diff[order(match(source, levels(g$Var2))),]

## make sure names in dt are the same as in diff 
dt <- data.table(Signif.NegBinom=g$Var1,  grob = gl )


#+ fig.width= 6, fig.height=7
ggplot(diff.all, aes(diff.with.trec, color=RefBiasCorrection)) +
    xlab("eQTL-effect") +
    ylab("Density") +
    geom_density()+
    scale_color_manual(values=pal) +
    theme_bw() +
    geom_vline(xintercept= 0, linetype = "dashed",  color = "grey")+
    facet_grid(rows= vars(Signif.NegBinom),  scales = "free_y")+
    theme(strip.background=element_rect(fill="white")) +
    ggtitle("Distribution of differences of uncorrected/corrected\n B-trecase relative to B-trec\n
with observed or hidden GT") #+
    #geom_custom(data=dt, aes(grob=grob), x = 0.25, y =0.75)

## get quantiles for distributions
probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)

qd <- lapply(unique(diff.all$Signif.NegBinom), function(i) sapply(unique(diff.all$RefBiasCorrection),
                                           function(j) quantile(diff.all[RefBiasCorrection == j & Signif.NegBinom == i, diff.with.trec],
                                                                probs=probs)
                                           ))

names(qd) <- unique(diff.all$Signif.NegBinom)

qd

## btrecase.plot(trec.rna.sub[null.95.norefbias=="no" & Signif.NegBinom=="Significant in Neg Binom model",], x1=paste0(c("log2_aFC_mean",  "log2_aFC_2.5%" , "log2_aFC_97.5%", "null.95" ), ".norefbias"),
##               x2=c(rep("log2_aFC_mean",3),"Signif") ,
##               xl="eQTL effect without correction",
##               yl="eQTL effect trec",
##               col=c("Without", "trec"),
##               title="eQTL estimates with trec or btrecase without reference panel bias correction"
##               )
## ## #######################################################################################################################################
############################### comparing the variability of running the same Baysian model twice ####################################


## With and without reference panel bias run trec model using same inputs, for those associations that dont have ASE info

trec.var <- lapply(btrec, function(i) {
    u <- unique(i$rbias)
    names(u) <- c("Without", "With")
    x <- paste("null.95", u, sep=".")
    i <- i[model == "trec",]
    dt <- merge(i[rbias == u[1],], i[rbias ==u[2], ], by=c("Gene_id", "tag", "tag.EAF", "gene.dist"), suffixes= paste0(".",u))
    dt <- add.signif(dt, x1=x[1], x2=x[2], col=names(u))
    return(dt)
    })

cols <- c("log2_aFC_mean", "log2_aFC_2.5%", "log2_aFC_97.5%", "null.95")


trec.var <- trec.var$GT
table1 <- tab2bplot(trec.var, tab.colors)



##' Variation in estimates for the same model run twice, negative binomial side is not affected by the refrence panel bias correction
btrecase.plot(dt=trec.var[Signif != "None"] , x1=paste0(cols, ".norefbias"),
              x2=paste0(cols, ".refbias"),
              xl="eQTL effect without correction",
              yl="eQTL effect with corrrection",
              col=c("Without", "With"),
              title="eQTL estimates for Trec with or without reference panel bias correction"
              ) +
    annotation_custom(tableGrob(table1[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 10,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(table1$color, rep("black", 4)))))), xmin=-5, xmax=-1.5, ymin=2, ymax=5)


############################################################################################################
########### Examples of reference panel bias #################

GT <- r.bias.99$GT

## merge GT with trec to select examples for refbias figure

GT.trec <- merge(GT, trec, by=c("Gene_id", "tag"))

ref.loss <- GT.trec[Signif=="Without" &  log2_aFC_null == "yes", .(Gene_id,tag, log2_aFC_mean.norefbias, log2_aFC_mean.refbias, log2_aFC_mean, log2_aFC_null)][31,]


gene <- ref.loss$Gene_id
tg <- ref.loss$tag

ref.correct <- readRDS(paste0(btrec_dir, "/GT/refbias.",gene,".GT.stan1.input.rds"))

inp.qc.gt(ref.correct[[tg]])


ref.loss

ref.gain <- GT.trec[Signif=="With" &  log2_aFC_null == "no", .(Gene_id,tag, log2_aFC_mean.norefbias, log2_aFC_mean.refbias, log2_aFC_mean, log2_aFC_null)][18,]

gene <- ref.gain$Gene_id
tg <- ref.gain$tag

ref.correct <- readRDS(paste0(btrec_dir, "/GT/refbias.",gene,".GT.stan1.input.rds"))

inp.qc.gt(ref.correct[[tg]])


ref.gain



##' # Conclusions

## There is good correlation of effects between obs-GT and hidden-GT.

## The false positives I have examined appears to be caused by a systematic
##increase in allelic imbalance from "wrong"genotype

## Re-running the model with observed genotypes after excluding ASE double counting
##resulted in very subtle difference with and without reference panel bias correction.

## In this setting, comparing Bayesian-trec with Bayesian trecase with or without
## reference panel bias correction didnt produce much evidence of reference panel
## bias. I need to re-think how I am going to motivate the reference panel bias correction. 

## Same analysis with RNA doesnt look good.

## The variability found when running the same model twice is low, I run it
## at 95% CI to have more associations.
