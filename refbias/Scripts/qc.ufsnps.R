source('/home/ev250/Bayesian_inf/trecase/Functions/aux_btrecase.R')

## get inputs

rbias <- snakemake@params[['rbias']]
##rbias=c("norefbias","rbias")

## btrecase_dir  <- paste0("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/", c("GT","RNA"))
btrecase_dir <- snakemake@params[['btrecase']]

ufsnp <- fread(snakemake@input[['ueSNPs']])
## ufsnp <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/inputs/fSNP/chr22.fSNP.unique.genes.txt")
ufsnp[,id:=paste(POS, REF,ALT, sep=":")]

fsnps <- fread(snakemake@input[['eSNPs']])
fsnps <- fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt')
fsnps[,id:=paste(POS, REF,ALT, sep=":")]

## Look at genes with fsnps shared

fsnps.shared <- fsnps[,.N, id][N>1,id]
genes.fshared <- unique(fsnps[id %in% fsnps.shared, gene_id])

## Check if were run with ASE in GT or run in rna. If so, check if the fSNPs that should have been removed were removed (fsnps.with.counts.rds or fisher in noGT)

tocheck <- lapply(btrecase_dir, function(i) {
    tmp <- lapply(rbias, function(j) {      
        tmp <- comb.files(path=i, pattern=paste0("^",j, "\\.ENSG[0-9]+.*stan.summary.txt"))      
                   })
    names(tmp) <- rbias
    dt <- rbindlist(tmp, idcol="rbias", fill=T)
    dt <- dt[Gene_id %in% genes.fshared,]
    return(dt)})

names(tocheck) <- basename(btrecase_dir)

lapply(tocheck, function(i) i[,.N, null.99])

lapply(tocheck, function(i) i[model != "trec",.N,.(rbias, Gene_id)])

## none significant associations, make sure aux.2 is working well 

genes.run <- lapply(tocheck, function(i) i[model != "trec",.N,Gene_id][, Gene_id])

## look at the fSNPs used in model and check if they were unique

fisher <- rbindlist(lapply(list.files(btrecase_dir[2], pattern=paste0("bias\\.ENSG[0-9]+.*fsnps.het.fisher.test.txt"), full.names=T), fread))
fisher <- fisher[gene_id %in% genes.fshared,]

## aux.2 okay here
fisher$fsnp %in% fsnps.shared

files <- list.files(btrecase_dir[1], pattern=paste0("bias\\.ENSG[0-9]+.*GT.fsnps.with.counts.rds"), full.names=T)
gt.frun <- lapply(files, readRDS)
names(gt.frun) <- basename(files)

gt.frun <- gt.frun[unlist(lapply(genes.fshared, function(i) grep(i, names(gt.frun), value=T)))]

## unlist and check that all are in ufsnp

gt.frun <- unique(unlist(gt.frun))

## remove ending .n .m from fsnps and then check if they are listed as unique
u <- lapply(unique(lapply(gt.frun, function(i) gsub(".[n:m]$", "", i))), function(j) ufsnp[id==j, gene_id])

## remove ending from fsnps run in model (ending .m .n) and check if listed for more than one gene in fsnps
lplus <- lapply(unique(lapply(gt.frun, function(i) gsub(".[n:m]$", "", i))), function(j) fsnps[id == j, gene_id])

## problematic are those elements in lplus with more than 1 gene, select them and identifiy source gene (should be in genes.run)

genes2rerun <- unique(lapply(lplus[unlist(lapply(lplus, function(i) length(i) >1))], function(j) genes.run$GT[j %in% genes.run$GT]))

write.table(data.table(Gene_id=genes2rerun), paste0(snakemake@params, "/GT.rerun.txt"), row.names=F)

## issue with  "ENSG00000100056", make sure now aux2 works well, no need to re-run as assoc are not significant, but should only run with 2 fsnps instead of 6. aux2 corrected.
