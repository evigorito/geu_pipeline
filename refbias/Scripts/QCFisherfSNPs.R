library(parallel)

library(data.table)
source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")
options(mc.cores = parallel::detectCores())


sample.file=snakemake@input[['sample']]
le.file=snakemake@input[['leRef']]
h.file=snakemake@input[['hapRef']]
fsnps=snakemake@input[['ueSNPS']]
vcf=snakemake@input[['vcf']]

population=snakemake@params[['population']]
rna_stan=snakemake@params[['rna_stan']]
dna_dir=snakemake@params[['dna_dir']]
dna_fsnps=snakemake@params[['dna_fsnps']]

out=snakemake@output[['out']]

## sample.file <- '/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3.sample'
## le.file <- '/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz'
## h.file <- '/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz'
## fsnps <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/inputs/fSNP/chr22.fSNP.unique.genes.txt'
## vcf <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/RNA/chr22.ASE.allsamples.vcf.gz'
## population="EUR"
## rna_stan <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/RNA'
## dna_dir <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/Btrecase/SpikeMixV3_2/GT'
## dna_fsnps <- "rbias.ENSG[0-9]+.*.\\.GT.fsnps.with.counts.rds"


#' Get Fisher test for all fSNPs called by RNA-seq and used in inference with GT
#'
#' @param sample.file sample file with reference panel samples, NULL if all samples are used
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity: AFR AMR EAS EUR SAS ALL, defaults to EUR. It is used to exclude fSNPs if frequency is different from reference panel when using no GT and selecting ex.fsnps argument
#' @param vcf path to vcf file with ASE and GT from RNA-seq
#' @param rna_stan path to dir with output for rna_stan to get genes run in model
#' @param dna_dir path to dir with output for stan run with GT to get fSNPs used in inference
#' @param dna_fsnps pattern to identify files with fSNPs used in inference with GT
#' @param fsnps file with fsnps witin each gene (only unique fSNPs)
#' @param out full name to output file
#' @export()
#' @keywords haplotype frequency
#' @return saves file with table rows fsnps and colnames OR (odds ratio) and pvalue for Fisher's exact test
#' fisher.het.test()

fisher.het.test <- function(sample.file=NULL, le.file, h.file, population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), vcf, rna_stan, dna_dir, dna_fsnps, fsnps, out){

    ## get GT called by RNA-seq
    rna <- vcf_w(vcf)
    
    ## get fSNPs run in inference in DNA
    fsnp.i <- lapply(list.files(dna_dir, pattern=dna_fsnps, full.names=T), readRDS)
    fsnp.i <- unique(unlist(fsnp.i))
  
    ##remove .n .m ending
    fsnp.i <- unique(unlist(lapply(fsnp.i, function(i) gsub(".[n:m]$", "", i))))

    ## link fSNPs to genes
    fsnp <- fread(fsnps)
    fsnp[,id:=paste(POS, REF,ALT, sep=":")]

    ## select genes run in RNA
    rna.genes <- list.files(path=rna_stan, pattern="^rbias\\.ENSG[0-9]+.*stan.summary.txt")
    rna.genes <- gsub(".*(ENSG[0-9]+).*", "\\1", rna.genes)

    ## select fSNPs used in GT inference
    fsnp.i <- fsnp[id %in% fsnp.i,][gene_id %in% rna.genes,]
    setkey(fsnp.i, POS)

    

    ## select samples to use

    sam <- fread(sample.file)
    ## get rows in sam  for population
    rows.panel <- sam[GROUP==population, which=T]
    ## samples in hap file are in columns and each sample corresponds to 2 cols
    ## I need to get for each sample in which col of hap file it is.
    cols <- sort(c(rows.panel*2, rows.panel*2-1))

    ## merge rna with fsnp.i to get genes 
    rna <- merge(rna, fsnp.i[,.(gene_id, id)], by="id")
    

    
    het.f <-mclapply(unique(rna$gene_id), function(u) {
        ## get reference panel GT
       
        rp <- haps.range(file1=le.file,file2=h.file,cw=c(min(rna[gene_id==u,POS]), max(rna[gene_id==u, POS])), population, maf=0)

        ## select fSNPs and samples to use
        rp <- rp[rna[gene_id ==u,id],cols, drop=F ]
        
        p <- prop_het(f.ase=rna[gene_id==u & id %in% rownames(rp),],
                 rp.f=rp,
                 gene=u)
        
        return(p)
        
    })
    
    write.table(rbindlist(het.f), out, row.names=F)
}
                                                               
fisher.het.test(sample.file, le.file, h.file, population, vcf, rna_stan, dna_dir, dna_fsnps, fsnps, out)
