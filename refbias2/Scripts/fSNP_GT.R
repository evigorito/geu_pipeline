library(data.table)
library("GenomicFeatures")
library("GenomicAlignments")

#' Map fSNPs to genes 
#'
#' @param ebg genomics ranges object exon by genes to extract exon coordinates
#' @param gt full name to rds file with SNPs and GT for samples to use, gt in 0,1,-1,2 scale but okay in 0,1,2 scale, extracted from vcf
#' @param het whether to select only snps with at least 1 ind, defaults to 1
#' @param chrom chromosome to extract snps
#' @param out full name to save text file with fSNP coordinates mapped to gene_id
#'
#' @return saves file with mapping info of fSNPs to genes
#' @export
#' gene.snps()

gene.snps <- function(ebg,het=1, gt, chrom, out){
   
    gt <- readRDS(gt)
    ## exclude snps if het >=1
    if(het>0){
        samp <- grep("_GT", names(gt), value=T)
        ## get rows with het  inds >=het
        r <- apply(gt[, samp,with=F], 1, function(i) sum(abs(i)==1))
        gt <- gt[which(r>=het),]
    }
    
    gt[ , ID:=paste(chrom, POS, REF, ALT, sep=":")][, POS:=as.numeric(POS)]
    rebg <- unlist(reduce(readRDS(ebg)))
    
    fsnp_granges = GenomicRanges::GRanges(seqnames =chrom, 
                                          ranges = IRanges::IRanges(start = gt$POS,
                                                                    end = gt$POS))
    values(fsnp_granges)  <- DataFrame(ID = gt[['ID']])
    feature_olaps = GenomicRanges::findOverlaps(rebg, fsnp_granges, ignore.strand=TRUE)
    ## extract matches and rename cols
    genes <- as.data.table(rebg[queryHits(feature_olaps)])
    genes[, gene_id:= rebg[queryHits(feature_olaps)]@ranges@NAMES ]
    snps <- as.data.table(fsnp_granges[subjectHits(feature_olaps)])
    ## add ref and alt alleles to snp
    snps <- merge(snps, gt, by="ID", all.x=T, sort=F)
    ## select columns
    comb <- cbind(snps[,.(seqnames, POS, ID, REF, ALT)], genes[, .(gene_id)] )

    ## keep SNPs only
    comb[,lenref:=nchar(REF)][,lenalt:=nchar(ALT)]
    comb <- comb[lenref==1 & lenalt==1,]

    ## remove header when saving and keep the requested columns by WASP: POS, REF, ALT
    comb[, c("seqnames", "ID", "gene_id","lenref","lenalt"):=NULL]
    comb <- unique(comb)
    setkey(comb, POS)
    comb[, POS:=format(POS, scientific = FALSE) ] ## avoid xe10 nomenclature
    
    write.table(comb, file=out, row.names=F, col.names=F, quote=F)
    
}

gene.snps(ebg=snakemake@input[['ebg']],
          het= as.numeric(snakemake@params[['het']]),
          gt=snakemake@input[['gt']],
          chrom=snakemake@wildcards[['chrom']],
          out=snakemake@output[['fsnps']] )
