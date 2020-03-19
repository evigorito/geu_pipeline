library(data.table)
library("GenomicFeatures")
library("GenomicAlignments")

#' Map fSNPs to genes 
#'
#' @param ebg genomics ranges object exon by genes to extract exon coordinates
#' @param snp full name to legend file from reference panel
#' @param chrom chromosome for legend file
#' @param maf maf cut-off to select SNPs
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL
#' @param out full name to save text file with fSNP coordinates mapped to gene_id
#'
#' @return saves file with mapping info of fSNPs to genes
#' @export
#' gene.snps()

gene.snps <- function(ebg, snp, chrom, maf=0.05, population=c("EUR", "AFR", "AMR", "EAS","SAS", "ALL"), out){
    population <- match.arg(population)
    eth <- 6:11 ## field number in legend file for ethniticy
    names(eth) <- c("AFR", "AMR", "EAS",  "EUR", "SAS", "ALL")
    
    ## get snp info filtering by maf from legend file
    et.field <-  unname(eth[names(eth)==population]) ## population field in leg file
    
    snp.i <- system(paste0("zcat ", snp, " | awk '{if ($", et.field, " >= ", maf, ") {print $2 \" \"$3 \" \"$4 } }' | awk '{if (NR!=1) {print}}' " ), intern=TRUE)  ## second field in legend file is POS,then ref then alt allele, exclude header

    if(length(snp.i)==0) return("no snps in reference panel")
    
    ## format snp.i as DT
     DT <- as.data.table(lapply(1:3, function(i) unlist(lapply(strsplit(snp.i, split=" ") , `[[`,i))))
    names(DT) <- c("POS","REF","ALT")
 
    rebg <- unlist(reduce(readRDS(ebg)))
    DT[ , ID:=paste(chrom, POS, REF, ALT, sep=":")][, POS:=as.numeric(POS)]
    fsnp_granges = GenomicRanges::GRanges(seqnames =chrom, 
                                          ranges = IRanges::IRanges(start = DT$POS,
                                                                    end = DT$POS))
    values(fsnp_granges)  <- DataFrame(ID = DT[['ID']])
    feature_olaps = GenomicRanges::findOverlaps(rebg, fsnp_granges, ignore.strand=TRUE)
    ## extract matches and rename cols
    genes <- as.data.table(rebg[queryHits(feature_olaps)])
    genes[, gene_id:= rebg[queryHits(feature_olaps)]@ranges@NAMES ]
    snps <- as.data.table(fsnp_granges[subjectHits(feature_olaps)])
    ## add ref and alt alleles to snp
    snps <- merge(snps, DT, by="ID", all.x=T, sort=F)
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
          snp=snakemake@input[['legend']],
          maf= as.numeric(snakemake@params[['maf']]),
          chrom=snakemake@wildcards[['chrom']],
          population=snakemake@params[['pop']],
          out=snakemake@output[['fsnps']] )
