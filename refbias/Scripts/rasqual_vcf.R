source("/home/ev250/Bayesian_inf/trecase/Functions/real.data.R")


#' Format vcf to accommodate rasqual AS field for hom individuals
#'
#' @param vcf full name with vcf file made for BaseQTL
#' @param allelic.files name of files with allelic counts per sample
#' @param out name for file to output
#' @export
#' @return saves a data table with vcf format to run rasqual
#' vcf.rasq()

vcf.rasq <- function(vcf, allelic.files,out){
    
    dna <- vcf_w(vcf, f.arg='"%CHROM %POS %ID %REF %ALT %QUAL %FILTER[ %GT] [ %AS]\\n"')
    ## dna2 = copy(dna)
    
    samples <-  gsub("\\.alleleCounts.txt", "", basename(allelic.files))

    ## get list with counts to replace
    ac <- lapply(allelic.files, function(i) {
        dt <- fread(i)
        dt[, id:=paste(POS,REF,ALT, sep=":")]
        ## RA col to add to vcf
        dt[, RA:=paste(NREF,NALT, sep=",")]
        dt[, sample:=gsub("\\.alleleCounts.txt", "", basename(i))]
    })
    names(ac) <-samples

    gt <- paste0(samples, "_GT")
    as <- paste0(samples,"_AS")
    tmp <- lapply(1:length(samples), function(i) {
        hom <- dna[get(gt[i]) == "0|0" | get(gt[i]) == "1|1", id]
        ## add RA col to homs
        rep <- ac[[i]][id %in% hom, .(id,RA)]
        dna[id %in% rep$id, eval(as[i]) := rep$RA]
        dna[, eval(gt[i]):=paste(get(gt[i]), get(as[i]), sep=":")]
    })

    ## format dna for vcf and save
    ## required cols:  CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO and FORMAT followed by samples, samples is GT:AS

    dna[, INFO:="."][, FORMAT:="GT:AS"]
    dna[, c("id", paste0(samples, "_AS")):=NULL]
    setcolorder(dna, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", gt))
     write.table(dna, out,row.names=F,col.names=F,quote=F,sep="\t")
}

vcf.rasq(vcf=snakemake@input[['vcf']], allelic.files=snakemake@input[['allelicCounts']] ,out=snakemake@output[['out']])

## vcf <- "/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/DNA/chr22.ID.ASE.allsamples.vcf.gz"#, f.arg='"%CHROM %POS %ID %REF %ALT[ %GT]\\n"')
## allelic.files <- list.files("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/allele_counts",full.names=T)


