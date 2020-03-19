library(data.table)



## open dna files

dna.m <- vcf_w('/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/inputs/corrected/rasq.vcf.gz')

dna.rasq <- vcf_w("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/rasqual/inputs/rasq_ase/chr22.86.vcf.gz")

fsnp <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt")

fsnp[,id:=paste(POS, REF,ALT, sep=":")]

f.m <- merge(dna.m, fsnp, by='id')

f.r <- merge(dna.rasq, fsnp, by='id')

## look for some snps with counts
r <- f.r[!HG01791_AS %in% "0,0", id]

## compare output for the same snps
f.m[id %in% head(r), ]
f.r[id%in% head(r),]

## The difference is that when making rasqual inputs I had only considered fSNPs with maf>0.05 to add hom counts.

