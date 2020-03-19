#' ---
#' title: Allelic imbalance estimates for GEUVADIS dataset
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

#' ## Allelic imbalance estimates for fSNPs
#'
#' Open file with estimates

AI <- fread(snakemake@input[['AI']])
## AI <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/pre_post_AI.txt")

## Calculate 99% confidence interval testing whether AI_post is different from 0.5

AI[Total_post >0 , c("lowCI", "hCI") := lapply(1:2, function (i) binom.test(x=NALT_post,
                                                                            n=Total_post,
                                                                            p = 0.5,
                                                                            alternative = "two.sided", conf.level = 0.99)$conf.int[i]),by=.(NALT_post,Total_post) ]

## Add 99% CI relative to 0.5

AI[ , c("null_lCI", "null_hCI") := lapply(c("lowCI", "hCI"), function(i)  0.5 + get(i) -AI_post) ]

AI[, Discard:= ifelse(AI_post > null_hCI, "yes", "no")][Total_pre<100, Discard := "yes"]

AI[,.N,Discard]

## Draw lines for 99% CIs from the null value
ggplot(AI[Total_post>0,] ) +
    geom_line(aes(x=log(Total_pre), y=null_lCI), colour="gray", linetype="dashed") +
    geom_line(aes(x=log(Total_pre), y=null_hCI),  colour="gray", linetype="dashed") +
    geom_point(aes(x=log(Total_pre), y=AI_post, colour=Discard)) +
    geom_vline(xintercept=log(100), linetype="dashed") +
    ylab("AI estimate") +
    ggtitle("Allelic imbalance estimates with \n lines for null 99% confidence intervals")

## get SNPs above null_hCI


AI[Total_pre >100 & Discard == "yes",]

hist(AI$AI_post)
summary(AI$AI_post)

nrow(AI)
nrow(AI[AI_post>0.5,])

## Percentage of reads with AI_post<=0.5

nrow(AI)*100/(nrow(AI)+nrow(AI[AI_post>0.5,]))

#' ## Investigating AI>0.5

## Look into one sample

ai <- fread(snakemake@input[['sample']])
#ai <- fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/post_remap/HG00096/Aligned.sortedByCoord.out.post_remapping_AI.txt")

ai[,AI:= NALT/(NREF+NALT)]

ai[AI==1,]


#' I selected the pair matching to 24000567. That pair overlapped 3 SNPs [[array([24000567, 24000569], dtype=int32)], [array([24000486], dtype=int32)]]

## Read1 is ALT for 24000567, 24000569 and read2 is REF for 24000486

ai[POS>=24000486 & POS<=24000569 ,]

## I generated a faked paired read as REF, REF, ALT (as in array coordinates).
## That pair was unmapped by aligner so it is missing in the output.

AI[POS>=24000486 & POS<=24000569 ,]

## Across samples AI>0.5 is not too bad when considering at least 100 starting reads

AI[AI_post>0.6 & Total_pre > 100,]


#' ## Look at imbalance of fSNPs run in model

stan = fread(snakemake@input[['stan']])

#stan=fread('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/RNAnoGTrsnp/compare.trec.gt.rna.wide.txt')

## Select potentially biased and unbiased fSNPs and compare AI estimates

pot.bias <- stan[Signif != "None" & log2_aFC_mean.gt < 0 & fSNP == "yes" &
                 log2.aFC_null.lm == "no", fSNP.id]
pot.nbias <- stan[Signif != "None" & fSNP == "yes" &
                  log2.aFC_null.lm == "yes", fSNP.id]


## get AI

AI[, id:=paste(POS,REF,ALT, sep=":")]

AI.bias <- lapply(list(pot.bias, pot.nbias),
                  function(i) AI[id %in% i , .(AI_post, Total_pre)])
names(AI.bias) <- c("pot.bias", "pot.unbias")

AI.bias <- rbindlist(AI.bias, idcol="stan_bias")


ggplot(AI.bias, aes(AI_post, fill=stan_bias)) +
    geom_histogram( position="dodge", alpha=1, binwidth=0.005) +
    geom_vline(xintercept=0.5, linetype="dashed") +
    ylab("Frequency") +
    ggtitle("Distribution of AI estimates")

ggplot(AI.bias, aes(x=AI_post, colour=stan_bias)) + geom_density()+
    geom_vline(xintercept=0.5, linetype="dashed") +
    ylab("Density") +
    ggtitle("Density of AI estimates")

ggplot(AI.bias, aes(log(Total_pre), fill=stan_bias)) +
    geom_histogram( position="dodge",binwidth=0.5) +
    geom_vline(xintercept=log(100), linetype="dashed") +
    ylab("Frequency") +
    ggtitle("Distribution of total reads pre-mapping")

mapk='ENSG00000100030'

mapk.snps=stan[Gene_id==mapk & fSNP=="yes" & Signif !="None" ,fSNP.id ]

AI[, id:=paste(POS,REF,ALT, sep=":")]

ggplot(AI[id %in% mapk.snps,], aes(x=AI_post)) +
    geom_histogram(binwidth=0.003, colour="black", fill="white") +    
    ggtitle("MAPK1 potentially biased fSNPs")



#' ## Conclusions

## Not sure about multi-SNP pair behaviour. In this case r2 is 0.2
## between distant SNPs and 1 between close ones. So I could just set AI>0.5 to 0.5
## as the AI estimate for the distant SNP is likely OK.

## I think ignoring AI>0.5 is OK.

## No gross difference between AI estimates from potentially biased vs unbiased fSNPs
## based on model results.
