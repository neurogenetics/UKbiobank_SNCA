#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# run like dis:
# Rscript --vanilla CNV_closer_inspection.R $sampleID
# or via loop of course
#cat SNCA_samples_moving_forward.txt.txt | while read line
#do 
#	Rscript --vanilla CNV_closer_inspection.R $line
#done
#
SAMPLENAME = args[1]
print(args[1])
print(SAMPLENAME)
### prep in terminal
#1) SNCA gene +/20 Mb(chr4:70645250-110759466, hg19)
#sed -n 18206,27440p L2R/ukb_l2r_chr4_v2.txt > SNCA_gene_region20MB_ukb_l2r_chr4_v2.txt
#sed -n 18206,27440p BAF/ukb_baf_chr4_v2.txt > SNCA_gene_region20MB_ukb_baf_chr4_v2.txt
#sed -n 18206,27440p BIM/ukb_snp_chr4_v2.bim > SNCA_gene_region20MB_ukb_bim_chr4_v2.txt
# sampleIDs extracted from order of the fam file... eg 7725 means line no. 7725 from fam file
#cut -d " " -f $LINE_NUMBERS_OF_INTEREST SNCA_gene_region20MB_ukb_baf_chr4_v2.txt > SNCA_gene_region20MB_ukb_baf_chr4_v2_SOI.txt
#cut -d " " -f $LINE_NUMBERS_OF_INTEREST SNCA_gene_region20MB_ukb_l2r_chr4_v2.txt > SNCA_gene_region20MB_ukb_l2r_chr4_v2_SOI.txt
#2) full chromosome 4 region....
#cut -d " " -f $LINE_NUMBERS_OF_INTEREST BAF/ukb_baf_chr4_v2.txt > FULL_CHR4_ukb_baf_chr4_v2_SOI.txt
#cut -d " " -f $LINE_NUMBERS_OF_INTEREST L2R/ukb_l2r_chr4_v2.txt > FULL_CHR4_ukb_l2r_chr4_v2_SOI.txt
#module load R
#R
require("data.table")
# BAF 20Mb +/- SNCA
BAF <- fread("SNCA_gene_region20MB_ukb_baf_chr4_v2_SOI.txt",header=F)
BIM <- fread("SNCA_gene_region20MB_ukb_bim_chr4_v2.txt",header=F)
BIM$V1 <- NULL
BIM$V2 <- NULL
BIM$V3 <- NULL
BIM$V5 <- NULL
BIM$V6 <- NULL
PLOT_BAF_small <- cbind(BIM,BAF)
colnames(PLOT_BAF_small)<- c("BP",$SAMPLENAMES_OF_INTEREST)
# BAF full chr 4
BAF <- fread("FULL_CHR4_ukb_baf_chr4_v2_SOI.txt",header=F)
BIM <- fread("BIM/ukb_snp_chr4_v2.bim",header=F)
BIM$V1 <- NULL
BIM$V2 <- NULL
BIM$V3 <- NULL
BIM$V5 <- NULL
BIM$V6 <- NULL
PLOT_BAF_large <- cbind(BIM,BAF)
colnames(PLOT_BAF_large)<- c("BP",$SAMPLENAMES_OF_INTEREST)
# L2R 20Mb +/- SNCA
L2R <- fread("SNCA_gene_region20MB_ukb_l2r_chr4_v2_SOI.txt",header=F)
BIM <- fread("SNCA_gene_region20MB_ukb_bim_chr4_v2.txt",header=F)
BIM$V1 <- NULL
BIM$V2 <- NULL
BIM$V3 <- NULL
BIM$V5 <- NULL
BIM$V6 <- NULL
PLOT_L2R_small <- cbind(BIM,L2R)
colnames(PLOT_L2R_small)<- c("BP",$SAMPLENAMES_OF_INTEREST)
# L2R full chr 4
L2R <- fread("FULL_CHR4_ukb_l2r_chr4_v2_SOI.txt",header=F)
BIM <- fread("BIM/ukb_snp_chr4_v2.bim",header=F)
BIM$V1 <- NULL
BIM$V2 <- NULL
BIM$V3 <- NULL
BIM$V5 <- NULL
BIM$V6 <- NULL
PLOT_L2R_large <- cbind(BIM,L2R)
colnames(PLOT_L2R_large)<- c("BP",$SAMPLENAMES_OF_INTEREST)
# BAF and L2R plots
options(scipen=20)
pdf(paste(SAMPLENAME,"_LARGER_SNCA_PLOT.pdf",sep=""),height=4, width=20)
# figure 1
plot(PLOT_BAF_large$BP,PLOT_BAF_large[[SAMPLENAME]],pch=20,ylab="B allele frequency",xlab="CHR 4 basepair",xlim=c(1,191154276),main="Full chromosome 4 BAF")
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0.66, col="blue")
abline(h=0.33, col="blue")
abline(h=0.5, col="blue", lwd = 3)
grid()
# figure 2
plot(PLOT_BAF_small$BP,PLOT_BAF_small[[SAMPLENAME]],pch=20,ylab="B allele frequency",xlab="CHR 4 basepair",xlim=c(70645250,110759466),main="SNCA region +/- 20Mb BAF")
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0.66, col="blue")
abline(h=0.33, col="blue")
abline(h=0.5, col="blue", lwd = 3)
grid()
# figure 3
plot(PLOT_L2R_large$BP,PLOT_L2R_large[[SAMPLENAME]],pch=20,ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(1,191154276),main="Full chromosome 4 L2R")
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0, col="blue", lwd = 3)
lines(lowess(PLOT_L2R_large$BP,PLOT_L2R_large[[SAMPLENAME]],f = 0.0001), col="red", lwd = 3)
grid()
# figure 4
plot(PLOT_L2R_small$BP,PLOT_L2R_small[[SAMPLENAME]],pch=20,ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(70645250,110759466),main="SNCA region +/- 20Mb BAF")
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0, col="blue", lwd = 3)
lines(lowess(PLOT_L2R_small$BP,PLOT_L2R_small[[SAMPLENAME]],f = 0.0001), col="red", lwd = 3)
grid()
dev.off()
# ALL DONE