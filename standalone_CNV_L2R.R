#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# run like dis:
# Rscript --vanilla standalone_CNV_L2R.R $sampleID

SAMPLENAME = args[1]
print(args[1])
print(SAMPLENAME)

###### 
# stand alone CNV plotter
# previously done:
# subset data # 2684 snps only (wider region)
# sed -n 21814,24186p L2R/ukb_l2r_chr4_v2.txt > SNCA_gene_region1_ukb_l2r_chr4_v2.txt
# sed -n 21814,24186p BAF/ukb_baf_chr4_v2.txt > SNCA_gene_region1_ukb_baf_chr4_v2.txt
# sed -n 21814,24186p BIM/ukb_snp_chr4_v2.bim > SNCA_gene_region1_ukb_bim_chr4_v2.txt
# R
# require("data.table")
# L2R <- fread("SNCA_gene_region1_ukb_l2r_chr4_v2.txt",header=F)
# FAM <- fread("ukb33601_cal_chr4_v2_s488264.fam",header=F)
# L2R2 <- t(L2R)
# L2R3 <- cbind(FAM,L2R2)
# fwrite(L2R3, file="SNCA_gene_region_1_ukb_l2r_chr4_v2_FINAL.txt", quote=FALSE,row.names=F,sep="\t")

require("data.table")
L2R <- fread("SNCA_gene_region_1_ukb_l2r_chr4_v2_FINAL.txt",header=F)
newdata <- subset(L2R, V1 == SAMPLENAME)
# V1 - V6 is crap
newdata$V1 <- NULL
newdata$V2 <- NULL
newdata$V3 <- NULL
newdata$V4 <- NULL
newdata$V5 <- NULL
newdata$V6 <- NULL
L2R <- t(newdata)
BIM <- fread("SNCA_gene_region1_ukb_bim_chr4_v2.txt",header=F)
PLOT <- cbind(L2R,BIM)
names(PLOT) <- c("L2R","CHR","RS","CRAP","BP","A1","A2")
options(scipen=20)
pdf(paste(SAMPLENAME,"_SNCA_REGIONAL_L2R_PLOT.pdf",sep=""),height=4, width=20)
# !!! note depending on which gene you are interested in you need to change the xlim and xlab
plot(PLOT$BP,PLOT$L2R,pch=20,ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(85645250,95759466))
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0, col="blue")
plot(PLOT$BP,PLOT$L2R,pch=20,ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(90145250,91259466))
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0, col="blue")
plot(PLOT$BP,PLOT$L2R,pch=20,ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(90645250,90759466))
abline(h=0, col="blue")
dev.off()

###
# DONE
