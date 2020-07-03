#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# run like dis:
# Rscript --vanilla standalone_CNV.R $sampleID

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
# BAF <- fread("SNCA_gene_region1_ukb_baf_chr4_v2.txt",header=F)
# FAM <- fread("ukb33601_cal_chr4_v2_s488264.fam",header=F)
# BAF2 <- t(BAF)
# BAF3 <- cbind(FAM,BAF2)
# fwrite(BAF3, file="SNCA_gene_region_1_ukb_baf_chr4_v2_FINAL.txt", quote=FALSE,row.names=F,sep="\t")

require("data.table")
BAF <- fread("SNCA_gene_region_1_ukb_baf_chr4_v2_FINAL.txt",header=F)
newdata <- subset(BAF, V1 == SAMPLENAME)
# V1 - V6 is crap
newdata$V1 <- NULL
newdata$V2 <- NULL
newdata$V3 <- NULL
newdata$V4 <- NULL
newdata$V5 <- NULL
newdata$V6 <- NULL
BAF <- t(newdata)
BIM <- fread("SNCA_gene_region1_ukb_bim_chr4_v2.txt",header=F)
PLOT <- cbind(BAF,BIM)
names(PLOT) <- c("BAF","CHR","RS","CRAP","BP","A1","A2")

# L2R and BAF plots
options(scipen=20)
pdf(paste(SAMPLENAME,"_SNCA_REGIONAL_BAF_PLOT.pdf",sep=""),height=4, width=20)
plot(PLOT$BP,PLOT$BAF,pch=20,ylab="B allele frequency",xlab="CHR 4 basepair",xlim=c(85645250,95759466))
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0.66, col="blue")
abline(h=0.33, col="blue")
plot(PLOT$BP,PLOT$BAF,pch=20,ylab="B allele frequency",xlab="CHR 4 basepair",xlim=c(90145250,91259466))
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0.66, col="blue")
abline(h=0.33, col="blue")
plot(PLOT$BP,PLOT$BAF,pch=20,ylab="B allele frequency",xlab="CHR 4 basepair",xlim=c(90645250,90759466))
abline(h=0.66, col="blue")
abline(h=0.33, col="blue")
dev.off()

###
# DONE

