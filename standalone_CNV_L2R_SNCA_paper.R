#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)

# run like dis:
# Rscript --vanilla standalone_CNV_L2R.R $sampleID

SAMPLENAME = args[1]
print(args[1])
print(SAMPLENAME)
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
# set scientific notation to off
options(scipen=20)
#remove NA's
PLOT$L2R <- sub("^$", "NA", PLOT$L2R)
PLOT$L2R <- as.numeric(PLOT$L2R)
PLOT <- na.omit(PLOT)
# start plotting
pdf(paste(SAMPLENAME,"_SNCA_REGIONAL_L2R_PLOT_NEW_LARGE.pdf",sep=""),height=4, width=15)
# !!! note depending on which gene you are interested in you need to change the xlim and xlab
plot(PLOT$BP,PLOT$L2R,pch=20, ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(87000000,95759466))
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0, col="blue", lwd = 3)
lines(lowess(PLOT$BP,PLOT$L2R,f = 0.01), col="red", lwd = 3) # lowess line (x,y)
grid()
dev.off()
pdf(paste(SAMPLENAME,"_SNCA_REGIONAL_L2R_PLOT_NEW_SMALL.pdf",sep=""),height=4, width=10)
plot(PLOT$BP,PLOT$L2R,pch=20,ylab="L2R ratio",xlab="CHR 4 basepair",xlim=c(89000000,93000000))
rect(xleft=90645250,xright = 90759447,ybottom=par("usr")[3], ytop=par("usr")[4], density=10, col = "blue")
abline(h=0, col="blue", lwd = 3)
lines(lowess(PLOT$BP,PLOT$L2R,f = 0.01), col="red", lwd = 3) # lowess line (x,y)
grid()
dev.off()
###
# DONE

