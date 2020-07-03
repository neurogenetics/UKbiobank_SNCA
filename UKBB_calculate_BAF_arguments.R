#!/bin/env Rscript

## DESCRIPTION
  # Author(s): Mary B. Makarious, Cornelis Blauwendraat
  # Project: UKBB CNVs - SNCA 
  # Collaborators: --
  # Date Last Updated: 27.04.2020

## GOAL
  # Generate counts of variants within a specified range for each sample 
  # In this case, counting variants >=0.65 and <=0.85 [OR] >=0.15 and <=0.35 for UKBB samples 

## PROPOSED WORKFLOW
  # 0. Getting Started
  # 1. Replace variants <0.15 and >0.85 with NA
  # 2. Replace variants between 0.35 and 0.65 with NA 
  # 3. Transpose the dataset
  # 4. Return counts of variants within a specified range for each sample 
  # 5. Sort in descending order
  # 6. Save out file with IDs and counts

## TO USE 
  # Rscript --vanilla UKBB_BAF_arguments.R SNCA_gene_region_ukb_baf_chr4_v2 sampleID_list
    # Argument 1: UKBB region with just B allele frequencies (no header, no row names)
    # Argument 2: Sample ID list as a .txt file (no header)

#########################################################################################################################
##### 0. GETTING STARTED ################################################################################################
#########################################################################################################################

# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(plyr)) install.packages('purrr')

# Load the necessary packages 
library(tidyverse)
library(data.table)
library(dplyr)
library(plyr)
library(purrr)

# Accepting arguments
args = commandArgs(trailingOnly=TRUE)
LR2_REGION = args[1]
SAMPLE_LIST = args[2]

# Read in BAF file 
UKB_BAF <- fread(paste(args[1], ".txt", sep=""), header=FALSE, sep=" ")

# Read in the sampleIDs 
sampleID_file <- read.table(paste(args[2], ".txt", sep=""))

# Make a list of the sample IDs 
  # This is so we can add this as a header to the L2R file later
sampleIDs_list <- unlist(as.character(sampleID_file$V1))

# Add column names to the L2R file
names(UKB_BAF) = c(sampleIDs_list)

#########################################################################################################################
##### 1. REPLACE VARIANTS <0.15 AND >0.85 WITH NA #######################################################################
#########################################################################################################################

# Replace all values <0.15 or >0.85 with NA 
sorted_bafs <- UKB_BAF %>% mutate_all(funs(ifelse(.<0.15 | .>0.85, NA, .)))

#########################################################################################################################
##### 2. REPLACE VARIANTS BETWEEN 0.35 AND 0.65 WITH NA #################################################################
#########################################################################################################################

# Then replace all values between 0.35 and 0.65 with NA 
sorted_bafs2 <- sorted_bafs %>% mutate_all(funs(ifelse(.>0.35 & .<0.65, NA, .)))

#########################################################################################################################
##### 3. TRANSPOSE THE DATASET ##########################################################################################
#########################################################################################################################

# Transpose the dataset
transposed <- as.data.frame(t(sorted_bafs2))

# Get rid of R's annoying "row names" and keep the sample ID as the first column instead 
transposed = rownames_to_column(transposed, "SampleID")

#########################################################################################################################
##### 4. RETURN COUNTS OF VARIANTS WITHIN A SPECIFIED RANGE FOR EACH SAMPLE #############################################
#########################################################################################################################

# Make a function that will return a vector of counts per row
nonNAs <- function(x) {
    as.vector(apply(x, 2, function(x) length(which(!is.na(x)))))
}

# Create a column of the counts 
transposed$COUNTS <- nonNAs(sorted_bafs2)

# Keep the 2 columns we are interested in
counts_per_sample <- transposed %>% select("SampleID", "COUNTS")

#########################################################################################################################
##### 5. SORT IN DESCENDING ORDER #######################################################################################
#########################################################################################################################

# Sort by descending order 
sorted_counts_per_sample <- counts_per_sample %>% arrange(desc(COUNTS))

#########################################################################################################################
##### 6. SAVE OUT FILE WITH IDS AND COUNTS ##############################################################################
#########################################################################################################################

# Save out the file keeping counts of samples with variants between >=0.65 and <=0.85 [OR] >=0.15 and <=0.35
write.table(sorted_counts_per_sample, file = paste("sorted_", args[1], "_BAFcounts_bw065_085and015_035.txt", sep=""), col.names = TRUE, row.names=FALSE, na="", quote = FALSE, sep="\t")


