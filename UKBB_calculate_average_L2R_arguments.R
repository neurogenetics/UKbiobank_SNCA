#!/bin/env Rscript

## DESCRIPTION
	# Author(s): Mary B. Makarious, Cornelis Blauwendraat
	# Project: UKBB CNVs - SNCA 
	# Collaborators: --
	# Date Last Updated: 23.04.2020

## GOAL
  # Generate averages for the log R ratios for each of the samples
  # Plot out these averages with a corresponding density+histogram plot 
  # In this case, doing this for a specified SNCA range and UKBB samples 
		
## PROPOSED WORKFLOW
	# 0. Getting Started
	# 1. Calculate the averages for each column 
	# 2. Transpose the dataframe 
	# 3. Sort in descending order and save
	# 4. Generating plots and saving out 

## TO USE 
  # Rscript --vanilla UKBB_L2R_arguments.R SNCA_gene_region_ukb_l2r_chr4_v2 sampleID_list
    # Argument 1: UKBB region with just log R ratios (no header, no row names)
    # Argument 2: Sample ID list as a .txt file (no header)

#########################################################################################################################
##### 0. GETTING STARTED ################################################################################################
#########################################################################################################################

# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(data.table)) install.packages('data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(gridExtra)) install.packages('gridExtra')

# Load the necessary packages 
library(tidyverse)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2) 
library(gridExtra)

# Accepting arguments
args = commandArgs(trailingOnly=TRUE)
LR2_REGION = args[1]
SAMPLE_LIST = args[2]

# Read in the L2R file
UKB_L2R <- fread(paste(args[1], ".txt", sep=""), header=FALSE, sep=" ")

# Read in the sampleIDs 
sampleID_file <- read.table(paste(args[2], ".txt", sep=""))

# Make a list of the sample IDs 
	# This is so we can add this as a header to the L2R file later
sampleIDs_list <- unlist(as.character(sampleID_file$V1))

# Add column names to the L2R file
names(UKB_L2R) = c(sampleIDs_list)


#########################################################################################################################
##### 1. CALCULATE THE AVERAGES FOR EACH COLUMN #########################################################################
#########################################################################################################################

# Calculate the means of each column
# Calculate even though there are NAs in the column (NA = variant doesn't exist but others do)

# Using full dataset
sample_l2r_averages <- colwise(mean)(UKB_L2R, na.rm = TRUE)

#########################################################################################################################
##### 2. TRANSPOSE THE AVERAGES AND SAMPLE IDS ##########################################################################
#########################################################################################################################

# Transpose into a dataframe of averages 
transposed_l2r_averages <- as.data.frame(t(sample_l2r_averages))

# Get rid of R's annoying "row names" and keep the sample ID as a column instead 
transposed_l2r_averages = rownames_to_column(transposed_l2r_averages, "SampleID")

# Give meaningful column names 
names(transposed_l2r_averages) = c("SampleID", "Average_L2R")

# See the dataframe
#head(transposed_l2r_averages)

#########################################################################################################################
##### 3. SORT IN DESCENDING ORDER AND SAVE ##############################################################################
#########################################################################################################################

# Sort averages from high to low 
sorted_l2r_averages <- transposed_l2r_averages %>% arrange(desc(Average_L2R))

# Save out the sorted averages 
write.table(sorted_l2r_averages, file = paste("sorted_", args[1], "_averageL2R.txt", sep=""), col.names = TRUE, row.names=FALSE, na="", quote = FALSE, sep="\t")

#########################################################################################################################
##### 4. GENERATING PLOTS AND SAVING OUT ################################################################################
#########################################################################################################################

# Formatting 
require(scales)

mult_format <- function() {
     function(x) format(10*x,digits = 2) 
}


# Use a scatterplot

l2r_averages <- ggplot(data = sorted_l2r_averages, aes(y = 1:nrow(sorted_l2r_averages), x=Average_L2R)) +
  geom_point(alpha=0.5, size=0.5, color="navy") +
  theme_light() +
  theme(axis.ticks = element_blank()) +
  ggtitle("Average log R Ratios") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_x_continuous(name = "Average L2R", limits = c(-0.6,0.4)) +
  scale_y_continuous(labels=comma, name = "Samples \n(Total)", limits = c(0,500000)) 
  
#l2r_averages

histogram_density <- ggplot(data = sorted_l2r_averages, aes(x=Average_L2R)) +
  theme_light() +
  theme(axis.ticks = element_blank()) +
  geom_histogram(aes(y=..density..), alpha=.5, color="orange", fill="orange") +
  geom_density(alpha=.5, fill="navy", color="blue") +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = mult_format(), name = "Samples \n(Total %)", limits = c(0,10)) 

#histogram_density

# Save out the plots 
ggsave(paste(args[1], "_L2R_AVERAGES_SCATTER_PLOT.jpg", sep=""), l2r_averages, width = 5, height = 3.5, units = "in")
ggsave(paste(args[1], "_L2R_AVERAGES_HISTDEN_PLOT.jpg", sep=""), histogram_density, width = 5, height = 3.5, units = "in")

# Use a scatterplot

#l2r_averages

# Plot both together 
both_plots <- grid.arrange(histogram_density, l2r_averages, ncol=1, nrow=2, widths=c(1), heights=c(2.5, 5))

# Save out 
ggsave(paste(args[1], "_L2R_AVERAGES_BOTH_PLOT.jpg", sep=""), both_plots, width = 5, height = 3.5, units = "in")


