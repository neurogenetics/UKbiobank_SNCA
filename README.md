# UKbiobank SNCA analysis
SNCA mutation (missense and copy number variant) analysis in the UK biobank targeted gene approach

July 2020

Code contributors -> Cornelis, Mary, Hampton, Mike and Andy

	LNG ❤️ Open science  😍

Goals: 

1) Understanding and checking how the CNV data looks like in UKbiobank...

2) Check for gene copy number variant in SNCA gene region...

3) Check for pathnogenic missense mutations in SNCA...

### Brief summary:
coming soon...

## Structure of Repo:
1. [Downloading data](#1-Downloading-data)
2. [UKbiobank information data structure](#2-UKbiobank-information-data-structure)
3. [Check files if the dimensions make sense](#3-Check-files-if-the-dimensions-make-sense)
4. [Diving into regions of interest](#4-Diving-into-regions-of-interest)

      A. [SNCA](#SNCA-region-subset)
     
      B. [Prioritization of sampleIDs](#Prioritization-of-sampleIDs)

5. [Standalone CNV plotting scripts](#5-Standalone-CNV-plotting-scripts)
6. [Results summary](#6-Results-summary)

      A. [SNCA results](#SNCA-results)
      
7. [Using Exome data as potential replication?](#7-Using-Exome-data-as-potential-replication)
8. [ICD10 code explorer of CNV carriers](#8-ICD10-code-explorer-of-CNV-carriers)
9. [Relatedness of CNV carriers](#9-Relatedness-of-CNV-carriers)
10. [Additional things to check](#10-Additional-things-to-check)


### 1. Downloading data
```
For sure not small files...
# size of folders...
1.4T	BAF/
2.3T	L2R/

module load ukbb/0.1
for CHROMOSOME_NUMBER in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
ukbgene l2r -c"$CHROMOSOME_NUMBER" -a/PATH/TO/KEY/KEYFILE.key 
ukbgene baf -c"$CHROMOSOME_NUMBER" -a/PATH/TO/KEY/KEYFILE.key 
done
```

### 2. UKbiobank information data structure
```
CNV files
===
These files contain the B-Allele-Frequency (baf) and Log2Ratio (log2r) 
transformed intensitiy values for performing CNV calling. 
There is a separate file for baf and log2r per chromosome. 
These are plaintext files with space separated columns. 
The rows corresond to markers (ordered as the calls BIM file) 
and the columns correspond to samples (ordered as the calls FAM file)
Missing values are represented by -1.

------------------> samples => 'ukbiobank_file'.fam
|
|
|
|
|
|
v
markers => .bim files
```


### 3. Check files if the dimensions make sense
```
------------------> samples => 'ukbiobank_file'.fam

head -1 ukb_l2r_chr14_v2.txt | tr ' ' '\n' | wc -l
head -1 ukb_baf_chr1_v2.txt | tr ' ' '\n' | wc -l
# 488377 samples corresponds to fam file

##### check files if the dimensions make sense
|
|
|
|
|
|
v
markers => .bim files

# number of markers is for each file different...
Markers	CHR
63487	1
61966	2
52300	3
47443	4
46314	5
53695	6
42722	7
38591	8
34310	9
38308	10
40824	11
37302	12
26806	13
25509	14
24467	15
28960	16
28835	17
21962	18
26186	19
19959	20
11342	21
12968	22
265	MT
18857	X
1357	XY
691	Y
```

### 4. Diving into regions of interest
Little bit if background... SNCA multiplications are associated Parkinson disease (PMID:14593171). Penetrance is often very high with these mutations... but large scale investigation of these genes hasnt been performed yet. 

#### SNCA region subset
```
# SNCA region 
1) SNCA gene +/5 Mb(chr4:85645250-95759466, hg19)
2) SNCA gene +/-0.5 Mb (chr4:90145250-91259466, hg19) 
3) SNCA gene body (chr4:90645250-90759466, hg19)

# check bim
cd /PATH/TO/UKBIOBANK/CNV_UKB/BIM
head ukb_snp_chr4_v2.bim
awk '$4 > 85645250' ukb_snp_chr4_v2.bim | awk '$4 < 95759466' | wc -l
# 2373
awk '$4 > 90145250' ukb_snp_chr4_v2.bim | awk '$4 < 91259466' | wc -l
# 391
awk '$4 > 90645250' ukb_snp_chr4_v2.bim | awk '$4 < 90759466' | wc -l
# 44
# save variants
awk '$4 > 85645250' ukb_snp_chr4_v2.bim | awk '$4 < 95759466' > SNCA_region1_variants.txt
awk '$4 > 90145250' ukb_snp_chr4_v2.bim | awk '$4 < 91259466' > SNCA_region2_variants.txt
awk '$4 > 90645250' ukb_snp_chr4_v2.bim | awk '$4 < 90759466' > SNCA_region3_variants.txt
scp *_variants.txt ../

# subsetting BAF and L2R
cd /PATH/TO/UKBIOBANK/CNV_UKB/
# region 1 = 2373 variants
# from rs62303078 = 21814 # position 1
# to rs79456659 = 24186 # position 2
sed -n 21814,24186p L2R/ukb_l2r_chr4_v2.txt > SNCA_gene_region1_ukb_l2r_chr4_v2.txt
sed -n 21814,24186p BAF/ukb_baf_chr4_v2.txt > SNCA_gene_region1_ukb_baf_chr4_v2.txt
# region 2 = 391 variants
# from rs73845905 = 22961 # position 1
# to rs79190057 = 23351 # position 2
sed -n 22961,23351p L2R/ukb_l2r_chr4_v2.txt > SNCA_gene_region2_ukb_l2r_chr4_v2.txt
sed -n 22961,23351p BAF/ukb_baf_chr4_v2.txt > SNCA_gene_region2_ukb_baf_chr4_v2.txt
# region 3 = 44 variants
# from rs17180453 = 23164 # position 1
# to rs1372519 = 23207 # position 2
sed -n 23164,23207p L2R/ukb_l2r_chr4_v2.txt > SNCA_gene_region3_ukb_l2r_chr4_v2.txt
sed -n 23164,23207p BAF/ukb_baf_chr4_v2.txt > SNCA_gene_region3_ukb_baf_chr4_v2.txt
```

#### Extract BAF and L2R from large data-frames and prioritize based on these values for plotting

Brief explanation:
BAF (= B allele frequency) values, when this value is 0 this mean genotype AA, when 0.5 this means genotype AB and when this is 1 this means genotype BB. If there is a duplication (which then means three genotypes) for example AAB then this value is 0.33 and for ABB 0.66 and for triplication (which then means four genotypes) for example AAAB then this value is 0.25 and for ABBB 0.75. 
Note that deletions typically do not have any 0.5 values because all A genotypes will be 0 and all B genotypes will 1, however this should not be confused with homozygous streches of DNA were all AA genotypes will be 0 and all B genotypes will 1.

L2R (=log R ratio) values, when this value is 0 this mean genotype AA, AB or BB (= normal). When this value is lower than 0 this indicates a deletion so genotype A or B eg -0.45. When this value is higher than 0 this indicates a duplication so genotype for example genotype ABB, AAB with values of 0.3 and higher, similar for triplication with values of 0.75 and higher.

See Table 2 of [here](https://www.illumina.com/documents/products/technotes/technote_cnv_algorithms.pdf) for more explanations 

So below we are prioritizing:

BAF values of => <0.85 and >0.65 + <0.35 and >0.15

L2R values of => <-0.20 and >0.20 <= note these are example values and depends on your region of interest, typically general outlier detection like couple SD from mean should work.

#### L2R Script 


Script can be found [here](https://github.com/neurogenetics/UKbiobank_SNCA/blob/master/UKBB_calculate_average_L2R_arguments.R)

name => UKBB_calculate_average_L2R_arguments.R

```
## GOAL
  # Generate averages for the log R ratios for each of the samples
  # Plot out these averages with a corresponding density+histogram plot 
  # In this case, doing this for a specified SNCA range and UKBB samples 
		
## WORKFLOW
	# 0. Getting Started
	# 1. Calculate the averages for each column 
	# 2. Transpose the dataframe 
	# 3. Sort in descending order and save
	# 4. Generating plots and saving out 

## TO USE 
  # Rscript --vanilla UKBB_L2R_arguments.R $SNCA_gene_region_ukb_l2r_chr4_v2 $sampleID_list
    # Argument 1: UKBB region with just log R ratios (no header, no row names)
    # Argument 2: UKBB Sample ID list as a .txt file (no header)
    
## OUTPUT
    # A text file with all average values per input region => *_averageL2R.txt
    # A couple images 

```

#### BAF Script

Script can be found here
names => UKBB_calculate_BAF_arguments.R

```
## GOAL
  # Generate counts of variants within a specified range for each sample 
  # In this case, counting variants >=0.65 and <=0.85 [OR] >=0.15 and <=0.35 for UKBB samples 

## WORKFLOW
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

## OUTPUT
    # A text file with the number of BAF values of interest per input region => *BAFcounts_bw015_035and065_085.txt
 
```

#### Prioritization of sampleIDs

Now we have 6 files per gene:

```
#L2R:
sorted_SNCA_gene_region1_ukb_l2r_chr4_v2_averageL2R.txt
sorted_SNCA_gene_region2_ukb_l2r_chr4_v2_averageL2R.txt
sorted_SNCA_gene_region3_ukb_l2r_chr4_v2_averageL2R.txt

#BAF:
sorted_SNCA_gene_region1_ukb_baf_chr4_v2_BAFcounts_bw065_085and015_035.txt
sorted_SNCA_gene_region2_ukb_baf_chr4_v2_BAFcounts_bw065_085and015_035.txt
sorted_SNCA_gene_region3_ukb_baf_chr4_v2_BAFcounts_bw065_085and015_035.txt

```

##### L2R

Lets start with L2R... We are looking for the very high and very low values.

Given that we dont expect a very high number of deletions and duplications we selected of each of the three above described regions:

- the 25 highest LR2 values 
- the 25 lowest LR2 values 


##### BAF

Lets continue with BAF... We are looking for the number of values <0.85 and >0.65 + <0.35 and >0.15.

Again given that we dont expect a very high number of deletions and duplications we selected

- samples with >=25 variants in region 1 (=SNCA gene +/5 Mb)
- samples with >=10 variants in region 2 (SNCA gene +/-0.5 Mb)
- samples with >=5 variants in region 3 (SNCA gene body)

These sample IDs we be plotted in step 5 for closer inspection.


### 5. Standalone CNV plotting scripts

#### Things to plot for SNCA region

```
cut -f 1 HIGH_BAF_variants.txt | grep -v FID > input_standalone_BAF.txt
cut -f 1 LOW_HIGH_L2R.txt | grep -v SampleID > input_standalone_L2R.txt

# example option how to run:
Rscript --vanilla standalone_CNV_BAF.R UKB_biobank_sample_ID
Rscript --vanilla standalone_CNV_L2R.R UKB_biobank_sample_ID

'Plotting BAF'
#!/bin/sh
# sbatch --cpus-per-task=10 --mem=40g --mail-type=END --time=24:00:00 input_standalone_BAF.sh
module load R
cat samples_of_interest.txt | while read line
do 
	Rscript --vanilla standalone_CNV_BAF.R $line
done
   
'Plotting L2R'
#!/bin/sh
# sbatch --cpus-per-task=10 --mem=40g --mail-type=END --time=24:00:00 input_standalone_L2R.sh
module load R
cat samples_of_interest.txt | while read line
do 
	Rscript --vanilla standalone_CNV_L2R.R $line
done
```


#### How does standalone_CNV_BAF.R and standalone_CNV_L2R.R looks like under the hood

##### standalone_CNV_BAF.R

Script can be found here
names => standalone_CNV_BAF.R 

add in pdf example

If you want to make it a bit fancier 
Script can be found here
names => standalone_CNV_BAF_SNCA_paper.R 

add in pdf example


##### standalone_CNV_L2R.R

Script can be found here
names => standalone_CNV_L2R.R 

add in pdf example

If you want to make it a bit fancier 
Script can be found here
names => standalone_CNV_L2R_SNCA_paper.R 

add in pdf example

### 6. Results summary

#### SNCA results:
Manual inspection of plots resulted in:

| Event type  | Count |
| ------------- | ------------- |
| duplications  | 6  |
| deletions | 6  |
| complex/mosaicism? | 14  |

#### Bit deeper dive into these...

What does each column mean again?
```
1 => SNCA gene +/5 Mb(chr4:85645250-95759466, hg19)
2 => SNCA gene +/-0.5 Mb (chr4:90145250-91259466, hg19) 
3 => SNCA gene body (chr4:90645250-90759466, hg19)

For BAF it is the number of variants with value <0.85 and >0.65 OR <0.35 and >0.15
For L2R it is the average L2R value in that region
```

#### Inspect wider region of "complex" hits

Script can be found here
names => CNV_closer_inspection.R 
Basically very similar as standalone_CNV_BAF.R and standalone_CNV_L2R.R but now assesses full chr4 region and wider SNCA region +/- 20mb

add in pdf example

If you want to make it a bit fancier 
Script can be found here
names => standalone_CNV_BAF_SNCA_paper.R 

add in pdf example

#### Inspect all chromosomes of "complex" hits

Script can be found here
names => full_chromome_plot.R 
Basically very similar as standalone_CNV_BAF.R and standalone_CNV_L2R.R but now assesses full all autosomes


#### OK coolio so now whats next?

- Try to find replication of CNVs in exome data

- Check ICD10 codes of CNVs carriers if anything stands out

- Check relatedness of CNV carriers


### 7. Using Exome data as potential replication?

One way of potentially "replicate" the findings of CNV is using the exome sequencing data. Specifically the allele depth field since you would expect an imbalance there.

Three samples of interest have exome data, besides these we will use 10 random other samples as "controls".

#### Downloading UKB exome data

```
# Getting 10 random samples from .fam file
sort -R ukb33601_efe_chr1_v1_s49959.fam | head -n 10
```

```
cd /PATH/TO/UKBIOBANK/CNV_UKB/EXOME
module load ukbb/0.1

cat UKB_ID_OF_INTEREST.txt  | while read UKBID
do 
	ukbfetch -e$UKBID -d23161_0_0 -a/PATH/TO/UKBIOBANK/UKBB_key_and_programs/application_number.key
done
```

#### Check allelic depth fields in exome data

```
module load vcftools
# reformat files to keep only chr4, only heterozygous variants and have it in the right format.

cat UKB_ID_OF_INTEREST.txt  | while read UKBID
do 
	vcftools --gzvcf "$UKBID"_23161_0_0.gvcf.gz --chr 4 --recode --recode-INFO-all --out "$UKBID"_CHR4_region
	grep "DP=" "$UKBID"_CHR4_region.recode.vcf | cut -f 1,2,9,10 | sed 's/:/\t/g' | sed 's/,/\t/g' | cut -f 1,2,9,10,11,13,14 | grep -v "1/1" | grep "0/1" > "$UKBID"_chr4.txt
done
```


```
# now plot in R
module load R

cat UKB_ID_OF_INTEREST.txt  | while read UKBID
do 
	Rscript --vanilla plot_AD_chr4.R $UKBID
done
```


```
How does plot_AD_chr4.R looks like under the hood?

#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# run like dis:
# Rscript --vanilla plot_AD_chr4.R $sampleID
SAMPLENAME = args[1]
print(args[1])
print(SAMPLENAME)
library(ggplot2)
data <- read.table(paste(SAMPLENAME,"_chr4.txt",sep=""),header=F)
colnames(data) <- c("CHR","BP","GT","AD1","AD2","DP","GQ")
# keep only high quality variants
newdata <- subset(data, DP > 19 & GQ > 80)
newdata$ADratio <- as.numeric(newdata$AD1)/as.numeric(newdata$AD2)
newdata$ADratio2 <- as.numeric(newdata$AD2)/as.numeric(newdata$AD1)
# set all value to >1
newdata$ADratio3 <- pmax(newdata$ADratio, newdata$ADratio2)
# remove extreme AD ratio's..
newdata2 <- subset(newdata, ADratio3 < 4)
# ggplot style # note SNCA coordinates are HG38
options(scipen=20)
newdata2$ADratio <- newdata2$ADratio3
p <- ggplot(newdata2, aes(BP,ADratio)) + geom_point() + geom_smooth()
p + annotate("rect", xmin=89724099, xmax=89838315, ymin=1, ymax= 4, 
             fill=NA, colour="red")
ggsave(paste(SAMPLENAME,"_chr4_new_v2.pdf",sep=""), width = 20, height = 4)

```


### 8. ICD10 code explorer of CNV carriers

#### Investigating ICD10 codes per sample of interest:

```
cd /path/to/UKBIOBANK/CNV_UKB/
mkdir phenotypes
# check in known disease group
grep -f SNCA_samples_moving_forward.txt /path/to/UKBIOBANK/PHENOTYPE_DATA/disease_groups/* > phenotypes/disease_groups.txt
# check in ICD10 codes
grep -f SNCA_samples_moving_forward.txt /path/to/UKBIOBANK/ICD10_UKBB/ICD10_codes/massive_ICD10_ALL_table.txt | sort -nk1 > phenotypes/ICD10.txt
# check European?
grep -f SNCA_samples_moving_forward.txt /path/to/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_to_use.txt > phenotypes/covariates.txt
# add in ICD10 annotation
module load R
R
ICD10 <- read.table("/path/to/UKBIOBANK/ICD10_UKBB/ICD10_codes/ICD10_coding.tsv",header=T)
pheno <- read.table("phenotypes/ICD10.txt",header=F)
MM <- merge(pheno,ICD10,by.x="V2",by.y="coding")
write.table(MM, file = "phenotypes/ICD10_annotated.txt",row.names=FALSE, quote = FALSE, sep="\t")
q()
n
# check for PD-ish ICD10 codes....
```

### 9. Relatedness of CNV carriers
#### Checking if any if the CNV carriers are related

```
module load plink
module load R
module load flashpca

cd /path/to/UKBIOBANK/GENOTYPE_DATA

# extract data from raw genotyes
for chnum in {1..22};
  do
  ./plink2 --bed ukb_cal_chr"$chnum"_v2.bed --bim ukb_snp_chr"$chnum"_v2.bim --fam ukb33601_cal_chr1_v2_s488363.fam \
  --keep /path/to/UKBIOBANK/CNV_UKB/SNCA_samples_moving_forward_plink.txt --make-bed \
  --out /path/to/UKBIOBANK/CNV_UKB/genotypes/SNCA_dups_"$chnum"
done

# merge data
cd /path/to/UKBIOBANK/CNV_UKB/genotypes/

ls | grep "bim" | sed 's/.bim//g' > merge_list.txt
module load plink
plink --merge-list merge_list.txt --make-bed --out MERGED

plink --bfile MERGED --maf 0.05 --geno 0.05 --hwe 5e-6 --make-bed --out MERGED_FILTER
plink --bfile MERGED_FILTER --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
plink --bfile MERGED_FILTER --extract pruned_data.prune.in --make-bed --out MERGED_FILTER_PRUNE
plink --bfile MERGED_FILTER_PRUNE --genome --min 0.05 --out MERGED_FILTER_PRUNE

# investigate .genome file for high PIHAT values ~0.5 means first degree, ~0.25 means second degree, etc
sort -gk 10 MERGED_FILTER_PRUNE.genome | head

# most noteworthy:
```

| Sample1  | Sample2 | PIHAT |
| ------------- | ------------- | ------------- |
| Dup1 | Dup6  | 0.4730 |
| Del5 | Del2  | 0.5637 |

```
# multiple others ~PIHAT of ~0.125

# also testing KING

module load king
king -b MERGED_FILTER_PRUNE.bed --related
# same result 2 first degrees...

```

### 10. Additional things to check

#### are SNCA mutations on the UKbiobank array (raw data)?
```
Check for SNCA mutations in Clinvar June 23 2020 => https://www.ncbi.nlm.nih.gov/clinvar => SNCA[gene]
Protein change variants only

# for sure pathogenic (hg19)
grep 90756731 ukb_snp_chr4_v2.bim # Ala30Pro
grep 90749321 ukb_snp_chr4_v2.bim # Glu46Lys
grep 90749305 ukb_snp_chr4_v2.bim # Gly51Asp
grep 90749299 ukb_snp_chr4_v2.bim # Ala53Glu
grep 90749300 ukb_snp_chr4_v2.bim # Ala53Thr

# debatable pathogenic (hg19)
grep 90749307 ukb_snp_chr4_v2.bim # His50Gln
grep 90650353 ukb_snp_chr4_v2.bim # Pro128Thr
grep 90650386 ukb_snp_chr4_v2.bim # Pro117Thr
grep 90743416 ukb_snp_chr4_v2.bim # Lys96Arg
grep 90743452 ukb_snp_chr4_v2.bim # Gly84Val
grep 90756775 ukb_snp_chr4_v2.bim # Val15Ala

# Nope none are on there...
```
#### are SNCA mutations identified in the UKbiobank exome sequencing data?
```
# for sure pathogenic (hg38
grep 89835580 ukb_fe_exm_chrall_v1.bim # Ala30Pro
grep 89828170 ukb_fe_exm_chrall_v1.bim # Glu46Lys
grep 89828154 ukb_fe_exm_chrall_v1.bim # Gly51Asp
grep 89828148 ukb_fe_exm_chrall_v1.bim # Ala53Glu
grep 89828149 ukb_fe_exm_chrall_v1.bim # Ala53Thr
# nope none are identified there

# debatable pathogenic (hg38)
grep 89828156 ukb_fe_exm_chrall_v1.bim # His50Gln => yes => 4:89828156:A:C
grep 89729202 ukb_fe_exm_chrall_v1.bim # Pro128Thr => yes => 4:89729202:G:T
grep 89729235 ukb_fe_exm_chrall_v1.bim # Pro117Thr => yes => 4:89729235:G:T
grep 89822265 ukb_fe_exm_chrall_v1.bim # Lys96Arg => yes => 4:89822265:T:C
grep 89822301 ukb_fe_exm_chrall_v1.bim # Gly84Val
grep 89835624 ukb_fe_exm_chrall_v1.bim # Val15Ala

# check number of alleles
module load plink
plink --bfile --bed ukb_efe_chr1_v1.bed --bim ukb_fe_exm_chrall_v1.bim \
--fam ukb33601_efe_chr1_v1_s49959.fam --snps 4:89828156:A:C,4:89729202:G:T,4:89729235:G:T,4:89822265:T:C \
--freq --out freq_SNCA --recodeA

 CHR              SNP   A1   A2          MAF  NCHROBS
   4   4:89729202:G:T    T    G    7.006e-05    99920 => 7 alleles => all heterozygous
   4   4:89729235:G:T    T    G    0.0001301    99920 => 13 alleles => all heterozygous
   4   4:89822265:T:C    C    T    1.001e-05    99920 => 1 allele => all heterozygous
   4   4:89828156:A:C    C    A    0.0002802    99920 => 28 alleles => all heterozygous

# any PD-ish phenotype for these carriers?
# checking PD parent or PD status 
sort -nk 7 freq_SNCA.raw | tail -8 | head -7 | cut -d " " -f 1 > SNCA_carriers.txt
sort -nk 8 freq_SNCA.raw | tail -14 | head -13 | cut -d " " -f 1 >> SNCA_carriers.txt
sort -nk 9 freq_SNCA.raw | tail -2 | head -1 | cut -d " " -f 1 >> SNCA_carriers.txt
sort -nk 10 freq_SNCA.raw | tail -29 | head -28 | cut -d " " -f 1 >> SNCA_carriers.txt

wc -l SNCA_carriers.txt # 49 so makese sense => 7+23+1+28

/path/to/UKBIOBANK/EXOME_DATA/PLINK_files/SNCA_carriers.txt

Three have a parent with PD:
sampleID1 => 4:89729235:G:T => Pro117Thr
sampleID2 => 4:89828156:A:C => His50Gln
sampleID3 => 4:89828156:A:C => His50Gln

```

## Done

![myImage](https://media.giphy.com/media/XRB1uf2F9bGOA/giphy.gif)

