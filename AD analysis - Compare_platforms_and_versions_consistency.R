# install.packages("reshape2")
library(reshape2)
library(dplyr)
library(stringr)
library(ggplot2)
source("~/capsule/code/microarray_transformation_evaluation/Analysis_Functions/variance_functions.R")

# calculate variance in Microarray / RNAseq datasets in CRC models built with CRB/CRS
# compare version 19 (CRS pan_cancer_v1) and version 34 (CRB + crc_v6 prior)

# read fdr data (downloaded from Tableau, we namually added the "platform" column)
V6_data <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v6.csv") # CRS
V7_data <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v7.csv") # CRB

# per version, per platform - generate_variance_table based on the significant results only (FDR <= 0.1)
V6_variance_microarray <- generate_variance_table(V6_data, "microarray") 
V7_variance_microarray <- generate_variance_table(V7_data, "microarray") 
V6_variance_rnaseq <- generate_variance_table(V6_data, "rnaseq") 
V7_variance_rnaseq <- generate_variance_table(V7_data, "rnaseq") 


#### ~~ compare variance of microarray platforms between v6 and v7 ~~ ###
compare_variance(V6_variance_microarray, "AD microarray V6 CRS", V7_variance_microarray, "AD microarray V7 CRB")

#### ~~ compare variance of rnaseqs platform between v6 and v7 ~~ ###
compare_variance(V6_variance_rnaseq, "AD rnaseq V6 CRS", V7_variance_rnaseq, "AD rnaseq V7 CRB")

#### ~~ compare variance of microarray and rnaseq platform within v7 ~~ ###
compare_variance(V7_variance_microarray, "AD microarray V7 CRB", V7_variance_rnaseq, "AD rnaseq V7 CRB")

#### ~~ compare variance of microarray and rnaseq platform within v6 ~~ ###
compare_variance(V6_variance_microarray, "AD microarray V6 CRS", V6_variance_rnaseq, "AD rnaseq V6 CRS")


# NOAM request to include all the results - so I changed the "significant" column (manual)
# so all the results will be included (so the "significant" column in the input csv are "TRUE" for all rows)

V7_data_all_sig <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v7 - ALL_sig.csv") # CRB
V7_variance_microarray_all_sig <- generate_variance_table(V7_data_all_sig, "microarray") 
V7_variance_rnaseq_all_sig <- generate_variance_table(V7_data_all_sig, "rnaseq") 
#### ~~ compare variance of microarray and rnaseq platform within v7 ~~ ###
compare_variance(V7_variance_microarray_all_sig, "AD microarray V7 CRB", V7_variance_rnaseq_all_sig, "AD rnaseq V7 CRB")
