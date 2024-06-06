library(reshape2)
library(dplyr)
library(stringr)
library(ggplot2)
source("~/capsule/code/microarray_transformation_evaluation/Analysis_Functions/variance_functions.R")
source("~/capsule/code/microarray_transformation_evaluation/Analysis_Functions/variance_by_avg_functions.R")


# calculate variance in Microarray / RNAseq datasets in CRC models built with CRB/CRS
# compare version 19 (CRS pan_cancer_v1) and version 34 (CRB + crc_v6 prior)

# read fdr data (downloaded from Tableau, we namually added the "platform" column)
V19_data <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/TvsTA_V19_Sig - ALL sig.csv")
V34_data <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/TvsTA_V34_Sig - ALL sig.csv")

# per version, per platform - generate_variance_table based on the significant results only (FDR <= 0.1)
V19_variance_microarray_avg <- generate_variance_table_by_avg(V19_data, "microarray") 
V34_variance_microarray_avg <- generate_variance_table_by_avg(V34_data, "microarray") 
V19_variance_rnaseq_avg <- generate_variance_table_by_avg(V19_data, "rnaseq") 
V34_variance_rnaseq_avg <- generate_variance_table_by_avg(V34_data, "rnaseq") 


#### ~~ compare variance of microarray platforms between v19 and v34 ~~ ###
compare_variance_by_avg_plot(V19_variance_microarray_avg, "microarray V19 CRC by AVG", V34_variance_microarray_avg, "microarray V34 CRC by AVG")

#### ~~ compare variance of rnaseqs platform between v19 and v34 ~~ ###
compare_variance_by_avg_plot(V19_variance_rnaseq_avg, "rnaseq V19 CRC by AVG", V34_variance_rnaseq_avg, "rnaseq V34 CRC by AVG")

#### ~~ compare variance of microarray and rnaseq platform within v19 ~~ ###
compare_variance_by_avg_plot(V19_variance_microarray_avg, "microarray V19 CRC by AVG", V19_variance_rnaseq_avg, "rnaseq V19 CRC by AVG")

#### ~~ compare variance of microarray and rnaseq platform within v34 ~~ ###
compare_variance_by_avg_plot(V34_variance_microarray_avg, "microarray V34 CRC by AVG", V34_variance_rnaseq_avg, "rnaseq V34 CRC by AVG")
