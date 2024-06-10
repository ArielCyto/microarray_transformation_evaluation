library(reshape2)
library(dplyr)
library(stringr)
library(ggplot2)

source("~/capsule/code/microarray_transformation_evaluation/Analysis_Functions/variance_functions.R")
source("~/capsule/code/microarray_transformation_evaluation/Analysis_Functions/variance_by_avg_functions.R")

# COMPARING WITHIN CRB VERSION - rnaseq vs microarray
# read the data into 2 files - one for platform microarray and one for rnaseq platform
ad_v7 = read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v7.csv")
v7_avg_microarray_all_sig <- generate_variance_table_by_avg(ad_v7,"microarray")
v7_avg_rnaseq_all_sig <- generate_variance_table_by_avg(ad_v7, "rnaseq")
# compare both files
comperison_by_avg <- compare_variance_by_avg_plot(v7_avg_microarray_all_sig,"AD CRB avg microarray", v7_avg_rnaseq_all_sig, "AD CRB avg rnaseq")


# COMPARING WITHIN CRS VERSION - rnaseq vs microarray
# read the data into 2 files - one for platform microarray and one for rnaseq platform
v6_avg_microarray_all_sig <- generate_variance_table_by_avg(read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v6.csv"),"microarray")
v6_avg_rnaseq_all_sig <- generate_variance_table_by_avg(read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v6.csv"), "rnaseq")
# compare both files
comperison_by_avg <- compare_variance_by_avg_plot(v6_avg_microarray_all_sig,"AD CRS avg microarray", v6_avg_rnaseq_all_sig, "AD CRS avg rnaseq")


# COMPARING between VERSIONS - rnaseq CRS vs rnaseq CRB
# compare both files
comperison_by_avg <- compare_variance_by_avg_plot(v6_avg_rnaseq_all_sig,"sig - AD CRS avg rnaseq", v7_avg_rnaseq_all_sig, "AD CRB avg rnaseq")


# COMPARING between VERSIONS - microarray CRS vs microarray CRB
# compare both files
comperison_by_avg <- compare_variance_by_avg_plot(v6_avg_microarray_all_sig,"sig - AD CRS avg microarray", v7_avg_microarray_all_sig, "AD CRB avg microarray")

