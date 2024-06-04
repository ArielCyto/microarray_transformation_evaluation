library(tidyverse)
library(reshape2)
library(dplyr)
library(stringr)
source("~/capsule/code/microarray_transformation_evaluation/variance_functions.R")
#install.packages("curl")
library(curl)
#install.packages('cytoreason.ccm.pipeline', repos = c(CRRAN = "https://crran.cytoreason.com/R/latest"))
library("cytoreason.ccm.pipeline")
library(ggplot2)
source("~/capsule/code/microarray_transformation_evaluation/estimation_functions.R")


# explore the correlation between abundace and variance in V34 CRC - DRAFT
res <- compare_variance(V34_variance_microarray, "microarray V34 CRC", V34_variance_rnaseq, "rnaseq V34 CRC")
CRC_abundance_table <- generate_estimate_1_boxplot(CRC_model_version,CRC_effect_id,CRC_comparison,CRC_model_groups)        

# check the median abundance per cell
median_abundance_CRC_table <- data.matrix(CRC_abundance_table[,-1])
rownames(median_abundance_CRC_table) <- CRC_abundance_table$cell_names
median_abundance_CRC_table <- cbind(median_abundance_CRC_table, "cell_median_abundance" =rowMedians(median_abundance_CRC_table))

#sort by cell_median_abundance
median_abundance_CRC_table <- data.frame(median_abundance_CRC_table)
median_abundance_CRC_table <- median_abundance_CRC_table %>%
  arrange(desc(cell_median_abundance))

# we can see that the top 7 HIGHEST abundanced cells are having the same direction in both platforms - rnaseq or microarray
top7_abundance_cells <- head(rownames(median_abundance_CRC_table),7)
res[res$Cell_name %in% top7_abundance_cells,] %>% arrange(desc(Cell_name))

# we can see that the top 7 LOWEST abundanced cells are having inconsistent direction in both platforms - sometimes the direction is same,
# sometimes we have NA for RNAseq (no significant results), sometimes the directions are opposite
bottom7_abundance_cells <- tail(rownames(median_abundance_CRC_table),7)
res[res$Cell_name %in% bottom7_abundance_cells,] %>% arrange(desc(Cell_name))

# TODO - to write a code that checks the correlation!

