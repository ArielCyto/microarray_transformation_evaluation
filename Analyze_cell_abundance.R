install.packages("curl")
library(curl)
install.packages('cytoreason.ccm.pipeline', repos = c(CRRAN = "https://crran.cytoreason.com/R/latest"))
library("cytoreason.ccm.pipeline")
library(ggplot2)
library(reshape2) ## for melt()
source("~/capsule/code/microarray_transformation_evaluation/estimation_functions.R")

### ~~ CRC ~~ ###
# check low abundance cells in CRC, V34
CRC_model_version <- as_ccm_fit("wf-eea4eae0a5") # CRC_v34
# set the effect_id and comparison to explore (for IO - "effect_id": "tumor","comparison": "tumor_vs_tumor_adjacent")
CRC_effect_id <- "tumor"
CRC_comparison <- "tumor_vs_tumor_adjacent"
CRC_model_groups <- "tumor_vs_tumor_adjacent"  # as taken from model_group column

# should return a boxplot and a df with the results per cell, per dataset
CRC_abundance_table <- generate_estimate_1_boxplot(CRC_model_version,CRC_effect_id,CRC_comparison,CRC_model_groups)        


### ~~ AD ~~ ###
# check low abundance cells in AD (daniel)
AD_model_version <- as_ccm_fit("wf-d605f95788")
# set the effect_id and comparison to explore (for IO - "effect_id": "tumor","comparison": "tumor_vs_tumor_adjacent")
AD_effect_id <- "AD"
AD_comparison <- "AD_vs_HC"
AD_model_groups <- "DZ_vs_HC"  # as taken from model_group column

# should return a boxplot and a df with the results per cell, per dataset
AD_abundance_table <- generate_estimate_1_boxplot(AD_model_version,AD_effect_id,AD_comparison,AD_model_groups)        




# TODO:
# find the median / other cutoff to use
# mark in the analyses of https://cytoreason.sharepoint.com/:p:/r/_layouts/15/Doc.aspx?sourcedoc=%7B0A4B6921-442F-4388-B325-E7BB2635D367%7D&file=V19vsV34_%20comparison.pptx&action=edit&mobileredirect=true who is under the cutoff 
# and who is above, check if it is correlate with the results
