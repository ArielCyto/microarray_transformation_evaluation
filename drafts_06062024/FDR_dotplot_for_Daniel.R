# AD estimate dot plot per platform - daniel request

V7_data <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/AD_CRB_Estimation_Data.csv") # CRB

V7_data_estimate <- V7_data
ds_list <- unique(as.character(V7_data_estimate$Experiment.ID)) # 11 ds
comparison_table <- data.frame('cell_names' = unique(V7_data_estimate$Cell.Type))

for (ds in ds_list){
  comparison_table[ds] <- NA
}

for (ds in ds_list){
  for (cell in comparison_table$cell_names){
    estimate_vals <- V7_data_estimate$Estimate[V7_data_estimate$Cell.Type == cell]
    comparison_table[comparison_table$cell_names == cell,2:length(ds_list)] <- estimate_vals
  }
}

# plot it
comparison_table$meta <- NULL
dm  <- melt(comparison_table,id.var=1)

# dot plot per dataset
g <- ggplot(data=dm,aes(cell_names,value,colour=variable))+ geom_hline(yintercept = 0) +
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggtitle(paste0("Cell Abundance for ", comparison, " comparison per dataset"))
print(g)
