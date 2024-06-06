# AD FDR dot plot per platform - daniel request

V7_data <- read.csv("~/capsule/code/microarray_transformation_evaluation/source_csvs/CT_Test_AD_v7.csv") # CRB
# change the direction of FDR for (-) for Down
V7_data_FDR <- V7_data
V7_data_FDR$fdr[V7_data_FDR$Estimate.Up.Down == 'Down'] <- -1*(V7_data_FDR$fdr[V7_data_FDR$Estimate.Up.Down == 'Down'])

ds_list <- unique(as.character(V7_data_FDR$Experiment.ID)) # 11 ds
comparison_table <- data.frame('cell_names' = unique(V7_data_FDR$Cell.Typeold))

for (ds in ds_list){
    comparison_table[ds] <- NA
}

for (ds in ds_list){
  for (cell in comparison_table$cell_names){
    FDR_vals <- V7_data_FDR$fdr[V7_data_FDR$Cell.Typeold == cell]
    comparison_table[comparison_table$cell_names == cell,2:length(ds_list)] <- FDR_vals
  }
}

# plot it
dm  <- melt(comparison_table,id.var=1)

# dot plot per dataset
g <- ggplot(data=dm,aes(cell_names,value,colour=variable))+ geom_hline(yintercept = 0) +
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggtitle(paste0("Cell Abundance for ", comparison, " comparison per dataset"))
print(g)
