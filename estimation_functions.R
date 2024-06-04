### function for estimating the new transformation results ###

generate_estimate_1_boxplot <- function(model_version, effect_id, comparison, model_group){
  library(ggplot2)
  library(reshape2)
  # generate results_cell for "__tumor_vs_tumor_adjacent" datasets
  ds_list <- names(model_version$datasets) # 12 ds
  comparison_table <- data.frame('cell_names' = rownames(model_version$datasets[[1]]$cell_contribution$eset@featureData@data))
  for (ds in ds_list){
    string_of_data <- paste0(effect_id,"__",ds,"__",comparison)
    if (!is.null(model_version$datasets[[ds]]$model[[string_of_data]]))
      comparison_table[ds] <- NA
  }
  # add the values of T_Tadj comparison
  for (ds in colnames(comparison_table)[2:ncol(comparison_table)]){
    string_of_data <- paste0(effect_id,"__",ds,"__",comparison)
    comparison_table[ds] <- model_version$datasets[[ds]]$model[[string_of_data]]@analysis$ct_test$fit[[model_group]]$pvalues$estimate_1
  }
  
  # plot it
  dm  <- melt(comparison_table,id.var=1)
  
  # dot plot per dataset
  g <- ggplot(data=dm,aes(cell_names,value,colour=variable))+
    geom_point(size=4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    ggtitle(paste0("Cell Abundance for ", comparison, " comparison per dataset"))
  print(g)

  # box plot
  p <- ggplot(dm, aes(x=cell_names, y=value,colour=cell_names)) + geom_boxplot()+
    theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("Cell Abundance for ", comparison, " comparison"))
  print(p)
  return(comparison_table)
  
}
