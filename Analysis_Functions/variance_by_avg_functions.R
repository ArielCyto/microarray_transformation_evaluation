### Functions to calculate variance by average of significant results and their direction (the 'new' format) ###


# function for calc data table and based on the "significant" column (all "TRUE" vals are included)
generate_variance_table_by_avg <- function(cttest_data_table, platform_type) {
  data_table <- cttest_data_table[cttest_data_table$Platform == platform_type,]
  cell_list <- unique(data_table$Cell.Typeold)
  final_cell_table <- data.frame("Cell_name" = cell_list, "directions/ds" = NA)
  for (cell in cell_list){
    temp_cell_table <- data_table[data_table$Cell.Typeold == cell,]
    number_of_ds <- dim(temp_cell_table)[1]
    UpDown_df <- as.data.frame(table(temp_cell_table$Estimate.Up.Down[temp_cell_table$significant == 'TRUE']))
    
    sumUp <- 1*(UpDown_df$Freq[UpDown_df$Var1 == 'Up'])
    if(length(sumUp)==0){sumUp = 0}
    sumDown <- -1*(UpDown_df$Freq[UpDown_df$Var1 == 'Down'])
    if(length(sumDown)==0){sumDown = 0}
    
    temp_cell_table$directions.ds <- (sumUp+sumDown)/number_of_ds
    if(sumDown == 0 & sumUp == 0){temp_cell_table$directions.ds <- NA}
      
    final_cell_table$directions.ds[final_cell_table$Cell_name == cell] <- unique(temp_cell_table["directions.ds"]) 
  }
  return(final_cell_table)
}


# A function that compare the variance between 2 model tables:
# arguments: 2 data tables, name per table
compare_variance_by_avg <- function(table_a, name_a, table_b, name_b){
  cells <- intersect(table_a$Cell_name, table_b$Cell_name)
  comparison_table <- data.frame("Cell_name" = cells,"direction_a" = NA, "direction_a_score" = NA,"direction_b" = NA, "direction_b_score" = NA)
  for (cell in cells){
    comparison_table$direction_a_score[comparison_table$Cell_name == cell] <- as.numeric(table_a$directions.ds[table_a$Cell_name == cell])
    #comparison_table$direction_a[comparison_table$Cell_name == cell] <- sign(as.numeric(comparison_table$direction_a_score[comparison_table$Cell_name == cell]))
    comparison_table$direction_b_score[comparison_table$Cell_name == cell] <- as.numeric(table_b$directions.ds[table_b$Cell_name == cell])
    #comparison_table$direction_b[comparison_table$Cell_name == cell] <- sign(as.numeric(comparison_table$direction_a_score[comparison_table$Cell_name == cell]))
  }
  
  return(comparison_table)
}


# A function that compare the variance between 2 model tables:
# arguments: 2 data tables, name per table
compare_variance_by_avg_plot <- function(table_a, name_a, table_b, name_b){
  comparison_table_avg <- compare_variance_by_avg(table_a,name_a, table_b, name_b)
  # prepare data to plot the scores per cell 
  a <- melt(
    comparison_table_avg %>% select(Cell_name, direction_a_score, direction_b_score), 
    id.vars =  "Cell_name",
    variable.name = "model.version")
  a$model.version <- str_replace(a$model.version, "direction_a_score", name_a)
  a$model.version <- str_replace(a$model.version, "direction_b_score", name_b)
  
  # split long cell names to fit into the facet title
  for (x in a$Cell_name){
    a$facet[a$Cell_name==x] <- split_string(x)
  }
  
  title_of_plot = paste0("Variance test: ", name_a," vs ", name_b)
  filename = paste0("~/capsule/code/microarray_transformation_evaluation/09_06_24/",title_of_plot,".png")
  
  # plot scores = sig.up/total_sig and sig.down/total_sig (up-regulated > 0, down-regulated <0)
  curr_plot <- a %>% ggplot() +
    geom_point(aes(x=model.version, y=value, color = model.version), shape="\u2605", size = 6) +
    facet_grid(~facet) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size=10),
          strip.text = element_text(angle = 90)) +
    geom_hline(yintercept = 0) +
    ylab("sig_sum(Up)-sig_sum(Down)/N_datasets") +
    labs(title = title_of_plot)
  
  print(curr_plot)
  
  ggsave(curr_plot, 
         filename = filename,
         device = "png",
         height = 6, width = 22, units = "in")
  
  print(paste0("The final plot is also saved in ", filename))
  
  return(a)
  
}
