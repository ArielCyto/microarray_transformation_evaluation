### Functions ###

# function for calc data table and to create the column of "significant" (<=0.1)
generate_variance_table <- function(cttest_data_table, platform_type) {
  data_table <- cttest_data_table[cttest_data_table$Platform == platform_type,]
  cell_list <- unique(data_table$Cell.Typeold)
  final_cell_table <- data.frame("Cell_name" = cell_list, "sig.up/ds" = NA, "sig.down/ds" = NA,
                                 "sig.up/sig" = NA, "sig.down/sig" = NA)
  for (cell in cell_list){
    temp_cell_table <- data_table[data_table$Cell.Typeold == cell,]
    number_of_ds <- dim(temp_cell_table)[1]
    sigUP <- dim(temp_cell_table[temp_cell_table$Estimate.Up.Down == 'Up' & temp_cell_table$significant == 'TRUE',])[1]
    sigDown <- dim(temp_cell_table[temp_cell_table$Estimate.Up.Down == 'Down' & temp_cell_table$significant == 'TRUE',])[1]
    sigTotal <- sigUP+sigDown
    
    temp_cell_table["sig.up/ds"] <- sigUP/number_of_ds
    temp_cell_table["sig.down/ds"] <- sigDown/number_of_ds
    temp_cell_table["sig.up/sig"] <- sigUP/sigTotal
    temp_cell_table["sig.down/sig"] <- sigDown/sigTotal
    
    final_cell_table$sig.up.ds[final_cell_table$Cell_name == cell] <- sigUP/number_of_ds
    final_cell_table$sig.down.ds[final_cell_table$Cell_name == cell] <- sigDown/number_of_ds
    final_cell_table$sig.up.sig[final_cell_table$Cell_name == cell] <- sigUP/sigTotal
    final_cell_table$sig.down.sig[final_cell_table$Cell_name == cell] <- sigDown/sigTotal
  }
  return(final_cell_table)
}

# A function that calculates the direction and score for a given cell 
calc_score_per_cell <- function(table_a, cell){
  if (cell %in% table_a$Cell_name){
    tmp_a <- table_a[table_a$Cell_name == cell,]
    if ( all(!is.na(tmp_a)) ){
      if (tmp_a$sig.down.sig > tmp_a$sig.up.sig){
        direction_a <- "Down"
        direction_a_score <- -1*(table_a$sig.down.sig[table_a$Cell_name == cell])
      } else{
        direction_a <- "Up"
        direction_a_score <- table_a$sig.up.sig[table_a$Cell_name == cell]
      }
    } else{
      direction_a <- NA
      direction_a_score <- NA
    } 
  } else{
    direction_a <- NA
    direction_a_score <- NA
  }
  return(c(direction_a, direction_a_score))
}

# function to create comparison_table (with the direction (up/down) and the score including the direction sign) 
calc_score_cells_two_tables <- function(table_a, table_b){
  cells <- intersect(table_a$Cell_name, table_b$Cell_name)
  comparison_table <- data.frame("Cell_name" = cells, "direction_a" = NA, "direction_a_score" = NA,
                                 "direction_b" = NA, "direction_b_score" = NA)
  for (cell in cells){
    out_score_a <- calc_score_per_cell(table_a, cell)
    out_score_b <- calc_score_per_cell(table_b, cell)
    comparison_table$direction_a[comparison_table$Cell_name == cell] <- out_score_a[1]
    comparison_table$direction_a_score[comparison_table$Cell_name == cell] <- as.numeric(out_score_a[2])
    comparison_table$direction_b[comparison_table$Cell_name == cell] <- out_score_b[1]
    comparison_table$direction_b_score[comparison_table$Cell_name == cell] <- as.numeric(out_score_b[2])
  }
  return(comparison_table)
}


# A function to split long string in the middle
# split_string <- function(long_string) {
#   middle_index <- ceiling(nchar(long_string) / 2) # Find the middle index of the string
#   space_index <- regexpr(" ", substr(long_string, middle_index, nchar(long_string))) # Find the index of the first space after the middle index
#   # If a space is found after the middle index, split the string
#   if (space_index > 0) {
#     split_index <- middle_index + space_index - 1
#     split_string <- paste0(substr(long_string, 1, split_index), "\n", substr(long_string, split_index + 1, nchar(long_string)))
#   } else {
#     # If no space is found after the middle index, return the original string
#     split_string <- long_string
#   }
#   return(split_string)
# }

# A function to split long string into x lines
split_string <- function(long_string, x=2) {
  split_string <- long_string
  middle_index <- ceiling(nchar(long_string) / x) # Find the middle index of the string
  for (i in array(1:x)[1:x-1]){
    space_index <- regexpr(" ", substr(long_string, i*middle_index, nchar(long_string))) # Find the index of the first space after the middle index
    if (space_index > 0) {
      split_index <- i*middle_index + space_index - 1
      split_string <- paste0(substr(split_string, 1, split_index-1), "\n", substr(split_string, split_index + 1, nchar(long_string)))
    } else {
      # If no space is found after the middle index, return the original string
      split_string <- long_string
    }
  }
  return(split_string)
}


# A function that splits a long list into x sublists
split_list <- function(lst, x) {
  sublist_length <- ceiling(length(lst) / x) # Calculate the length of each sublist
  sublists <- split(lst, rep(1:x, each = sublist_length, length.out = length(lst))) # Split the list into sublists
  return(sublists)
}


# A function that compare the variance between 2 model tables:
# arguments: 2 data tables, name per table
compare_variance <- function(table_a, name_a, table_b, name_b){
  comparison_table_microarray <- calc_score_cells_two_tables(table_a, table_b)
  
  # prepare data to plot the scores per cell 
  a <- melt(
    comparison_table_microarray %>% select(Cell_name, direction_a_score, direction_b_score), 
    id.vars =  "Cell_name",
    variable.name = "model.version")
  a$model.version <- str_replace(a$model.version, "direction_a_score", name_a)
  a$model.version <- str_replace(a$model.version, "direction_b_score", name_b)
  
  # split long cell names to fit into the facet title
  for (x in a$Cell_name){
    a$facet[a$Cell_name==x] <- split_string(x)
  }
  
  title_of_plot = paste0("Variance test: ", name_a," vs ", name_b)
  filename = paste0("~/capsule/code/microarray_transformation_evaluation/",title_of_plot,".png")
  
  # plot scores = sig.up/total_sig and sig.down/total_sig (up-regulated > 0, down-regulated <0)
  curr_plot <- a %>% ggplot() +
    geom_point(aes(x=model.version, y=value, color = model.version), shape="\u2605", size = 6) +
    facet_grid(~facet) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size=10),
            strip.text = element_text(angle = 90)) +
    geom_hline(yintercept = 0) +
    ylab("sig(up or down)/sig.total") +
    labs(title = title_of_plot)
  
  print(curr_plot)
  
  ggsave(curr_plot, 
         filename = filename,
         device = "png",
         height = 6, width = 22, units = "in")
  
  print(paste0("The final plot is also saved in ", filename))
  
  
}





