### Set of functions that allow permutations for MICOM procedure

#Note:Some functions can work on more than 2 groups, while some are limited to working on 2 ONLY, which makes it possible to evaluate Measurement Invariance on ONLY 2 groups as for now.

permute_and_group <- function(data, groups) {
  permuted <- data[sample(nrow(data), replace = FALSE), ]
  grouped <- split(permuted, groups) #as.character(groups)
  return(grouped)
}

extract_outer_weights <- function(pls_models, viable_segments_ID) {
  outer_weights_list <- lapply(1:length(pls_models), function(i) {
    model <- pls_models[[i]]
    model$outer_weights
  })
  names(outer_weights_list) <- paste0("segment_", viable_segments_ID, "_weights")
  return(outer_weights_list)
}

compute_composite_scores <- function(data, outer_weights_list, viable_segments_ID) {
  # Convert full dataset into matrix
  data_matrix <- as.matrix(data) #[, -1]
  composites_list <- lapply(1:length(outer_weights_list), function(i) {
    composite_scores <- data_matrix %*% outer_weights_list[[i]]
  })
  names(composites_list) <- paste0("segment_", viable_segments_ID, "_scores")
  return(composites_list)
}


compute_correlations <- function(composites_list) {
  # Generate all combinations of the groups taken two at a time
  combinations <- combn(names(composites_list), 2)
  cor_list <- apply(combinations, 2, function(c) {
    group1 <- composites_list[[c[1]]]
    group2 <- composites_list[[c[2]]]
    cor_matrix <- cor(group1, group2)
    diag(cor_matrix)
  })
  cor_list <- as.data.frame(cor_list)
  colnames(cor_list) <- apply(combinations, 2, paste, collapse = "_")
  return(cor_list)
}


compute_means <- function(composites_list, construct_names) {
  # Compute the means
  mean_composites_list <- lapply(composites_list, function(group) {
    sapply(1:ncol(group), function(i) {
      mean(group[, i])
    })
  })
  
  # Name the elements of mean_composites_list
  mean_composites_list <- lapply(seq_along(mean_composites_list), function(i) {
    names(mean_composites_list[[i]]) <- paste0("mean_", construct_names)
    mean_composites_list[[i]]  # return the modified group
  })
  
  names(mean_composites_list) <- names(composites_list)  # restore the original group names
  
  return(mean_composites_list)
}


compute_mean_diffs <- function(mean_composites_list) {
  # Compute the mean differences
  mean_diffs <- mean_composites_list[[1]] - mean_composites_list[[2]]
  return(mean_diffs)
}


compute_variances <- function(composites_list, construct_names) {
  # Compute the variances
  var_composites_list <- lapply(composites_list, function(group) {
    sapply(1:ncol(group), function(i) {
      var(group[, i])
    })
  })
  
  # Name the elements of var_composites_list
  var_composites_list <- lapply(seq_along(var_composites_list), function(i) {
    names(var_composites_list[[i]]) <- paste0("var_", construct_names)
    var_composites_list[[i]]
  })
  
  names(var_composites_list) <- names(composites_list)  # restore the original group names
  
  return(var_composites_list)
}


compute_log_var_ratios <- function(var_composites_list) {
  # Compute log ratios of the variances
  log_var_ratios <- log(var_composites_list[[1]] / var_composites_list[[2]])
  return(log_var_ratios)
}


### Permutations

permute_all <- function(data, original_model, labels) {
  
  grouped <- permute_and_group(data, labels$segment_labels) #call the permutation and grouping function
  
  pls_models <- suppressMessages(lapply(grouped, function(data) {
    rerun(original_model, data = data)
  })) # rerunning original_model on the new permuted data groups
  
  outer_weights_list <- extract_outer_weights(pls_models, labels$viable_segments_ID) #get outer weights
  
  composites_list <- compute_composite_scores(data, outer_weights_list,  labels$viable_segments_ID) #compute composite scores
  
  cor_list <- compute_correlations(composites_list) #compute correlations
  
  mean_composites_list <- compute_means(composites_list, labels$construct_names) #compute means
  mean_diffs <- compute_mean_diffs(mean_composites_list) #compute mean differences
  
  var_composites_list <- compute_variances(composites_list, labels$construct_names) #compute variances
  log_var_ratios <- compute_log_var_ratios(var_composites_list) #compute log-ratio of vars
  
  permutation_report <- cbind(cor_list, mean_diffs, log_var_ratios) #bind results into a dataframe
  
}


### Data manipulation - Function to extract measurement items' names, constructs' names, and segment names

# data_manipulation <- function(original_data, original_model, segments) {
#   # Creating variables that contain the names of measurement items and constructs
#   item_names <- seminr:::all_items(original_model$measurement_model)
#   construct_names <- seminr:::construct_names(original_model$structural_model)
# 
#   # Extracting viable segment IDs
#   viable_segments_ID <- as.character(segments$viable)
# 
#   # Creating a list containing the node cases of the viable segments
#   viable_segments <- sapply(viable_segments_ID,function(id) node_cases = segments$node_cases[[id]])
# 
#   # Dividing the original pooled data into separate datasets for each segment
#   segment_dataset <- function(dataset, segments) {
#     lapply(segments, function(segment) {
#       dataset[segment, ]
#     })
#   }
# 
#   segmented_datasets <- segment_dataset(original_data, viable_segments)
# 
#   # Creating a vector with group names
#   segment_labels <- rep(NA, nrow(original_data)) #create a new vector to store the segment_id
# 
#   # Add the value of segment_id to each observation (row)
#   sapply(seq_along(viable_segments), function(i) {
#     segment_labels[as.numeric(viable_segments[[i]])] <<- paste0("segment_", viable_segments_ID[i])
#   })
#   
#   return(list(item_names = item_names, construct_names = construct_names, viable_segments_ID = viable_segments_ID, segment_labels = segment_labels))
# }


data_manipulation <- function(original_data, original_model, segments_or_condition) {
  # Creating variables that contain the names of measurement items and constructs
  item_names <- seminr:::all_items(original_model$measurement_model)
  construct_names <- seminr:::construct_names(original_model$structural_model)
  
  if (is.list(segments_or_condition)) { # Case where segments_or_condition is a list from model_segments()
    # Extracting viable segment IDs
    viable_segments_ID <- as.character(segments_or_condition$viable) 
    
    # Creating a list containing the node cases of the viable segments
    viable_segments <- sapply(viable_segments_ID,function(id) node_cases = segments_or_condition$node_cases[[id]])
    
  } else { # Case where segments_or_condition is a condition vector
    viable_segments_ID <- c("satisfy", "not_satisfy")
    
    # Creating a list of indices that meet and do not meet the condition
    satisfy <- which(segments_or_condition)
    not_satisfy <- which(!segments_or_condition)
    viable_segments <- list(satisfy, not_satisfy)
  }
  
  # Dividing the original pooled data into separate datasets for each segment
  segment_dataset <- function(dataset, segments) {
    lapply(segments, function(segment) {
      dataset[segment, ]
    })
  }
  
  segmented_datasets <- segment_dataset(original_data, viable_segments)
  
  # Creating a vector with group names
  segment_labels <- rep(NA, nrow(original_data)) #create a new vector to store the segment_id
  
  # Add the value of segment_id to each observation (row)
  sapply(seq_along(viable_segments), function(i) {
    segment_labels[as.numeric(viable_segments[[i]])] <<- paste0("segment_", viable_segments_ID[i])
  })
  
  return(list(item_names = item_names, construct_names = construct_names, viable_segments_ID = viable_segments_ID, segment_labels = segment_labels))
}

