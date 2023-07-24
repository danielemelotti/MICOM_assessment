permute_all <- function(data, original_model, labels) {
  
  data_labeled <- cbind(labels$segment_labels, data)
  colnames(data_labeled)[1] <- "indicator"
  grouped <- permute_and_group(data_labeled, labels$segment_labels) #call the permutation and grouping function
  
  pls_models <- suppressMessages(lapply(grouped, function(data) {
    rerun(original_model, data = data)
  })) # rerunning original_model on the new permuted data groups
  
  ##########################    ##########################    ##########################
  # Split the scores
  original_perm_scores <- lapply(seq_along(pls_models), function(i) {
    df <- as.data.frame(pls_models[[i]]$construct_scores)
    df <- cbind(pls_models[[i]]$rawdata$indicator, df)
    colnames(df)[1] <- "indicator"
    return(df)
  })
  
  splitted_perm_scores <- lapply(original_perm_scores, function(df) {
    split(df, df$indicator)
  })
  
  # Remove 'indicator' variable from each data frame in each sublist
  splitted_perm_scores_no_ind <- lapply(splitted_perm_scores, function(sublist) {
    lapply(sublist, function(df) {
      df[, !names(df) %in% "indicator"]
    })
  })
  ##########################    ##########################    ##########################
  
  outer_weights_list <- extract_outer_weights(pls_models, labels$viable_segments_ID) #get outer weights
  
  composites_list <- compute_composite_scores(data, outer_weights_list,  labels$viable_segments_ID) #compute composite scores
  
  cor_list <- compute_correlations(composites_list) #compute correlations
  
  ##########################    ##########################    ##########################
  ## Compute the observed mean differences
  # Compute observed_means
  perm_means <- lapply(splitted_perm_scores_no_ind, function(x) compute_means(x, construct_names = labels$construct_names))
  
  # Compute observed_mean_diffs
  perm_mean_diffs <- lapply(perm_means, compute_mean_diffs)
  
  ## Compute the log ratios of observed variances
  # Compute observed_vars
  perm_vars <- lapply(splitted_perm_scores_no_ind, function(x) compute_variances(x,construct_names = labels$construct_names)) 
  
  # Compute observed_log_var_ratios
  perm_log_var_ratios <- lapply(perm_vars, compute_log_var_ratios)
  ##########################    ##########################    ##########################
  
  #mean_composites_list <- compute_means(composites_list, labels$construct_names) #compute means
  #mean_diffs <- compute_mean_diffs(mean_composites_list) #compute mean differences
  
  #var_composites_list <- compute_variances(composites_list, labels$construct_names) #compute variances
  #log_var_ratios <- compute_log_var_ratios(var_composites_list) #compute log-ratio of vars
  
  #permutation_report <- cbind(cor_list, mean_diffs, log_var_ratios) #bind results into a dataframe
  ##########################    ##########################    ##########################
  permutation_report <- cbind(cor_list, perm_mean_diffs, perm_log_var_ratios) #bind results into a dataframe
  ##########################    ##########################    ##########################
  
}
