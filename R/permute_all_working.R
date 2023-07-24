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
