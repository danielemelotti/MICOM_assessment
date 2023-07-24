### Set of functions that allow permutations for MICOM procedure


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


### Dealing with missing values - WORK IN PROGRESS


### Permutations

permute_all <- function(data, original_model, labels) {

  grouped <- permute_and_group(data, labels$segment_labels) # for step 2

  permuted_pooled_data <- data[sample(nrow(data), replace = FALSE), ] # for step 3
  
  pooled_pls <- suppressMessages(estimate_pls(data = permuted_pooled_data,
                                              measurement_model = original_model$measurement_model,
                                              structural_model = original_model$structural_model))
  
  pls_models <- suppressMessages(lapply(grouped, function(data) {
    rerun(original_model, data = data)
  })) # rerunning original_model on the new permuted data groups
  
  outer_weights_list <- extract_outer_weights(pls_models, labels$viable_segments_ID) #get outer weights

  composites_list <- compute_composite_scores(data, outer_weights_list,  labels$viable_segments_ID) #compute composite scores

  permuted_cor_list <- compute_correlations(composites_list) #compute correlations

  ## For equality of means and variances
  # Split the scores
  permuted_scores <- as.data.frame(pooled_pls$construct_scores) 
  
  labeled_permuted_dataset <- cbind(labels$segment_labels, permuted_scores)
  colnames(labeled_permuted_dataset)[1] <- "indicator" #add an indicator column to allow for splitting
  
  splitted_permuted_scores <- split(permuted_scores, labeled_permuted_dataset$indicator)
  
  #splitted_permuted_scores <- permute_and_group(pooled_pls$construct_scores, labels$segment_labels) # call the permutation and grouping function for original scores

  # Compute permuted_means
  permuted_means <- compute_means(splitted_permuted_scores, construct_names = labels$construct_names)

  # Compute permuted_mean_diffs
  permuted_mean_diffs <- compute_mean_diffs(permuted_means)

  ## Compute the log ratios of permuted variances
  # Compute permuted_vars
  permuted_vars <- compute_variances(splitted_permuted_scores, labels$construct_names)

  # Compute permuted_log_var_ratios
  permuted_log_var_ratios <- compute_log_var_ratios(permuted_vars)

  permutation_report <- cbind(permuted_cor_list, permuted_mean_diffs, permuted_log_var_ratios) #bind results into a dataframe

}


### Data manipulation - Function to extract measurement items' names, constructs' names, and segment names

data_manipulation <- function(original_data, original_model, segments_or_condition) {
  # Creating variables that contain the names of measurement items and constructs
  item_names <- seminr:::all_items(original_model$measurement_model)
  construct_names <- seminr:::construct_names(original_model$structural_model)
  
  if (is.list(segments_or_condition)) { # Case where segments_or_condition is a list from model_segments()
    # Extracting viable segment IDs
    viable_segments_ID <- as.character(segments_or_condition$viable) 
    
    # Creating a list containing the node cases of the viable segments
    viable_segments <- sapply(viable_segments_ID,function(id) node_cases = segments_or_condition$node_cases[[id]])
    
  } else { # Case where segments_or_condition is a condition string
    condition_string <- segments_or_condition
    viable_segments_ID <- c(condition_string, paste0("not (", condition_string, ")"))
    
    # Creating a list of indices that meet and do not meet the condition
    satisfy <- which(eval(parse(text=condition_string), envir = original_data))
    not_satisfy <- which(!eval(parse(text=condition_string), envir = original_data))
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

### Reporting functions
# Special function for compositional invariance
calculate_compositional_values <- function(observed_values, permuted_values, confidence_level) {
  # Loop over the 11 constructs
  sapply(1:length(observed_values), function(idx) {
    
    observed <- observed_values[idx]
    perms <- unlist(lapply(permuted_values, function(x) x[idx, 1])) # 1st column is for correlations
    
    mean_perm <- mean(perms)
    p5 <- quantile(perms, 1 - confidence_level)
    invariant <- ifelse(observed >= p5, "Yes", "No")
    
    # p-value calculation: proportion of permuted values less than the observed
    p_value <- sum(perms < observed) / length(perms)
    
    c(round(c(observed, mean_perm, p5), 4), round(p_value, 4), invariant)
  })
}


# This function will calculate the values for each construct
calculate_values <- function(observed_values, permuted_values, perm_column, confidence_level, lower_quantile, upper_quantile) {
  
  sapply(1:length(observed_values), function(idx) {
    
    observed <- observed_values[idx]
    perms <- unlist(lapply(permuted_values, function(x) x[idx, perm_column]))
    
    mean_perm <- mean(perms)
    ci <- quantile(perms, c(lower_quantile, upper_quantile))
    invariant <- ifelse(observed >= ci[1] & observed <= ci[2], "Yes", "No")
    
    # p-value calculation: proportion of permuted values outside the confidence interval
    #p_value <- sum(perms > ci[2] | perms < ci[1]) / length(perms)
    p_value <- sum(abs(perms) > abs(observed)) / length(perms)
    
    c(round(c(observed, mean_perm, ci), 4), round(p_value, 4), invariant)
  })
}


### For Wilcoxon tests in step 2 - Paired
perform_wilcox_step2 <- function(scores1, scores2) {
  # Adding very small random noise
  scores1 <- jitter(scores1, amount = 1e-10)
  scores2 <- jitter(scores2, amount = 1e-10)

  test_result <- wilcox.test(scores1, scores2, paired = TRUE)
  return(round(test_result$p.value, 3))
}


### For Wilcoxon tests in step 3.1 - Unpaired
perform_wilcox_step3 <- function(scores1, scores2) {
  # Adding very small random noise
  scores1 <- jitter(scores1, amount = 1e-10)
  scores2 <- jitter(scores2, amount = 1e-10)
  
  test_result <- wilcox.test(scores1, scores2, paired = FALSE)
  return(round(test_result$p.value, 3))
}

### For Levene tests in step 3.2
perform_levene <- function(var1, var2) {
  stacked_vars <- stack(list(group1 = var1, group2 = var2))
  test_result <- leveneTest(values ~ ind, data = stacked_vars)

  return(round(test_result$`Pr(>F)`[1], 3))
}

### For evaluate_multiMICOM:

multi_data_manipulation <- function(original_data, original_model, segments_or_conditions) {
  # Creating variables that contain the names of measurement items and constructs
  item_names <- seminr:::all_items(original_model$measurement_model)
  construct_names <- seminr:::construct_names(original_model$structural_model)
  
  if (is.list(segments_or_conditions)) { # Case where segments_or_condition is a list from model_segments()
    # Extracting viable segment IDs
    viable_segments_ID <- as.character(segments_or_conditions$viable) 
    
    # Creating a list containing the node cases of the viable segments
    viable_segments <- sapply(viable_segments_ID,function(id) node_cases = segments_or_conditions$node_cases[[id]])
    
  } else if (is.character(segments_or_conditions)) { # Case where segments_or_condition is an array of condition strings
    # Creating a list of indices that meet each condition
    viable_segments <- lapply(segments_or_conditions, function(condition_string) {
      return(which(eval(parse(text=condition_string), envir = original_data)))
    })
    viable_segments_ID <- segments_or_conditions
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
  
  return(list(item_names = item_names, construct_names = construct_names, viable_segments_ID = viable_segments_ID, segment_labels = segment_labels, segmented_datasets = segmented_datasets))
}

# Function to compute pairwise differences of means and return as dataframe
compute_pairwise_mean_diffs_df <- function(means_list) {
  construct_names <- names(means_list[[1]])
  diff_df <- data.frame(matrix(ncol = choose(length(means_list), 2), nrow = length(construct_names)))
  colnames(diff_df) <- combn(paste0("Condition_", 1:length(means_list)), 2, paste, collapse = " - ")
  rownames(diff_df) <- construct_names
  
  for (construct in construct_names) {
    group_means <- sapply(means_list, function(group) group[[construct]])
    diffs <- combn(group_means, 2, diff)
    diff_df[construct, ] <- diffs
  }
  diff_df
}

# Function to compute pairwise log ratios of variances and return as dataframe
compute_pairwise_var_ratios_df <- function(vars_list) {
  construct_names <- names(vars_list[[1]])
  ratio_df <- data.frame(matrix(ncol = choose(length(vars_list), 2), nrow = length(construct_names)))
  colnames(ratio_df) <- combn(paste0("Condition_", 1:length(vars_list)), 2, paste, collapse = " / ")
  rownames(ratio_df) <- construct_names
  
  ratio_df <- apply(ratio_df, 1, function(row) {
    group_vars <- sapply(vars_list, function(group) group[[rownames(ratio_df)[row]]])
    ratios <- combn(group_vars, 2, function(x) log(x[1] / x[2]))
    return(ratios)
  })
  ratio_df
}

# Function to compute pairwise log ratios of variances and return as dataframe
compute_pairwise_log_var_ratios_df <- function(vars_list) {
  construct_names <- names(vars_list[[1]])
  ratio_df <- data.frame(matrix(ncol = choose(length(vars_list), 2), nrow = length(construct_names)))
  colnames(ratio_df) <- combn(paste0("Condition_", 1:length(vars_list)), 2, paste, collapse = " / ")
  rownames(ratio_df) <- construct_names
  
  for (construct in construct_names) {
    group_vars <- sapply(vars_list, function(group) group[construct])
    ratios <- combn(group_vars, 2, function(x) log(x[1] / x[2]))
    ratio_df[construct, ] <- ratios
  }
  ratio_df
}

# Function to compute Friedman test
compute_friedman <- function(observed_composites, construct_names) {
  # Prepare and run the Friedman test for each construct
  friedman_tests <- lapply(construct_names, function(construct) {
    # Preparing data for the Friedman test: each construct's scores by segment
    friedman_data <- do.call(cbind, lapply(observed_composites, function(scores) {
      scores[, construct]
    }))
    
    # Run the Friedman test
    test_result <- friedman.test(as.matrix(friedman_data))
    
    # Return the full test result, the test statistic and the p-value
    list('full' = test_result, 
         'statistic' = round(test_result$statistic, 4), 
         'p.value' = round(test_result$p.value, 4))
  })
  names(friedman_tests) <- construct_names
  return(friedman_tests)
}

# Function to compute Kruskal-Wallis test
compute_kruskal <- function(splitted_composites, construct_names) {
  # Prepare and run the Kruskal-Wallis test for each construct
  kruskal_tests <- lapply(construct_names, function(construct) {
    # Prepare data for the Kruskal-Wallis test: each construct's scores by segment
    kruskal_data <- unlist(lapply(splitted_composites, function(scores) {
      scores[, construct]
    }))
    
    # Prepare groups factor for the Kruskal-Wallis test
    kruskal_groups <- rep(names(splitted_composites), sapply(splitted_composites, nrow))
    
    # Run the Kruskal-Wallis test
    test_result <- kruskal.test(kruskal_data, g = kruskal_groups)
    
    # Return the full test result, the test statistic and the p-value
    list('full' = test_result, 
         'statistic' = round(test_result$statistic, 4), 
         'p.value' = round(test_result$p.value, 4))
  })
  names(kruskal_tests) <- construct_names
  return(kruskal_tests)
}

# Function to compute Levene's test
compute_levenes <- function(observed_composites, construct_names) {
  
  # Prepare and run Levene's test for each construct
  levene_tests <- lapply(construct_names, function(construct) {
    
    # Preparing data for the Levene's test: each construct's scores by segment
    levene_data <- do.call(rbind, lapply(seq_along(observed_composites), function(i) {
      data.frame(score = observed_composites[[i]][, construct], group = as.factor(i))
    }))
    
    # Run the Levene's test
    test_result <- car::leveneTest(score ~ group, data = levene_data)
    
    # Return the full test result, the test statistic and the p-value
    list('full' = test_result, 
         'statistic' = round(test_result[1, "F value"], 4), 
         'p.value' = round(test_result[1, "Pr(>F)"], 4))
  })
  names(levene_tests) <- construct_names
  return(levene_tests)
}


### Post-hoc test functions

compute_conover_posthoc <- function(full_results, constructs) {
  # Compute the Conover post-hoc test for each construct
  conover_tests <- lapply(constructs, function(construct) {
    # Prepare data for the Conover test: residuals and grouping factor
    test_data <- full_results$Kruskal_Wallis[[construct]]$data
    group <- full_results$Kruskal_Wallis[[construct]]$group
    
    # Run the Conover test
    test_result <- PMCMRplus::frdAllPairsConoverTest(test_data, group)
    test_result
  })
  names(conover_tests) <- constructs
  conover_tests
}

compute_dunn_posthoc <- function(full_results, constructs) {
  # Compute the Dunn post-hoc test for each construct
  dunn_tests <- lapply(constructs, function(construct) {
    # Prepare data for the Dunn test: residuals and grouping factor
    test_data <- full_results$Kruskal_Wallis[[construct]]$data
    group <- full_results$Kruskal_Wallis[[construct]]$group
    
    # Run the Dunn test
    test_result <- PMCMRplus::kwAllPairsDunnTest(test_data, group)
    test_result
  })
  names(dunn_tests) <- constructs
  dunn_tests
}

compute_gameshowell_posthoc <- function(full_results, constructs) {
  # Compute the Games-Howell post-hoc test for each construct
  gameshowell_tests <- lapply(constructs, function(construct) {
    # Prepare data for the Games-Howell test: residuals and grouping factor
    test_data <- full_results$Levene[[construct]]$residuals
    group <- full_results$Levene[[construct]]$group
    
    # Run the Games-Howell test
    test_result <- PMCMRplus::gamesHowellTest(test_data, group)
    test_result
  })
  names(gameshowell_tests) <- constructs
  gameshowell_tests
}
