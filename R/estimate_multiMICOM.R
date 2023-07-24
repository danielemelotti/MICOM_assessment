### Loading necessary packages
require(car)

### Sourcing from relevant scripts
source(file = "R/library-micom.R") #necessary for computation and permutations


evaluate_multiMICOM <- function(original_model, segments_or_condition, std = TRUE) {
  
  data <- original_model$data
  
  ## Splitting the dataset and estimating PLS models for each viable segment or condition
  labels <- multi_data_manipulation(data, original_model, segments_or_condition)
  
  # Split the dataset using the labels
  splitted_groups <- split(data, labels$segment_labels)
  
  # Split the scores
  original_observed_scores <- as.data.frame(original_model$construct_scores)
  
  labeled_dataset <- cbind(labels$segment_labels, original_observed_scores)
  colnames(labeled_dataset)[1] <- "indicator" #add an indicator column to allow for splitting
  
  splitted_observed_scores <- split(original_observed_scores, labeled_dataset$indicator)
  
  # Re-estimating original_model on the separate groups of data
  estimated_models <- lapply(splitted_groups, function(data) {
    rerun(original_model, data = data)
  })
  
  # Standardize the data if required (but leave binary variables unstandardized)
  data_or_std_data <- if(std) {
    data.frame(apply(data, 2, function(column) {
      # Check if the column is binary
      if (length(unique(column)) == 2) {
        return(column)  # Return original values if binary
      } else {
        return(scale(column))  # Return standardized values otherwise
      }
    }))
  } else {
    data  # Return original data if standardization is not requested
  }
  
  ## Retrieve the weights for each segment into a list
  observed_outer_weights <- extract_outer_weights(estimated_models, labels$viable_segments_ID)
  
  ## Compute the observed correlation c between the composite scores
  # Compute the composite scores
  observed_composites <- compute_composite_scores(data = data_or_std_data,
                                                  outer_weights_list = observed_outer_weights, 
                                                  viable_segments_ID = labels$viable_segments_ID)
  
  # Compute observed_c
  observed_c <- compute_correlations(observed_composites) 
  
  ## Compute the observed mean differences
  # Compute observed_means
  observed_means <- compute_means(splitted_observed_scores, construct_names = labels$construct_names)
  
  # Compute pairwise differences of means
  observed_mean_diffs <- compute_pairwise_mean_diffs_df(observed_means)
  
  ## Compute the log ratios of observed variances
  # Compute observed_vars
  observed_vars <- compute_variances(splitted_observed_scores, labels$construct_names) 
  
  # Compute pairwise ratios of variances
  observed_log_var_ratios <- compute_pairwise_log_var_ratios_df(observed_vars)
  
  ## Step 2: Friedman test
  friedman_tests <- compute_friedman(observed_composites, labels$construct_names)
  
  ## Step 3.1: Kruskal-Wallis test
  kruskal_tests <- compute_kruskal(splitted_observed_scores, labels$construct_names)
  
  ## Step 3.2: Levene test
  levene_tests <- compute_levenes(splitted_observed_scores, labels$construct_names)
  
  # Get only the statistics and p-values
  friedman_tests_stats <- do.call(rbind, lapply(friedman_tests, function(x) data.frame('statistic' = x$statistic, 'p.value' = x$p.value)))
  kruskal_tests_stats <- do.call(rbind, lapply(kruskal_tests, function(x) data.frame('statistic' = x$statistic, 'p.value' = x$p.value)))
  levene_tests_stats <- do.call(rbind, lapply(levene_tests, function(x) data.frame('statistic' = x$statistic, 'p.value' = x$p.value)))
  
  tests_results <- list('Friedman' = friedman_tests_stats, 'Kruskal-Wallis' = kruskal_tests_stats, 'Levene' = levene_tests_stats)
  
  # Get the full results
  full_results <- list('Friedman' = friedman_tests, 'Kruskal-Wallis' = kruskal_tests, 'Levene' = levene_tests)
  
  results <- list("Results" = tests_results,
                  "Full Results" = full_results,
                  "Observed Composites" = observed_composites,
                  "Splitted Observed Scores" = splitted_observed_scores)
  
  return(results)
}

