### Loading necessary packages
require(seminr)
require(car)

### Sourcing from relevant scripts
source(file = "R/library-micom.R") #necessary for computation and permutations


### Evaluating Measurement Invariance using MICOM

## Step 1 --- Configural Invariance ---

# Confirmed because we are using exactly the same measurement model for both groups

## Step 2 & 3 --- Compositional and Scalar Invariance ---

evaluate_robustMICOM <- function(original_model, segments_or_condition, std = TRUE){ #segments_or_condition can also take a condition as a string e.g. "PERF2 < 6"
  
  ## Extracting the data from the original_model
  data <- original_model$data
  
  ## Splitting the dataset and estimating PLS models for each viable segment
  labels <- data_manipulation(data, original_model, segments_or_condition)
  
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
  
  # Compute observed_mean_diffs
  observed_mean_diffs <- compute_mean_diffs(observed_means)
  
  ## Compute the log ratios of observed variances
  # Compute observed_vars
  observed_vars <- compute_variances(splitted_observed_scores, labels$construct_names) 
  
  # Compute observed_log_var_ratios
  observed_log_var_ratios <- compute_log_var_ratios(observed_vars)
  
  ## Step 2
  wilcox_observed_composites <- sapply(1:ncol(observed_composites[[1]]), function(i) {
    perform_wilcox_step2(observed_composites[[1]][,i], observed_composites[[2]][,i])
  })
  
  ## Step 3.1
  wilcox_means <- sapply(1:ncol(splitted_observed_scores[[1]]), function(i) {
    perform_wilcox_step3(splitted_observed_scores[[1]][,i], splitted_observed_scores[[2]][,i])
  })
  
  ## Step 3.2
  # Compute Levene's test for each construct
  levene <- mapply(perform_levene, splitted_observed_scores[[1]], splitted_observed_scores[[2]])
  
  # Prepare results
  step2 <- data.frame(observed_c = observed_c,
                    Wilcoxon_composites_p_value = wilcox_observed_composites)

  step3_1 <- data.frame(observed_mean_diffs = observed_mean_diffs,
                    Wilcoxon_means_p_value = wilcox_means)

  step3_2 <- data.frame(observed_log_var_ratios = observed_log_var_ratios,
                    Levene_vars_p_value = levene)

  # Combine them into a list
  observed_values_report <- list(step2, step3_1, step3_2)

  names(observed_values_report) <- c("Step 2: Compositional Invariance",
                                     "Step 3.1: Equality of Means",
                                     "Step 3.2: Equality of Variances")
  
  # Return the list of data frames
  return(observed_values_report)
  
}
