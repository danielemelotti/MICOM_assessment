# AD-HOC MICOM IMPLEMENTATION

# ***Add description***

### Loading necessary packages
require(seminr)


### Sourcing from relevant scripts
source(file = "R/library-micom.R") #necessary for computation and permutations


### Evaluating Measurement Invariance using MICOM

## Step 1 --- Configural Invariance ---

# Confirmed because we are using exactly the same measurement model for both groups

## Step 2 & 3 --- Compositional and Scalar Invariance ---

evaluate_simpleMICOM <- function(original_model, segments_or_condition, nperms, std = TRUE){ #segments_or_condition can also take a condition as a string e.g. "PERF2 < 6"
  
  ## Extracting the data from the original_model
  data <- original_model$data
  
  ## Splitting the dataset and estimating PLS models for each viable segment
  labels <- data_manipulation(data, original_model, segments_or_condition)
  
  # Split the dataset using the labels
  splitted_groups <- split(data, labels$segment_labels)
  
  # Split the scores
  original_observed_scores <- as.data.frame(original_model$construct_scores) ### IS IT?
  
  labeled_dataset <- cbind(labels$segment_labels, original_observed_scores)
  colnames(labeled_dataset)[1] <- "indicator" #add an indicator column to allow for splitting
  
  splitted_observed_scores <- split(original_observed_scores, labeled_dataset$indicator)
  
  # Re-estimating original_model on the separate groups of data
  estimated_models <- lapply(splitted_groups, function(data) {
    rerun(original_model, data = data)
  })
  
  # Convert pooled dataset into a matrix
  #original_data_matrix <- as.matrix(original_data)
  
  # Standardize the data if required (but leave binary variables unstandardized)
  #data_or_std_data <- if (std) scale(data) else data
  
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
  
  # Reporting observed values
  observed_values_report <- cbind(observed_c, observed_mean_diffs, observed_log_var_ratios)
  
  # Running permutations with permute_all function
  permuted_items <- replicate(nperms,
                              permute_all(data = data, ##########
                                          original_model = original_model,
                                          labels = labels),
                              simplify = FALSE)
  
  # Reporting observed values and permuted values
  combined_results <- list(observed_values_report = observed_values_report,
                           permuted_values_report = permuted_items)
  
  return(combined_results)
  
}
