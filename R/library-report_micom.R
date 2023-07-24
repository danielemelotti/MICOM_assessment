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
