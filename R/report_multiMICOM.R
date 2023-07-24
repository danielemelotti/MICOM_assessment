report_multiMICOM <- function(multiMICOM_object, significance_level = 0.05) {
  
  model_results <- multiMICOM_object$Results
  
  names(model_results) <- c("Step 2: Compositional Invariance", 
                            "Step 3.1: Equality of Scores' Means", 
                            "Step 3.2: Equality of Scores' Variances")
  
  report <- lapply(names(model_results), function(name) {
    test_result <- model_results[[name]]
    test_result$Decision <- ifelse(is.nan(test_result$p.value), 
                                   "Not estimated",
                                   ifelse(test_result$p.value < significance_level, "No", "Yes"))

    names(test_result) <- c("Statistic", "P-Value", "Invariance between groups?")  # Rename the columns
    return(test_result)
  })
  
  names(report) <- names(model_results)
  
  return(report)
}

