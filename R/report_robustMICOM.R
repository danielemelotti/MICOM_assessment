report_robustMICOM <- function(micom_results, confidence = 0.95) {
  
  significance <- 1 - confidence
  
  report_results <- lapply(micom_results, function(results) {
    results$Invariance <- ifelse(results[,2] < significance, "No", "Yes")
    rownames(results) <- sub("mean_|var_", "", rownames(results))
    return(results)
  })
  
  # Rename column names
  names(report_results[[1]])[1] <- "c value"
  names(report_results[[2]])[1] <- "Difference of the composite’s mean value"
  names(report_results[[3]])[1] <- "Logarithm of the composite’s variances ratio"
  
  names(report_results) <- paste(names(report_results), "?", sep = " ")
  
  return(report_results)
}
