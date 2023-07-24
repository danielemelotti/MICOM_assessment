### Reporting

source(file = "R/library-micom.R")


report_simpleMICOM <- function(micom_results, confidence_level = 0.95) {
  # Calculate the quantiles for the confidence interval
  lower_quantile <- (1 - confidence_level) / 2
  upper_quantile <- 1 - lower_quantile
  
  # Obtaining construct names
  construct_names <- rownames(micom_results$observed_values_report)
  
  # Step 2: Compositional Invariance
  compositional_invariance <- calculate_compositional_values(micom_results$observed_values_report[, 1],
                                                             micom_results$permuted_values_report,
                                                             confidence_level = confidence_level)
  compositional_invariance <- cbind(construct_names, data.frame(t(compositional_invariance)))
  colnames(compositional_invariance) <- c("Composite", "Observed c Value", "Mean Permuted c", "5%", "p-value", "Compositional invariance?")
  
  # Step 3.1: Equality of Mean Values
  equal_mean_values <- calculate_values(observed_values = micom_results$observed_values_report[, 2],
                                        permuted_values = micom_results$permuted_values_report,
                                        perm_column = 2,
                                        confidence_level = confidence_level,
                                        lower_quantile = lower_quantile,
                                        upper_quantile = upper_quantile)
  equal_mean_values <- cbind(construct_names,data.frame(t(equal_mean_values)))
  colnames(equal_mean_values) <- c("Composite", "Observed Mean Difference", "Mean Permuted Difference", "2.5%", "97.5%","p-value", "Equal mean values?")
  
  # Step 3.2: Equality of Variances
  equal_variances <- calculate_values(observed_values = micom_results$observed_values_report[, 3],
                                      permuted_values = micom_results$permuted_values_report, 
                                      perm_column = 3,
                                      lower_quantile = lower_quantile,
                                      upper_quantile = upper_quantile)
  equal_variances <- cbind(construct_names, data.frame(t(equal_variances)))
  colnames(equal_variances) <- c("Composite", "Observed Log of Variance Ratio", "Mean Permuted Logs", "2.5%", "97.5%", "p-value", "Equal variances?")
  
  # Combine into a list
  report <- list(
    "Step 2: Compositional Invariance" = compositional_invariance,
    "Step 3.1: Equality of Mean Values" = equal_mean_values,
    "Step 3.2: Equality of Variances" = equal_variances
  )
  
  return(report)
}
