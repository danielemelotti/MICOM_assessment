### Loading necessary packages
require(PMCMRplus)
require(conover.test)
require(FSA)

source(file = "R/library-micom.R")

multiMICOM_posthoc <- function(data, constructs = NULL, alpha = 0.05) {
  # Step 1: Check inputs
  #if (!is.null(constructs)) {
  #  if (!all(constructs %in% colnames(data))) {
  #    stop("All construct names must match column names in the data frame.")
  #  }
  #} else {
  #  constructs <- colnames(data)[!(colnames(data) %in% grouping_variable)]
  #}
  
  # Step 2: Post-hoc tests
  dunn_tests <- compute_dunn_posthoc(data, constructs)
  conover_tests <- compute_conover_posthoc(data, constructs)
  gameshowell_tests <- compute_gameshowell_posthoc(data, constructs)
  
  posthoc_results <- list("Dunn" = dunn_tests, "Conover" = conover_tests, "Games-Howell" = gameshowell_tests)
  
  return(posthoc_results)
}
