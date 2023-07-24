### VIABLE SEGMENTS EXTRACTION

# ***Add description***

### Loading necessary packages
require(seminr)
require(semcoa)
library(rpart)


### Sourcing from other scripts

# Segmentation scripts
source(file = "R/semcoa-library/pls_predict.R")
source(file = "R/semcoa-library/segmentation.R")
source(file = "R/semcoa-library/tree_extract.R")
source(file = "R/semcoa-library/tree_grow.R")
source(file = "R/semcoa-library/unstable.R")

# MICOM scripts
source(file = "R/estimate_simpleMICOM.R")
source(file = "R/estimate_robustMICOM.R")
source(file = "R/estimate_multiMICOM.R")
#source(file = "R/report_micom.R")
source(file = "R/report_simpleMICOM.R")
source(file = "R/report_robustMICOM.R")
source(file = "R/report_multiMICOM.R")

### Models declaration

# Measurement model- mode_A
utaut_mm <- constructs(
  composite("PE", multi_items("PERF", 1:4)),
  composite("EE", c("PEOU1","PEOU3","PEOU5","PEOU6")),
  composite("SI", c(multi_items("NORM", 1:2),"INFL3")),
  composite("FC", multi_items("FACL", 1:4)),
  composite("HM", multi_items("MOTIV", 1:3)),
  composite("PV", multi_items("VALUE", 1:3)),
  composite("HAB", multi_items("HAB", 1:4)),
  composite("BI", multi_items("INT", 1:3)),
  composite("Exp", single_item("Experience")),
  composite("Age", single_item("age")),
  composite("Gender", single_item("gender"))
)

# Measurement model - mode_B
utaut_mm_mode_b <- constructs(
  composite("PE", multi_items("PERF", 1:4), weights = mode_B),
  composite("EE", c("PEOU1","PEOU3","PEOU5","PEOU6"), weights = mode_B),
  composite("SI", c(multi_items("NORM", 1:2),"INFL3"), weights = mode_B),
  composite("FC", multi_items("FACL", 1:4), weights = mode_B),
  composite("HM", multi_items("MOTIV", 1:3), weights = mode_B),
  composite("PV", multi_items("VALUE", 1:3), weights = mode_B),
  composite("HAB", multi_items("HAB", 1:4), weights = mode_B),
  composite("BI", multi_items("INT", 1:3), weights = mode_B),
  composite("Exp", single_item("Experience")),
  composite("Age", single_item("age")),
  composite("Gender", single_item("gender"))
)

# Structural Model
utaut_sm <- relationships(
  paths(from = c("PE","EE","SI","FC","HM","PV","HAB","Exp","Age","Gender"), 
        to = "BI")
)

# Creating a variable that contain the names of measurement items
item_names <- seminr:::all_items(utaut_mm)
construct_names <- seminr:::construct_names(utaut_sm)

# Loading the dataset, making sure to import only the items used in the measurement model. This is a technical requirement that allows the MICOM calculations to be executed later in the script (matrix multiplication between pooled data and outer_weights)
utaut_data <- read.csv(file = "data/trello_utaut.csv")[, item_names]

### Running COA Framework

# Estimating PLS model on the pooled data - segmentation purpose
utaut_model <- estimate_pls(data = utaut_data,
                            measurement_model = utaut_mm,
                            structural_model = utaut_sm)

# Running COA
utaut_overfit <- coa(pls_model = utaut_model, 
                     focal_construct = "BI",
                     params = "path_coef")

# Finding viable segments with model_segments()
utaut_segments <- model_segments(utaut_overfit)

# Estimating PLS model on the pooled data with mode B - MICOM purpose
utaut_model_mode_b <- estimate_pls(data = utaut_data,
                            measurement_model = utaut_mm_mode_b,
                            structural_model = utaut_sm)


### Running MICOM

# Executing MICOM evaluation (1000 permutations)
micom_results <- evaluate_micom(original_model = utaut_model, # choose modeB or modeA
                                segments_or_condition = utaut_segments, # any other condition is accepted by segments_or_condition
                                nperms = 5000)

# Obtaining MICOM reports
micom_reports <- report_micom(micom_results = micom_results, confidence_level = 0.95)

### Running robustMICOM

robustMICOM_a_raw <- evaluate_micom_wilcox(original_model = utaut_model,
                      segments_or_condition = utaut_segments, std = FALSE)

robustMICOM_a_std <- evaluate_micom_wilcox(original_model = utaut_model,
                                           segments_or_condition = utaut_segments)

robustMICOM_b_raw <- evaluate_micom_wilcox(original_model = utaut_model_mode_b,
                                           segments_or_condition = utaut_segments, std = FALSE)

robustMICOM_b_std <- evaluate_micom_wilcox(original_model = utaut_model_mode_b,
                                           segments_or_condition = utaut_segments)


raw_robust_A <- report_micom_wilcox(robustMICOM_a_raw)
std_robust_A <- report_micom_wilcox(robustMICOM_a_std)
raw_robust_B <- report_micom_wilcox(robustMICOM_b_raw)
std_robust_B <- report_micom_wilcox(robustMICOM_b_std)


raw_robust_A
std_robust_A
raw_robust_B
std_robust_B 

### Running multiMICOM
#2-groups benchmark
results_multi_A_raw <- evaluate_multiMICOM(original_model = utaut_model, segments_or_condition = utaut_segments, std = FALSE) # even if the conditions are just two, we should write them both

results_multi_A_std <- evaluate_multiMICOM(original_model = utaut_model, segments_or_condition = utaut_segments, std = TRUE) # even if the conditions are just two, we should write them both

results_multi_B_raw <- evaluate_multiMICOM(original_model = utaut_model_mode_b, segments_or_condition = utaut_segments, std = FALSE) # even if the conditions are just two, we should write them both

results_multi_B_std <- evaluate_multiMICOM(original_model = utaut_model_mode_b, segments_or_condition = utaut_segments, std = TRUE) # even if the conditions are just two, we should write them both

raw_multi_A <- report_multiMICOM(results_multi_A_raw, significance_level = 0.05)
std_multi_A <- report_multiMICOM(results_multi_A_std, significance_level = 0.05)
raw_multi_B <- report_multiMICOM(results_multi_B_raw, significance_level = 0.05)
std_multi_B <- report_multiMICOM(results_multi_B_std, significance_level = 0.05)

#3-groups demo
results_multi3_A_std <- evaluate_multiMICOM(original_model = utaut_model, segments_or_condition = c("HAB4 < 3", " HAB4 >= 3 & HAB4 <= 5", "HAB4 > 5"))

std_multi3_A <- report_multiMICOM(results_multi3_A_std, significance_level = 0.05)
