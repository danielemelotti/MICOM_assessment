### Comparing outer_weights produced by SEMinR vs smartPLS

# Importing smartPLS weights for pooled data, segment 2, segment 3
pooled_reflective_smartpls_weights <- read.csv(file = "data/smartpls/pooled_reflective_smartpls_weights.csv", row.names = 1)[item_names, construct_names]

pooled_reflective_smartpls_weights[is.na(pooled_reflective_smartpls_weights)] <- 0

seg2_reflective_smartpls_weights <- read.csv(file = "data/smartpls/seg2_reflective_smartpls_weights.csv", row.names = 1)[item_names, construct_names]

seg2_reflective_smartpls_weights[is.na(seg2_reflective_smartpls_weights)] <- 0

seg3_reflective_smartpls_weights <- read.csv(file = "data/smartpls/seg3_reflective_smartpls_weights.csv", row.names = 1)[item_names, construct_names]

seg3_reflective_smartpls_weights[is.na(seg3_reflective_smartpls_weights)] <- 0

pooled_reflective_smartpls_weights
seg2_reflective_smartpls_weights
seg3_reflective_smartpls_weights

# SEMinR results for pooled data, segment 2, segment 3
segment_labels <- data_manipulation(utaut_data, utaut_model, utaut_segments) #obtaining labels to divide dataset #####################

labeled_data <- cbind(segment_labels$segment_labels, utaut_data)
colnames(labeled_data)[1] <- "indicator" 

seg2_data <- labeled_data[labeled_data$indicator == "segment_2", item_names]
seg3_data <- labeled_data[labeled_data$indicator == "segment_3", item_names]

seg2_seminr <- rerun(utaut_model, data = seg2_data)
seg3_seminr <- rerun(utaut_model, data = seg3_data)

# Eventually, the weights are exactly the same as smartPLS'
pooled_seminr_weights <- as.data.frame(utaut_model$outer_weights)
seg2_seminr_weights <- as.data.frame(seg2_seminr$outer_weights)
seg3_seminr_weights <- as.data.frame(seg3_seminr$outer_weights)



### Tests to find out the incongruous observed_c between our technique and smartPLS

# Extracting the FC weights from seg2 and seg3
seg2_seminr_FC_weights <- as.data.frame(seg2_seminr_weights[c("FACL1", "FACL2", "FACL3", "FACL4"), "FC"])
seg3_seminr_FC_weights <- as.data.frame(seg3_seminr_weights[c("FACL1", "FACL2", "FACL3", "FACL4"), "FC"])

# Adjusting the data to allow the composite computation -> reduce to only having 
pooled_FC_data <- as.matrix(utaut_data[, c("FACL1", "FACL2", "FACL3", "FACL4")])

# Standardize the data
standardized_pooled_FC_data <- scale(pooled_FC_data)

# Computing composite scores
seg2_FC_scores <- pooled_FC_data %*% as.matrix(seg2_seminr_FC_weights)
seg3_FC_scores <- pooled_FC_data %*% as.matrix(seg3_seminr_FC_weights)

# And computing composite scores using standardized data
standardized_seg2_FC_scores <- standardized_pooled_FC_data %*% as.matrix(seg2_seminr_FC_weights)
standardized_seg3_FC_scores <- standardized_pooled_FC_data %*% as.matrix(seg3_seminr_FC_weights)

# Checking observed_c
cor(seg2_FC_scores, seg3_FC_scores) # c = 0.1717, no change if we use regular data...
cor(standardized_seg2_FC_scores, standardized_seg3_FC_scores) # c = 0.30789 if we use standardized data!!!

# Trying the whole procedure on the whole dataset and all constructs
standardized_pooled_test <- scale(utaut_data, center = TRUE, scale = TRUE) # standardize the whole data

seg_2_test_scores <- standardized_pooled_test %*% as.matrix(seg2_seminr_weights) 
seg_3_test_scores <- standardized_pooled_test %*% as.matrix(seg3_seminr_weights)

standardized_test_results <- round(diag(cor(seg_2_test_scores, seg_3_test_scores)), 6) #0.308
our_micom <- micom_reports$`Step 2: Compositional Invariance`[, 2] #0.1717
smart_micom<- smartpls_micom$`Step 2: Compositional Invariance`[, 1] #0.308

test_results_comparison <- as.data.frame(cbind(standardized_test_results,
                                               our_micom,
                                               smart_micom))

# For correlations, the problem seems solved, but not for means and variance ratios
colMeans(seg_2_test_scores) - colMeans(seg_3_test_scores)


### Comparing MICOM results produced by evaluate_micom() vs smartPLS

# Define function to import data
import_smartPLS_data <- function(mode, type) {
  cors <- read.csv(paste0("data/smartpls/micom_cors_", mode, type, ".csv"), row.names = 1)[construct_names, ]
  means <- read.csv(paste0("data/smartpls/micom_means_", mode, type, ".csv"), row.names = 1)[construct_names, ]
  vars <- read.csv(paste0("data/smartpls/micom_vars_", mode, type, ".csv"), row.names = 1)[construct_names, ]
  set <- list(cors, means, vars)
  
  names(set) <- c("Step 2: Compositional Invariance", "Step 3.1: Equality of Mean Values", "Step 3.2: Equality of Variances")
  set
}

# Importing ALL smartPLS micom results
smartpls_micom_modeA_raw <- import_smartPLS_data("modeA_", "raw")
smartpls_micom_modeA_std <- import_smartPLS_data("modeA_", "std")
smartpls_micom_modeB_raw <- import_smartPLS_data("modeB_", "raw")
smartpls_micom_modeB_std <- import_smartPLS_data("modeB_", "std")

### smartPLS results
smartpls_micom_modeA_raw 
smartpls_micom_modeA_std 
smartpls_micom_modeB_raw 
smartpls_micom_modeB_std 

### Our micom results
modeA_raw 
modeA_std
modeB_raw 
modeB_std 

### cSEM results
csem_micom_modea_raw
csem_micom_modea_std
csem_micom_modeb_raw
csem_micom_modeb_std
