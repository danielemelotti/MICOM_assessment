### Performing MICOM using csem package

#install.packages("cSEM")
require(cSEM)

# My two segments
seg2_data
seg3_data

# Making the model
lavaan_model <- '
    # structural
    BI ~ PE + EE + SI + FC + HM + PV + HAB + Exp + Age + Gender
    
    # measurement
    PE <~ PERF1 + PERF2 + PERF3 + PERF4
    EE <~ PEOU1 + PEOU3 + PEOU5 + PEOU6
    SI <~ NORM1 + NORM2 + INFL3
    FC <~ FACL1 + FACL2 + FACL3 + FACL4
    HM <~ MOTIV1 + MOTIV2 + MOTIV3
    PV <~ VALUE1 + VALUE2 + VALUE3
    HAB <~ HAB1 + HAB2 + HAB3 + HAB4
    BI <~ INT1 + INT2 + INT3
    Exp <~ Experience
    Age <~ age
    Gender <~ gender
'

# Random seed
set.seed(121514)

## Raw

# Making a list of the 2 datasets
dat <- list(seg2_data, seg3_data)

#Running csem
csem_modeb <- csem(dat, lavaan_model, .resample_method = "bootstrap", .PLS_modes = "modeB")
csem_modea <- csem(dat, lavaan_model, .resample_method = "bootstrap", .PLS_modes = "modeA")

# MICOM
csem_micom_modeb_raw <- testMICOM(csem_modeb, .R = 5000)
csem_micom_modea_raw <- testMICOM(csem_modea, .R = 5000)


## Standardized

standardized_utaut_data <- as.data.frame(utaut_data)
standardized_utaut_data[ , names(standardized_utaut_data) != "gender"] <- 
  scale(utaut_data[ , names(utaut_data) != "gender"])

segment_labels_std <- data_manipulation(utaut_data, utaut_model, utaut_segments) #obtaining labels to divide dataset

labeled_data_std <- cbind(segment_labels_std$segment_labels, standardized_utaut_data)
colnames(labeled_data_std)[1] <- "indicator" 

seg2_data_std <- labeled_data_std[labeled_data_std$indicator == "segment_2", item_names]
seg3_data_std <- labeled_data_std[labeled_data_std$indicator == "segment_3", item_names]

stand_dat <- list(seg2_data_std, seg3_data_std)

csem_modea_std <- csem(stand_dat, lavaan_model, .resample_method = "bootstrap", .PLS_modes = "modeA")
csem_modeb_std <- csem(stand_dat, lavaan_model, .resample_method = "bootstrap", .PLS_modes = "modeB")

csem_micom_modea_std <- testMICOM(csem_modea_std, .R = 5000)
csem_micom_modeb_std <- testMICOM(csem_modeb_std, .R = 5000)
