# Install cSEM
install.packages("cSEM")
require(cSEM)

# Models - MM (Composite) and SM
model <- "
# Measurement model 
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
  
# Structural Model
BI ~ PE + EE + SI + FC + HM + PV + HAB + Exp + Age + Gender
"

# Importing data
g1 <- read.csv(file = "data/trello_utaut.csv")[,-66]
g2 <- read.csv(file = "data/segment2.csv")[,-66]
g3 <- read.csv(file = "data/segment3.csv")[,-66]


# Estimating model
utaut_cSEM_3g <- csem(.data = list("group1" = g1, "group2" = g2, "group3" = g3), .model = model)

# Verifying groups' admissibility (admissibility of results)
verify(utaut_cSEM$group1)
verify(utaut_cSEM$group2)
verify(utaut_cSEM$group3)

# MICOM
MICOM_3g <- testMICOM(.object = utaut_cSEM_3g, .R = 1000, .seed = 123)
