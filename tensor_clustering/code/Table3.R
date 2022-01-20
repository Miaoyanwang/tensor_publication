# Code for Table 3. Peru Legislation data analysis

# load dependencies & data ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)
library(RSKC)

true_z = factor(peru$attr_data$party)
levels(true_z) = 1:5

# start analysis ------------------------------------------------------------------------------------------------

result = dtbm(peru$network_data, r = rep(5,3), max_iter = 20, asymm = F)
CER(result$z[[1]], true_z)
