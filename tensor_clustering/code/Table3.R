# Code for Table 3. Peru Legislation data analysis

# load dependencies & data ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)
library(RSKC)

# start analysis ------------------------------------------------------------------------------------------------
result = dtbm(peru$network_data, r = rep(5,3), max_iter = 20, asymm = F)