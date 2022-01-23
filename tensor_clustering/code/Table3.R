# Code for Table 3. Peru Legislation data analysis

# load dependencies & data ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)

source("software/HOLloyd.R")
source("software/tensor_score.R")
source("software/hosvd_related.R")

library(RSKC)

# true assignment
true_z <- factor(peru$attr_data$party)
levels(true_z) <- 1:5

# start analysis ------------------------------------------------------------------------------------------------

r <- 5

# dTBM result
result <- dtbm(peru$network_data, r = rep(r, 3), max_iter = 20, asymm = F)
CER(result$z[[1]], true_z)

set.seed(10)
# HLloyd result
z.HOSC <- HO.SC(as.tensor(peru$network_data), r)
z.Lloyd.HOSC <- HO.Lloyd(as.tensor(peru$network_data), z.HOSC)
min(CER(z.Lloyd.HOSC[[1]], true_z), CER(z.Lloyd.HOSC[[2]], true_z), CER(z.Lloyd.HOSC[[3]], true_z))

# SCORE result
ke_result <- tensor_score_adj(peru$network_data, r, rm_diag = T, hooi = T, reg_hooi = T, score_reg = T, max_iter = 20)
CER(ke_result, true_z)

set.seed(50)
# HOSVD result
hosvd1 <- HOSVD_kmeans(peru$network_data, r)
CER(hosvd1, true_z)

# HOSVD+ result
hosvd2 <- HOSVD_score(peru$network_data, r)
CER(hosvd2, true_z)
