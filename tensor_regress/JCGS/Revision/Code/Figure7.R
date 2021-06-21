# Code for Figure 7. HCP data analysis

# load dependencies ---------------------------------------------------------------------------------------------
install.packages("software/tensorregress.tar.gz", repos = NULL, type = "source")
library(tensorregress)

# load data -----------------------------------------------------------------------------------------------------
data <- load("rawdata/HCP.RData")

# start analysis ------------------------------------------------------------------------------------------------
tsr <- tensor
X <- attr[, 4:5]
# three groups
levels(X[, 2]) <- c("22-25", "26-30", "31+", "31+")
contrasts(X[, 1]) <- contr.sum
contrasts(X[, 2]) <- contr.sum
X_covar3 <- model.matrix(~ as.factor(X[, 1]) + as.factor(X[, 2]))

set.seed(0)
core_shape <- c(10, 10, 4)
# The following line takes quite a while to finish. We do not suggest to run on personal laptop. We implement it in the server.
result <- tensor_regress(tsr, X_covar1 = NULL, X_covar2 = NULL, X_covar3, core_shape, niter = 50, cons = "non", lambda = 0.1, alpha = 10, solver = "CG", dist = "binary", initial = "random")

save(result, file = "output_HCP.RData")
# the matrix for plotting brain connectivity is saved in presaved/data_for_Figure7
