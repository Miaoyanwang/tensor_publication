#### Algorithm for SBM in Gao ####

library(plot.matrix)

source("/Users/March/Desktop/myresearch/code/sbm_gao/bricks_sbm.R")
source("/Users/March/Desktop/myresearch/code/sbm_gao/get_sbm_data.R")
source("/Users/March/Desktop/myresearch/code/sbm_gao/greedy_initial_clustering.R")
source("/Users/March/Desktop/myresearch/code/sbm_gao/refinement_sbm.R")


## generate data
test_data = get_sbm_data(50,3,15,5)
# plot(test_data$x)
# plot(test_data$truthX)
# test_data$truthC

## initialization clustering \sigma^0
test_ini = greedy_initial_clustering(test_data$x,3,mu = 0.5)
test_ini

svd(test_data$truthX)$d
# when lambda k of truth x is larger, the test ini is more accurate.

## refinement
test_refine = refinement_sbm(test_data$x,3,mu = 0.5)
test_refine
test_data$truthC
