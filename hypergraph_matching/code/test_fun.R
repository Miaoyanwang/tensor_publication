###### tests function for the function 

##### Load dependencies

library(RcppHungarian) # HungarianSolver
library(gtools) # permute
library(combinat) # permn
library(MASS) # mvrnorm

source("functions.R")

# dist = "Ldist"
# dist = "supnorm"
# dist = "W1"

# small test

n = 30
test_dat = sim_corr_tensor(n = n, rho = 0.995)
test_dat$pi_star
# test_dat$pi_star_pair
# L = ceiling(3*log(n)) 
L = 15

A = test_dat$A
B = test_dat$B
summary(as.vector(B))
test_res = tensor_matching_unseed(A,B, dist = "Ldist", L)
perm_error(test_res, test_dat$pi_star)

test_res = tensor_matching_unseed(A,B, dist = "supnorm", L)
perm_error(test_res, test_dat$pi_star)

test_res = tensor_matching_unseed(A,B, dist = "W1", L)
perm_error(test_res, test_dat$pi_star)

