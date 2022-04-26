############### Order-3 Gaussian Tensor Matching ###############
############### Jiaxin Hu  Apr 25 22 #################

########### Previous version: functions_old.R, bricks_old.R
########### Major revision for notation consistency and algorithm modification


##### Load dependencies

library(RcppHungarian) # HungarianSolver
library(gtools) # permute
library(combinat) # permn
library(MASS) # mvrnorm

source("bricks.R")


####### Unseeded algorithm 

tensor_matching_unseed = function(A, B, dist = c("Ldist", "supnorm", "W1"), L = NULL){
  
  n = dim(A)[1]
  
  # calculate the distance matrix D
  d_collect = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (k in 1:n) {
      
      if(dist == "Ldist"){
        d_collect[i,k] = Ldist(as.vector(A[i,,]), as.vector(B[k,,]), L)
      }else if(dist == "supnorm"){
        d_collect[i,k] = supnorm(as.vector(A[i,,]), as.vector(B[k,,]))
      }else if(dist == "W1"){
        d_collect[i,k] = W1dist(as.vector(A[i,,]), as.vector(B[k,,]))
      }
      
    }# end k
  }# end i
  
  # Hungrarian algorithm for linear assignment 
  S = HungarianSolver(d_collect)$pairs
  pi_hat = S[order(S[,2]), 1]
  
  return(pi_hat)
}



######## Generator 

sim_corr_tensor = function(seed = NA,n,rho){
  
  if (is.na(seed) == FALSE) set.seed(seed)
  
  A = array(0,dim = rep(n,3))
  B = array(0,dim = rep(n,3))
  
  for (i in 1:n) {
    for (j in i:n) {
      for (k in j:n) {
        
        # correlated data for true pair
        #pair = mvrnorm(n = 1, mu = rep(0,2), Sigma = matrix(c(0.5,rho*0.5,rho*0.5,0.5), nrow = 2))
        pair = mvrnorm(n = 1, mu = rep(0,2), Sigma = matrix(c(1,rho,rho,1), nrow = 2))
        
        # generate super-symmetric pair
        #perm = permutations(n = 3, r = 3, v = c(i,j,k), repeats.allowed = TRUE)
        perm = unique(permn(c(i,j,k)))
        perm = matrix(unlist(perm), ncol= 3,byrow = T)
        
        A[perm] = pair[1]
        B[perm] = pair[2]
      } # i
    } # j
  } # k 
  
  # generate permutation
  pi_star = permute(1:n)
  
  B = B[pi_star, pi_star, pi_star]
  pi_star_pair= cbind(pi_star,c(1:n))
  
  return(list(A = A, B = B, pi_star = pi_star, pi_star_pair = pi_star_pair))
}
