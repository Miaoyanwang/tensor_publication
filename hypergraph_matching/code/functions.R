############### Order-3 Gaussian Tensor Matching ###############
############### Jiaxin Hu  Feb 26 22 #################

##### load package 

library(RcppHungarian) # HungarianSolver
library(gtools) # permute
library(combinat) # permn
library(MASS) # mvrnorm
#library(LaplacesDemon) # KLD
#library(seewave) # kl.dist

source("bricks.R")



###### unseeded matching ##### 

tensor_matching_unseed = function(A,B, p, test = FALSE){
  
  d = dim(A)[1]
  
  # calculate distance
  d_collect = matrix(0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (k in 1:d) {
      # i = 1; k = 15
      d_collect[i,k] = dist_p(as.vector(A[i,,]), as.vector(B[k,,]),p)
    }
  }
  
  if(test == FALSE){
    # sort distance
    S = which(d_collect <= sort(d_collect)[d], arr.ind = T)
    
    # check S is a perfect matching
    if(length(unique(S[,1])) != d | length(unique(S[,2])) != d){
      warning("No perfect matching output!")
    }
    
    S=S[order(S[,2]),]
    return(S)
  }else if(test == TRUE) {
    S = HungarianSolver(d_collect)$pairs
    return(S)
  }
}

###### seeded matching ######

# d = 50
# rho = 0.97
# p = 1
# test = sim_corr_tensor(seed = 6, d, rho)
# A = dat$A
# B = dat$B
# p = 1
# # 
#xi = 2*sqrt((log(d))^(1/2))
# #zeta = quantile(as.vector(d_collect), 0.05)
# zeta = 3*sqrt(sqrt(1 - rho^2)/d^2)

tensor_matching_seed = function(A,B,p,xi,zeta, clean_up){
  
  d = dim(A)[1]
  
  # calculate distance
  d_collect = matrix(0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (k in 1:d) {
      d_collect[i,k] = dist_p(as.vector(A[i,,]), as.vector(B[k,,]),p)
    }
  }

  
  low_dist = which(d_collect <= zeta, arr.ind = T)
  
  # high-degree seed set
  # calculate ai bk
  ai = apply(A, 1, function(x) sum(x)/sqrt(d))
  bk = apply(B, 1, function(x) sum(x)/sqrt(d))
  
  seed_set = low_dist[ low_dist[,1] %in% which(ai >= xi) & low_dist[,2] %in% which(bk >= xi), ]
  #seed_set1 =  clean_seed(seed_set,d_collect)
  count = 0
  while (! perfect_matching(seed_set) & count < 10) {
    seed_set = clean_seed(seed_set,d_collect)
    count = count + 1
    
    if(count == 10){
      warning("No clear seed.")
      return(NULL)
    }
  }
  
  if(is.null(dim(seed_set)) | dim(matrix(seed_set, ncol = 2))[1] == 0){
    warning("No seed. Choose a smaller xi or larger zeta!")
    return(NULL)
  }
  
  # seeded algorithm
  hat_pi = seeded_bipartite_matching(A,B,seed_set,clean_up)
  
  if(length(unique(hat_pi[,1])) != d | length(unique(hat_pi[,2])) != d){
    warning("No perfect matching output!")
  }
  
  hat_pi = hat_pi[order(hat_pi[,2]),]
  return(hat_pi)
}



seeded_bipartite_matching = function(A,B,seed_set, clean_up){
  
  d = dim(A)[1]
  s1_c = c(1:d)[-seed_set[,1]]
  s2_c = c(1:d)[-seed_set[,2]]
  
  # generate H
  H = matrix(0, nrow = length(s1_c), ncol =  length(s1_c))
  rownames(H) = s1_c
  colnames(H) = s2_c
  
  for (i in 1:length(s1_c)) {
    for (k in 1:length(s2_c) ) {
      
      #i = 1; k = 1
      #cat(s1_c[i],s2_c[k], "\n")
      
      A_iw = A[s1_c[i], seed_set[,1], seed_set[,1] ]
      B_kw = B[s2_c[k], seed_set[,2], seed_set[,2] ]
      
      H[i,k] = sum(A_iw*B_kw)
    }
  }
  
  # bipartite matching with H
  pi1_tilde = HungarianSolver(-H)$pairs
  pi1_tilde[,1] = s1_c[pi1_tilde[,1]]
  pi1_tilde[,2] = s2_c[pi1_tilde[,2]]
  
  pi1 = rbind(seed_set, pi1_tilde)
  
  if(clean_up == T){
    # clean up
    
    W = matrix(0, nrow = d, ncol = d)
    for (i in 1:d) {
      for (k in 1:d) {
        
        A_iw = A[i,pi1[,1],pi1[,1] ]
        B_kw = B[k, pi1[,2],pi1[,2] ]
        
        W[i,k] = sum(A_iw*B_kw)
      }
    }
    
    hat_pi = which(W >= sort(W, decreasing = T)[d], arr.ind = T)
  }else if(clean_up == F){
    hat_pi = pi1
  }
  
  
  return(hat_pi)
}


####### generator ########



sim_corr_tensor = function(seed = NA,d,rho){
  
  if (is.na(seed) == FALSE) set.seed(seed)
  
  A = array(0,dim = rep(d,3))
  B = array(0,dim = rep(d,3))
  
  for (i in 1:d) {
    for (j in i:d) {
      for (k in j:d) {
        
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
  pi_star = permute(1:d)
  
  B = B[pi_star, pi_star, pi_star]
  pi = cbind(pi_star,c(1:d))
  
  return(list(A = A, B = B, pi = pi))
}
