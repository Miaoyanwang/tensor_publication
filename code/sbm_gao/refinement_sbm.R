### refinement in sbm
source("/Users/March/Desktop/myresearch/code/sbm_gao/bricks_sbm.R")

library(matrixcalc)


refinement_sbm = function(A, k, trim = TRUE, tau = Inf,mu){
  # A is adjacency
  # k is # of clusters
  # trim, tau, mu is for initialization
  
  # A = test_data$x
  # k = 3
  
  n = dim(A)[1]
  
  sigma_t = matrix(0,nrow = n,ncol = n)
  
  for (i in 1:n) {
    
    #i = 1
    
    #get A-u
    A_sub = A[-i,]; A_sub = A_sub[,-i]
    
    
    sigmau = rep(0,n)
    #need to change
    sigmau[-i] = greedy_initial_clustering(A_sub,k,trim,tau,mu)
    
    B_est = matrix(0,nrow = k,ncol = k)
    for (j in 1:(k-1)) {
      
      Cj = which(sigmau == j)
      B_est[j,j] = sum(A[Cj,Cj])/(2*length(Cj)*(length(Cj) - 1))
      
      for (l in (j+1):k) {
        Cl = which(sigmau == l)
        B_est[l,j] = 2*sum(A[Cl,Cj])/(length(Cj)*length(Cl))
      }
    }
    
    Ck = which(sigmau == k)
    B_est[k,k] = sum(A[Ck,Ck])/(2*length(Ck)*(length(Ck) - 1))
    
    # now B_est is a lower triangle matrix
    #B_est[upper.tri(B_est)] = t(B_est)[upper.tri(B_est)]
    
    a_est = n*min(diag(B_est))
    #B_off_est = B_est; diag(B_off_est) = 0
    b_est = n*max(B_est[lower.tri(B_est)])
    
    #calculate rho
    tu = (1/2)*log((a_est*(1-b_est/n))/(b_est*(1-a_est/n))) 
    rhou = (-1/(2*tu))*log( (a_est/n*exp(-tu) + 1 - a_est/n  ) / (b_est/n*exp(tu) + 1 - b_est/n ) )
    
    usum = rep(0,k)
    for (j in 1:k) {
      Cj = which(sigmau == j)
      usum[j] = sum(A[i,Cj]) - rhou*length(Cj)
    }
    
    sigmau[i] = which(usum == max(usum))
    
    sigma_t[i,] = sigmau
  }
  
  #consensus
  sigma_est = rep(0,n)
  sigma_est[1] = sigma_t[1,1]

  for (i in 2:n) {
    Cu = which(sigma_t[i,] == sigma_t[i,i])
    inter = sigma_t[1,][Cu]
    sigma_est[i] = Mode(inter)
  }
  
  return(sigma_est)
}

