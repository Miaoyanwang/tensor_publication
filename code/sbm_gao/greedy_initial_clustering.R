## initialization in sbm
library(pracma)

source("/Users/March/Desktop/myresearch/code/sbm_gao/bricks_sbm.R")

greedy_initial_clustering = function(A,k,trim = TRUE, tau = Inf,mu){
  # USC
  #A is adjacency
  #k is # of clusters
  #tau is trimming para
  #r is cluster radius
  
  n = dim(A)[1]
  r = mu*sqrt(k/n)
  
  if(trim == TRUE){
    dvec = colSums(t(A))
    index_trim = dvec >= tau
    A[index_trim,] = 0
    A[,index_trim] = 0
  }
  
  U = eigen(A)$vectors[,1:k]
  
  # U = matrix(1:20,ncol = 2)
  # n = 10
  # k = 3
  # r = 5*sqrt(2/10)
  
  S = 1:n
  
  est_C = rep(0,n)
  # generate norm matrix
  norm_mat = gene_norm_mat(U)
  norm_mat0 = as.matrix(norm_mat)
  
  for (i in 1:k) {
    
    index_norm_mat = norm_mat0 < r
    if(dim(index_norm_mat)[1] == 1) {rsum = norm_mat0
    } else {rsum = colSums(index_norm_mat)}
    ti = which(rsum == max(rsum))[1]
    Ci = which((index_norm_mat[,ti]) == TRUE)
    est_C[S[Ci]]= i
    
    S = S[-Ci]
    
    if(i < k & length(S) == 0){
      print("please take smaller mu")
      return(est_C)
    }
    
    norm_mat0 = norm_mat0[-Ci,]
    if(class(norm_mat0) == "numeric"){norm_mat0 = matrix(norm_mat0,nrow = 1)}
    norm_mat0 = (norm_mat0)[,-Ci]
  }
  
  est_C
  
  if(length(S) != 0){
    for (i in 1:length(S)) {
      usum = rep(0,k)
      for (j in 1:k) {
        index_c = which(est_C == j)
        usum[j] = sum(norm_mat[index_c,S[i]])/length(index_c)
      }
      est_C[S[i]] = which(usum == min(usum))
    }
  }
  est_C
  return(est_C)
}


