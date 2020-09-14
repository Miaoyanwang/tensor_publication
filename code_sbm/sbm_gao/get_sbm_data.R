## get data function for SBM in Gao

source("/Users/March/Desktop/myresearch/code/sbm_gao/bricks_sbm.R")

get_sbm_data = function(n,k,a,b){
  # n is the matrix demsion
  # k is the # of clusters
  # a is the signal level of diagonal elements 
  # b is the signal level of off-diagonal elements
  
  # n = 50; k = 4
  # a = 25; b = 25
  
  prob = matrix(runif(k^2),nrow = k)
  # symmetric
  prob[lower.tri(prob)] = t(prob)[lower.tri(prob)]  
  # signal contrain
  # min Bii > a/n; max Bij < b/n
  diag(prob)[which(diag(prob) <= a/n)] = a/n
  index_mat = prob >= b/n
  diag(index_mat) = FALSE
  prob[index_mat]  = b/n 
  
  if (k == n) truthC = 1:n
  else if (k == 1) truthC = rep(1,n)
  else {truthC = sort(ReNumber(sample(1:k,n,replace=TRUE)) )}
  
  x = matrix(1,nrow = n,ncol = n)
  truthX =  matrix(0,nrow = n,ncol = n)
  for(i in 1:max(truthC)){
    for(j in 1:max(truthC)){
        x[truthC==i, truthC==j] = rbinom(sum(x[truthC==i, truthC==j]),1,prob[i,j])
        truthX[truthC==i, truthC==j] =  prob[i,j]
    }
  }
  
  result = list("x"=x,"truthX"=truthX,"truthC"=truthC,"prob"=prob)
  # x is the observation
  # truthX is true connection
  # truthC is true partition
  # prob is true adjacency
  return(result)
}


