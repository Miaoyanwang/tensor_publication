### bricks for sbm
ReNumber = function (Cs) 
{
  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  for (i in 1:length(uniq)) {
    newCs[Cs == uniq[i]] <- i
  }
  return(newCs)
}


norm_cal = function(v1,v2){
  return(Norm(v1-v2))
}

gene_norm_mat = function(U){
  # U \in R^n \times k
  # cal every || Ui. - Uj. ||
  
  #U = matrix(1:10,nrow = 5,ncol = 2)
  n = dim(U)[1]
  
  norm_mat = matrix(0,nrow = n,ncol = n)
  for (i in 1:(n-2)) {
    norm_mat[(i+1):n,i] = apply(U[(i+1):n,], 1, norm_cal,v2 = U[i,])
  }
  norm_mat[n,(n-1)] = norm_cal(U[n,],U[(n-1)])
  
  norm_mat[upper.tri(norm_mat)] = t(norm_mat)[upper.tri(norm_mat)]
  #norm_mat is a symmetric matrix
  return(norm_mat)
}
