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

## cal loss

gene_conf_mat_one = function(r,true_assign, est_assign){
  #assigns are vectors
  
  conf_mat = matrix(0,nrow = r,ncol = r)
  d = length(true_assign)
  
  for (i in 1:d) {
    
    x1 = true_assign[i] 
    x2 = est_assign[i] 
    conf_mat[x1,x2] = conf_mat[x1,x2] + 1
    
  }
  
  return(conf_mat)
}

cal_loss_sbm = function(k,true_assign,est_assign){
  conf_mat = gene_conf_mat_one(k,true_assign,est_assign)
  #conf_mat =  matrix(c(10,8,11,1,5,4,2,9,5,7,4,9,5,6,8,6),nrow = 4)
  
  conf_mat0 = conf_mat
  S = 1:k
  for (i in 1:(k-1)) {
    
    index_max = which(conf_mat0 == max(conf_mat0),arr.ind = T)[1,]
    
    if(index_max[1] != index_max[2]){
      indr = S[index_max[1]]; indc = S[index_max[2]]
      #swap colmuns indr and indc
      est_col = conf_mat[,indc]
      conf_mat[,indc] = conf_mat[,indr]
      conf_mat[,indr] = est_col
      
      #also swap conf_mat0
      est_col0 = conf_mat0[,index_max[2]]
      conf_mat0[,index_max[2]] = conf_mat0[,index_max[1]]
      conf_mat0[,index_max[1]] = est_col0
    }
    
    if(i == (k-1)) break
    
    conf_mat0 = conf_mat0[-index_max[1],]
    
    #if(class(conf_mat0) == "numeric") conf_mat0 = matrix(conf_mat0,nrow = 1)
    conf_mat0 = conf_mat0[,-index_max[1]]
    S = S[-index_max[1]]
  }
  
  n = length(true_assign) 
  MCR_sbm = (sum(conf_mat[lower.tri(conf_mat)]) + sum(conf_mat[upper.tri(conf_mat)]) )/n
  
  return(MCR_sbm)
}
