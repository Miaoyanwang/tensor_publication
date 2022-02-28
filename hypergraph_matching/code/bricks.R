### bricks  for Gaussian tensor matching
### Feb 26

dist_p = function(Ai,Bi,p){
  
  if(p == 1){
    # 1-Wasserstein 
    # Ai, Bi are two vectors of dim d^2
    dp = sum(abs(sort(Ai) - sort(Bi)))
  }
  if(p == 2){
    cdf1 = ecdf(Ai)
    cdf2 = ecdf(Bi)
    
    knots = c(Ai,Bi)
    
    dp = sqrt(sum( (cdf1(knots) - cdf2(knots))^2 ))
  }
  
  return(dp)
}

match_error = function(m1, m2){
  # m1, m2 two-column matrices of the same dimension
  # m2 true
  
  d = dim(m2)[1]
  
  if(length(unique(m2[,1]))!= d |length(unique(m2[,2]))!= d){
    warning("m2 should be the true one-to-one matching!" )
    return()
  }
  
  err = 0
  for (i in 1:d) {
    
    tm = m2[i,1]
    m1_cut = m1[m1[,2] == i, 1]
    
    err = err + sum(m1_cut!= tm)
  }
  
  err = err/d
  
  # m1 = m1[order(m1[,2]),]
  # m2 = m2[order(m2[,2]),]
  
  return(err)
}


