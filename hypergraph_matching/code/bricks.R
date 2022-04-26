### bricks functions.R
### Apr 25 22 Jiaxin Hu

###### Sub-functions for main functions
# calculate L distance with partition over [-1,1]

Ldist = function(X, Y, L){
  
  cdf_x = ecdf(X)
  cdf_y = ecdf(Y)
  
  # uniform L-partition over [-1, 1]
  knots = seq(min(c(X,Y)), max(c(X,Y)), length.out = (L+1))

  # calculate L dist
  Fn = cdf_x(knots[2:(L+1)]) - cdf_x(knots[1:L])
  Gn = cdf_y(knots[2:(L+1)]) - cdf_y(knots[1:L])
  L_dist = sum(abs(Fn - Gn))
  # 
  # knots = c(X, Y)
  # Fn = cdf_x(knots[2:length(knots)]) - cdf_x(knots[1:(length(knots) - 1)])
  # Gn = cdf_y(knots[2:length(knots)]) - cdf_y(knots[1:(length(knots) - 1)])
  # L_dist = sum(abs(Fn - Gn))
  
  return(L_dist = L_dist)
}

# sup-norm 

supnorm = function(X,Y){
  
  cdf_x = ecdf(X)
  cdf_y = ecdf(Y)
  
  knots = c(X, Y)
  
  sup_norm = max(abs(cdf_x(knots) -cdf_y(knots) ))
 
  
  return(sup_norm = sup_norm)
}

# 1-Wasserstein distance

W1dist = function(X,Y){
  
  W1_dist = sum(abs(sort(X) - sort(Y)))/length(X)
  
  return(W1_dist = W1_dist)
}

####### Tools 

# measure the error between permutation

perm_error = function(m1, m2){
  
  # m1, m2 are two permutation of the same dimension n
  # m1 = (pi1(1), ... ,pi1(n)); m2 = (pi2(1), ... ,pi2(n)); 
  
  err = sum(m1 != m2)/length(m1)
  return(err)
}


