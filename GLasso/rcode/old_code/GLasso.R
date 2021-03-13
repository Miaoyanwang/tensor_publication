GLasso = function(S,n){
  
  # this function tunes the choice of rho via BIC
  
  # Input:
  # S -- sample covariance matrix (symmetric matrix)
  # n -- sample size
  
  # Output:
  # Omega -- estimated inverse covariance matrix
  
  tmp = max(abs(S))
  tune = seq(tmp/10,tmp,tmp/10)
  BIC = rep(0,length(tune))
  
  for (i in 1:length(tune)) {
    rho = tune[i]
    Omega = graphicalLasso(S,rho)$Theta
    num = sum(upper.tri(Omega) != 0)
    BIC[i] = 2*n*(sum(diag(S%*%Omega)) - log(det(Omega))) + log(n)*num
    # trace = sum(diag())
  }# end for
  
  plot(tune,BIC,"b",xlab = "tune",ylab = "BIC")
  
  ind = order(BIC)[1]
  Omega = graphicalLasso(S, tune[ind])$Theta
  
  return(list(Omega = Omega, rho_opt = tune[ind]))
  
}


