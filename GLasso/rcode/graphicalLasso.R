library(glmnet)

graphicalLasso = function(S, rho, maxIt = 1e2, tol = 1e-6){
  
  # solving the graphical lasso 
  # problem (2.1): 
  # min_{Theta > 0} tr(S*Theta) - log det(Theta) + rho*||Theta||_1
  # or the equivalent problem (2.4):
  # min_{beta} || W_11^1/2 beta - W_11^-1/2 s_12 ||^2 + rho*||beta||_1
  
  # need a Lasso algorithm in the penalized form
  # Here we use the function glmnet() in the package glmnet
  
  # Input
  # S -- sample covariance matrix (symmetric matrix)
  # rho -- the given regularization parameter
  # maxIt -- maximum number of iterations
  # tol -- convergence tolerance level
  
  # Output
  # Theta -- estimated inverse covariance matrix
  # W -- estimated regularized covariance matrix
  # Theta = W^-1
  
  p = dim(S)[1] 
  
  # initialization
  W = S + rho*diag(rep(1,p))
  W_old = W
  
  i = 0
  
  # Graphical Lasso loop
  while (i < maxIt) {
    i = i+1
    
    for (j in p:1) {
      #jminus = setdiff(1:p,j)
      eigen_W11 = eigen(W[-j,-j],symmetric = T)
      V = eigen_W11$vectors
      D = eigen_W11$values
      
      X = V%*%diag(sqrt(D))%*%t(V) # W_11^1/2, W is symmetric, thus V is orthonormal
      Y = V%*%diag(1/sqrt(D))%*%t(V)%*%S[-j,j] # W_11^{1/2}*s_12
      
      #Apply Lasso algorithm to solve beta
      fit_lasso = glmnet(X,Y,family = "gaussian", alpha = 1, lambda = rho, intercept = FALSE,
                    thresh = tol,maxit = maxIt)
      beta = as.matrix(coef(fit_lasso)[-1])
      
      # plug w12 in W
      W[-j,j] = W[-j,-j]%*%beta
      W[j,-j] = t(W[-j,j])
      
      #print(norm(W-W_old,"1"))
    }# end for
  
    # stop criterion
    if( norm(W-W_old,"1") < tol){
      break
    }
    
    W_old = W
  }#end while
  
  if(i == maxIt){
    cat("Maximum number of iteration reached, glasso may not converge.")
  }
  
  Theta = solve(W)
  
  return(list(W = W, Theta = Theta))
}


