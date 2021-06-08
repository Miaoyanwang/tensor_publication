####### Useful functions #####
library(rTensor)
library(pracma)
source("netglm_code_adjust.R") 

# logistic = function(x){
#   return(exp(x)/(1 + exp(x)))
# }

######## generate data with orthogonal columns ########

tensor_trans = function(tsr){
  # Y \in R^{m x l x n} --> Y \in R^{n x m x l}
  
  mode = dim(tsr)
  new_tsr = array(0,c(mode[3], mode[1:2]))
  for (i in 1:mode[3]) {
    new_tsr[i,,] = tsr[,,i]
  }
  return(new_tsr)
}

####### grid search for GLSNet ########

search_GLM = function(data, X, true_U,rank_range, sparsity_range, dist){
  nr = length(rank_range); ns = length(sparsity_range)
  cor_mreg = err_mreg = array(0, dim = c(nr, ns))
  
  for (r in 1:nr) { # rank
    for (s in 1:ns) { # sparisty
      
      R = rank_range[r]; sparsity = sparsity_range[s]
      
      cat("rank = ", R, "sparsity = ", sparsity, "\n")
      if(dist == "normal"){
        mreg_res = mySymGLM_normal(data, X,  R, sparsity, niter=50)
        
        Theta = mreg_res$A%*%diag(mreg_res$w)%*%t(mreg_res$A)
        mreg_B = ttm(mreg_res$BB,X,m = 3)@data
        mreg_U = mreg_B
        for (k in 1:dim(mreg_B)[3]) {
          mreg_U[,,k] = mreg_B[,,k] + Theta
        }
        err_mreg[r,s] =  mean((mreg_U - true_U)^2)
        cor_mreg[r,s] = cor(mreg_U, true_U)
      }else if(dist == "binary"){
        mreg_res = mySymGLM(data, X,  R, sparsity, niter=50)
        
        Theta = mreg_res$A%*%diag(mreg_res$w)%*%t(mreg_res$A)
        mreg_B = ttm(mreg_res$BB,X,m = 3)@data
        mreg_U = mreg_B
        for (k in 1:dim(mreg_B)[3]) {
          mreg_U[,,k] = mreg_B[,,k] + Theta
        }
        err_mreg[r,s] =  mean((logistic(mreg_U) - logistic(true_U))^2)
        cor_mreg[r,s] = cor(as.vector(logistic(mreg_U)), as.vector(logistic(true_U)))
      }
      
    } # end s
  }# end r

  err_mreg[which(is.na(err_mreg))] = 0
  cor_mreg[which(is.na(cor_mreg))] = 0
  
  err_index = which(err_mreg == min(err_mreg), arr.ind = T)
  cor_index = which(cor_mreg == max(cor_mreg), arr.ind = T)

  r_opt = rank_range[err_index[1]]; s_opt = sparsity_range[err_index[2]]
  r_opt_cor = rank_range[cor_index[1]]; s_opt_cor = sparsity_range[cor_index[2]]
  
  cat("cor: rank optimal = ", r_opt_cor, "sparsity optimal = ",s_opt_cor,"\n")
  
  return(list(para_opt = data.frame(r = r_opt,s = s_opt),para_opt_cor = data.frame(r = r_opt_cor,s = s_opt_cor),
            err_opt = min(err_mreg), cor_opt = max(cor_mreg)))
}


env_sim = function(D, u, p, n, signal_level, cor_level, dup){
  # d = c(d1,d2) dim of matrix response and sample size 
  # r = c(r1, r2) envelope dim
  # p feature dim 
  # n sample size, n = d3
  
  # signal_level: max abs value of core tensor
  # cor_level: max abs value of A to generate Omega = A%*%t(A), 
  # the larger the cor_level is, the covariance matrix is futher away from diagonal matrix
  
  #D = c(20,20); u = c(3,3); p = 5; n = 20; signal_level = 3; cor_level = 0; dup = 1
  
  C = as.tensor(array(runif(u[1]*u[2]*p, min = -signal_level, max = signal_level), dim = c(u[1],u[2],p))) # core tensor
  
  gamma = list(); gamma0 = list() # orthogonal factor matrices
  Omega = list(); Omega0 = list(); # list of Omega
  Sigma = list() # list of covariance matrices on different modes, Sigma = Gamma Omega Gamma^T
  Sigmasqrtm = list() # Sigma^{1/2}
  
  for (i in 1:2) {
    
    d = D[i]; r = u[i]
    
    g_othro = randortho(d, "orthonormal")
    gamma[[i]] = g_othro[,1:r]; gamma0[[i]] = g_othro[,(r+1):d]
    
    A = array(runif(r*r, min = -cor_level, max = cor_level), dim = c(r,r))
    Omega[[i]] = A %*% t(A)
    
    A = array(runif((d-r)*(d-r), min = -cor_level, max = cor_level), dim = c(d-r,d-r))
    Omega0[[i]] = A %*% t(A)
    
    Sigma[[i]] = gamma[[i]] %*% Omega[[i]] %*% t(gamma[[i]]) + gamma0[[i]] %*% Omega0[[i]] %*% t(gamma0[[i]]) 
    Sigma[[i]] = round(Sigma[[i]], digits = 3) + diag(rep(1,d)) # sigma is non identity matrix
    
    Sigmasqrtm[[i]] = sqrtm(Sigma[[i]])$B
  }
  
  B = ttl(C, gamma, c(1,2)) # coefficient tensor
  X = array(rnorm(n*p),dim = c(n,p)) # covariates
  U = ttm(B, X, 3) 
  
  Y = list()
  for (k in 1:dup) {
    # generate noise
    Epsilon = as.tensor(array(rnorm(D[1]*D[2]*n), dim = c(D[1],D[2],n)))
    Epsilon = ttl(Epsilon, Sigmasqrtm, c(1,2)) # tE = Sigma \otimes Sigma
    
    Y[[k]] = Epsilon + ttm(B, X, 3)
  }
  
  return(list(tsr = Y, U = U, X = X, B = B, Sigma = Sigma, gamma = gamma))
}

sparse_sim = function(d, d3, p, R, r3, s, signal, dup) {
  # d: dimension of the network 
  # d3: sample size (dimension of the third mode)

  # p: dimension of the feature
  # R: rank of the intercept
  # r3: rank of the third mode

  # s: sparsity of coefficient tensor (except the intercept)
  # signal: max abs value of the linear predictor entries
  # dup: duplication

  #d = 20; d3 = 20; p = 8; R =4; r3 = 4; s = 0.9; signal = 2; dup = 2

  # mode 3 covariates
  X = cbind(rep(1,d3), randortho(d3)[,1:p])
  
  # intercept matrix
  A = matrix(rnorm(d*R), ncol = R, nrow = d)
  B0 = A%*%t(A)

  G = array(runif(R^2*r3, min = -1, max = 1), dim = c(R,R,r3))
  W1 = randortho(d)[,1:R]
  W2 = randortho(d)[,1:R]
  W3 = randortho(p)[,1:r3]
  
  B = array(0, dim = c(d,d,(p+1)))
  B[,,1] = B0

  B1 = ttl(as.tensor(G), list(W1,W2,W3), c(1,2,3))@data
  
  U1 = ttm(as.tensor(B1), X[,2:(p+1)], 3)@data

  G=G/max(abs(U1))*signal ## rescale subject to entrywise constraint
  U1=U1/max(abs(U1))*signal

  B1 = ttl(as.tensor(G),list(W1,W2,W3),ms = c(1,2,3))@data

  a = 1:(d^2*p); b = round(((1-s)*d^2*p))
  sparse_index = sample(a, b , replace = FALSE)
  B1[sparse_index] = 0

  B[,,2:(p+1)] = B1

  U = ttm(as.tensor(B), X, 3)@data

  tsr = lapply(seq(dup), function(x) array(rbinom(d*d*d3,1,prob = as.vector( 1/(1 + exp(-U)))),dim = c(d,d,d3)))

  return(list(tsr = tsr, X = X, U = U, B = B, G = G))
}
