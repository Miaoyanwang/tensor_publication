########## Tensor-SCORE (order 3) ###########

library(rTensor)

######## thresholding
thresholding = function(U, delta){
  
  U = as.matrix(U)
  l2 = apply(U, 1, function(x) sqrt(sum(x^2)))
  index = which(l2 > delta)
  
  U_tilde = U
  U_tilde[index,] = delta * diag(l2[index]^(-1)) %*% U[index,]
  
  return(U_tilde)
}


######### algorithms
hosvd_diag_rm = function(Y, r){
  M1 = unfold(as.tensor(Y), 1, c(3,2))@data
  G = M1 %*% t(M1)
  G = G - diag(G)
  
  U = svd(G)$u[,1:r]
  
  return(U)
}


hooi_reg = function(Y, U0, delta,  max_iter){
  r = dim(U0)[2]
  
  U = U0
  for (iter in 1:max_iter) {
    cat("Reg HOOI iter = ", iter, "\n")
    
    # obtain regularized factor
    U_tilde = thresholding(U, delta)
    P_delta = svd(U_tilde)$u
    
    # update factor
    U_new = svd(unfold(ttl(as.tensor(Y), list(t(P_delta), t(P_delta)), ms = c(2,3)), 1, c(3,2))@data)$u[,1:r]
    
    if(identical(U, U_new)){
      break
    }
    
    U = U_new
  }
  
  return(U)
}


hooi0 = function(Y, U0, max_iter){
  # no regularization
  
  r = dim(U0)[2]
  
  U = U0
  for (iter in 1:max_iter) {
    cat("HOOI iter = ", iter, "\n")
    
    # update factor
    U_new = svd(unfold(ttl(as.tensor(Y), list(t(U), t(U)), ms = c(2,3)), 1, c(3,2))@data)$u[,1:r]
    
    if(identical(U, U_new)){
      break
    }
    
    U = U_new
  }
  
  return(U)
}

####### tensor score

tensor_score_adj = function(Y,r, rm_diag = T, hooi = T,
                            reg_hooi = T,delta = NULL,
                            score_reg = T, t = NULL, max_iter){
  
  if(rm_diag == T){
    # Diagonal removed HOSVD
    U0 = hosvd_diag_rm(Y,r)
  }else if(rm_diag == F){
    # regular HOSVD
    cat("no diagonal remove \n")
    U0 = hosvd(as.tensor(Y), rep(r,3))$U[[1]]
  }
  
  if(hooi == T){
    if(reg_hooi == T){
      # Reg HOOI
      
      if(is.null(delta)){
        li = apply(Y,1, mean)
        delta = 2*sqrt(r)*max(li)/sqrt(sum(li^2))
        cat("use default delta = ", delta,"\n")
      }
      
      Uhat = hooi_reg(Y,U0,delta, max_iter)
    }else if(reg_hooi == F){
      cat("no regularized HOOI \n")
      Uhat = hooi0(Y,U0, max_iter)
    }
  }else if(hooi == F){
    cat("no HOOI \n")
    Uhat = U0
  }
  
  # SCORE 
  R = Uhat/Uhat[,1]
  R = R[,-1]
  
  if(score_reg == T){
    
    if(is.null(t)){
      t = sqrt(log(dim(Y)[1]))
      cat("use default T = ", t,"\n")
    }
    
    R_tilde = thresholding(R, t)
  }else if(score_reg == F){
    cat("no regularized score \n")
    R_tilde = R
  }
  
  
  # kmeans
  z = kmeans(R_tilde, r)$cluster
  
  return(z)
}


# tensor_score = function(Y,r, delta,t,max_iter){
#   
#   # Diagonal removed HOSVD
#   U0 = hosvd_diag_rm(Y,r)
#   
#   # Reg HOOI
#   Uhat = hooi_reg(Y,U0,delta,  max_iter)
#   
#   # SCORE 
#   R = Uhat/Uhat[,1]
#   R = R[,-1]
#   
#   R_tilde = thresholding(R, t)
#   
#   # kmeans
#   z = kmeans(R_tilde, r)$cluster
#   
#   return(z)
# }




