### bricks for new graph

#library(rTensor)
library(pracma)
#library(gtools)
library(MASS)
library(speedglm)
#library(fastglm)
#library(rgl)
#library("RColorBrewer")



#######################

######----------   Tensor regression

## This function is indeed the same as tensor_regress in the package tensorregress
## except the iteration number is forced to be larger than 8
## even though the algorithm has already converged.
## This modification is made for better visualization.

tensor_regress1 = function(tsr,X_covar1 = NULL, X_covar2 = NULL,X_covar3 = NULL, core_shape, Nsim=20, cons = c("non","vanilla","penalty"), lambda = 0.1, alpha = 1, solver ="CG",dist = c("binary", "poisson","normal")){
  
  
  tsr = as.tensor(tsr)
  Y_1 = unfold(tsr, row_idx = 1, col_idx = c(2,3))@data
  Y_2 = unfold(tsr, row_idx = 2, col_idx = c(1,3))@data
  Y_3 = unfold(tsr, row_idx = 3, col_idx = c(1,2))@data
  d1 = dim(tsr)[1] ; d2 = dim(tsr)[2] ; d3 = dim(tsr)[3]
  r1 = core_shape[1] ; r2 = core_shape[2] ; r3 = core_shape[3]
  
  ###  check whether unsupervised on each mode
  un_m1 = FALSE ; un_m2 = FALSE ; un_m3 = FALSE
  if(is.null(X_covar1)|(identical(X_covar1,diag(d1)))) {X_covar1 = diag(d1) ; un_m1 = TRUE}
  if(is.null(X_covar2)|(identical(X_covar2,diag(d2)))) {X_covar2 = diag(d2) ; un_m2 = TRUE}
  if(is.null(X_covar3)|(identical(X_covar3,diag(d3)))) {X_covar3 = diag(d3) ; un_m3 = TRUE}
  p1 = dim(X_covar1)[2] ; p2 = dim(X_covar2)[2] ; p3 = dim(X_covar3)[2]
  
  
  if(dist=="binary"){
    tsr.transform=as.tensor(2*tsr@data-1)
  }else if(dist=="poisson"){
    tsr.transform=as.tensor(log(tsr@data+0.1))###?? new initilization
  }else if (dist=="normal"){
    tsr.transform=tsr
  }
  
  C_ts=ttl(tsr.transform,list(ginv(X_covar1),ginv(X_covar2),ginv(X_covar3)),ms=c(1,2,3))
  
  tckr = tucker(C_ts, ranks = core_shape)
  #W1 = tckr$U[[1]] ; W2 = tckr$U[[2]] ; W3 = tckr$U[[3]] ## initilized tucker factors
  #G = tckr$Z
  W1=randortho(p1)[,1:core_shape[1]];W2=randortho(p2)[,1:core_shape[2]];W3=randortho(p3)[,1:core_shape[3]]
  G=ttl(C_ts,list(t(W1),t(W2),t(W3)),ms=1:3)
  
  A = X_covar1%*%W1
  B = X_covar2%*%W2
  C = X_covar3%*%W3
  
  core=update_core(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist)
  G=core$G
  lglk=core$lglk
  violate=core$violate
  iter = 1
  for(n in 1:Nsim){
    
    
    ## parameter from previous step
    
    W10 = W1 ; W20 = W2 ; W30 = W3 ; G0=G; A0=A;B0=B;C0=C;lglk0=tail(lglk,1);
    
    ###### update W1
    G_BC = ttl(G, list(B,C), ms = c(2,3))
    G_BC1 = unfold(G_BC, row_idx = 1, col_idx = c(2,3))@data
    
    if(un_m1) {re = glm_mat(t(Y_1),t(G_BC1),dist=dist) ## no covariate
    } else {re = glm_two(Y = Y_1, X1 = X_covar1, X2 = G_BC1, dist=dist)}
    
    if(dim(re[[1]])[1]==1) W1=t(re[[1]]) else W1 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W1*
    qr_res=qr(W1)
    W1=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),1)
    #print("W1 Done------------------")
    
    ##### calculate A
    A = X_covar1%*%W1;
    ##### update W2
    G_AC = ttl(G, list(A,C), ms = c(1,3))
    G_AC2 = unfold(G_AC, row_idx = 2, col_idx = c(1,3))@data
    
    if(un_m2) {re = glm_mat(t(Y_2),t(G_AC2),dist=dist)
    } else {re = glm_two(Y_2, X_covar2, G_AC2, dist=dist)}
    
    if(dim(re[[1]])[1]==1) W2=t(re[[1]]) else W2 = as.matrix(re[[1]])
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W2*
    qr_res=qr(W2)
    W2=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),2)
    #print("W2 Done------------------")
    
    ##### calculate B
    B = X_covar2%*%W2;
    
    
    ###### update W3
    G_AB = ttl(G, list(A,B), ms = c(1,2))
    G_AB3 = unfold(G_AB, row_idx = 3, col_idx = c(1,2))@data
    
    if(un_m3) {re = glm_mat(t(Y_3),t(G_AB3),dist=dist)
    } else {re = glm_two(Y_3, X_covar3, G_AB3,dist=dist)}
    
    if(dim(re[[1]])[1]==1) W3=t(re[[1]]) else W3 = as.matrix(re[[1]])
    
    lglk = c(lglk,re[[2]])
    
    ## orthogonal W3*
    qr_res=qr(W3)
    W3=qr.Q(qr_res)
    G=ttm(G,qr.R(qr_res),3)
    #print("W3 Done------------------")
    
    ##### calculate C
    C = X_covar3%*%W3;
    
    #########-----------------------------------------------
    ###  obtain core tensor under constraint
    core=update_core(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist)
    G=core$G
    lglk=c(lglk,core$lglk)
    violate=c(violate,core$violate)
    
    #print("G Done------------------")
    
    message(paste(n,"-th  iteration -- when dimension is",d1,d2,d3,"- rank is ",r1,r2,r3," -----------------"))
    #print(paste(n,"-th  iteration"))
    
    if(iter < 8){
      iter = iter+1
      next
    }
    
    
    if ((tail(lglk,1)-lglk0)/abs(lglk0)<= 0.0001 & tail(lglk,1)>= lglk0 ){
      message(paste(n,"-th iteration: convergence"))
      break
    } else if (tail(lglk,1)-lglk0 < 0) {
      W1 = W10 ; W2 = W20 ; W3 = W30; G=G0; lglk=lglk[-c((length(lglk)-3):length(lglk))];
      A=A0;B=B0;C=C0;
      break
    }
    
  }
  
  U=ttl(G,list(A,B,C),ms = c(1,2,3))@data
  
  sigma_est=mean((tsr@data-U_to_mean(U,dist))^2)
  
  return(list(W = list(W1 = W1,W2 = W2,W3 = W3),G = G@data,U=U, C_ts=ttl(G,list(W1,W2,W3),ms = c(1,2,3))@data,lglk = lglk, sigma=sigma_est,violate = violate))
}




######----------   functions for update

## This part used for penalty (Conjugate gradient method)

loss = function(beta,y,X,lambda,alpha,dist){
  U = X %*% beta
  L= - loglike(y,U,dist)
  if(max(abs(U))>alpha) return(Inf)
  else{
    L2 =L - lambda * sum(log(1 - (U/ alpha)^2))  ## object function
    return(c(L2))
  }
}

loglike=function(data,linearpre,dist){
  if(dist=="binary"){
    p=plogis(linearpre)## log-likelihood
    L=dbinom(data,1,p,log=TRUE)
    return(sum(L))
  }
  else if (dist=="normal"){
    sigma_est=mean((data-linearpre)^2)
    L=dnorm(data,linearpre,sqrt(sigma_est),log=TRUE)
    sum(L)##log-likelihood
  }
  else if (dist=="poisson"){
    lambda=exp(linearpre)
    L=dpois(data,lambda, log = TRUE) ## log-likelihood for Poisson
    return(sum(L))
  }
}

loss_gr <- function(beta,y,X,lambda,alpha,dist){
  U = X %*% beta
  if(dist=="binary") p=plogis(U)
  else if (dist=="normal") p=U
  else if (dist=="poisson") p=exp(U)
  L_g=t(X) %*% (p - y)
  penal_g = 2 * lambda * t(X) %*% (U/(alpha^2 - U^2))
  return(c(L_g) + c(penal_g))
}


glm_modify=function(y,x,dist){
  if(dist=="binary"){
    fit1=suppressWarnings(speedglm(y~-1+x,family=binomial(link="logit")))
    return(list(coef(fit1), logLik(fit1)))
  }
  else if (dist=="normal"){
    fit1 =speedlm(y~-1+x)
    
    return(list(coef(fit1), logLik(fit1)))
  }
  else if (dist=="poisson"){
    fit1=suppressWarnings(speedglm(y~-1+x,family=poisson(link="log")))
    return(list(coef(fit1), logLik(fit1)))
  }
}

############################################################

#    form U = X*Beta ## in parallel
glm_mat = function(Y,X,dist){
  R = dim(X)[2]   ## R in note
  p = dim(Y)[2]
  ma = mapply(glm_modify, y =  as.data.frame(Y),MoreArgs = list(X),dist=dist)
  re = t(matrix(unlist(ma), nrow = p, byrow=T))
  beta = re[1:R,]
  
  lglk = sum(re[R + 1,])
  return(list(t(beta),lglk))
}

###########---------  GLM on two modes
glm_two = function(Y, X1, X2, dist){
  ## Y_size = m * n
  # logit(E(Y)) = X1 %*% coe %*% X2
  m = dim(Y)[1] ; n = dim(Y)[2]
  q1 = dim(X1)[2] ; q2 = dim(X2)[1]
  
  N_long = kronecker_list(list(t(X2),X1))
  
  mod_re=glm_modify(as.vector(Y),N_long,dist)
  coe = mod_re[[1]]
  coe = matrix(coe, nrow = q1, ncol = q2)
  lglk= mod_re[[2]]
  return(list(coe=coe,lglk=lglk))
}



update_core=function(tsr,G,A,B,C,core_shape,cons,lambda,alpha,solver,dist){
  
  M_long = kronecker_list(list(C,B,A))
  U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
  violate=0
  
  if (cons == 'penalty'){
    if(max(abs(U))>=alpha){
      G=G/max(abs(U))*(alpha-0.01)
      U = U/max(abs(U))*(alpha-0.01)
      message("Violate constrain ------------------")
      violate = 1
    }
    mod_re = optim(par = as.vector(G@data),loss,loss_gr,y = as.vector(tsr@data),X =M_long, lambda = lambda, alpha = alpha, method = solver,dist=dist)
    coe = mod_re$par
    G = as.tensor(array(data = coe,dim = c(core_shape)))
    lglk = -mod_re$value
  }else {
    mod_re = glm_modify(as.vector(tsr@data), M_long,dist=dist)
    coe = mod_re[[1]]
    G = as.tensor(array(data = coe,dim = core_shape))
    U = ttl(G,list(A,B,C),ms = c(1,2,3))@data
    lglk = loglike(tsr@data,U,dist)
    
    if(cons== 'vanilla' & max(abs(U))>=alpha){
      G=G/max(abs(U))*(alpha-0.01)
      U = U/max(abs(U))*(alpha-0.01)
      message("Violate constrain ------------------")
      violate = 1
    }
    lglk = loglike(tsr@data,U,dist)}
  
  
  return(list(G=G,violate=violate,lglk=lglk))
}



U_to_mean=function(U,dist){
  if(dist=="normal") return(U)
  else if (dist=="binary") return(plogis(U))
  else if (dist=="poisson") return(exp(U))
}
###
###----  This function shows how we select rank
##   recommend to use non-constrain verison to select rank



massive_glm=function(response,X,dist){
  d1=dim(response)[1]
  d2=dim(response)[2]
  coe=array(0,dim=c(d1,d2,dim(X)[2]))
  if(dist=="normal"){
    for(i in 1:d1){
      for(j in 1:d2){
        fit=lm(response[i,j,]~-1+X)
        coe[i,j,]=coef(fit)
      }
    }
  }
  else if(dist=="binary"){
    for(i in 1:d1){
      for(j in 1:d2){
        fit=suppressWarnings(glm(response[i,j,]~-1+X,family=binomial("logit")))
        coe[i,j,]=coef(fit)
      }
    }
  }
  else if(dist=="poisson"){
    for(i in 1:d1){
      for(j in 1:d2){
        fit=suppressWarnings(glm(response[i,j,]~-1+X,family=poisson("log")))
        coe[i,j,]=coef(fit)
      }
    }
  }
  return(coe)
}
