#############################################
########## this is the code #################
########## for binary undirected  ###########
#############################################


library(rTensor)
#library(MASS)
#library(gplots)
#library(glmnet)
#library(igraph)
#library(logistf)
#library(mgcv)
#library(RMTstat)
#library(reliaR)

EPS = 1e-3

logit<-function(x){
	log(x/(1-x))
}


logistic<-function(x){
	1-1/(1+exp(x))
}


norm_vec <- function(x) sqrt(sum(x^2))
scale_vec <- function(x) x/sqrt(sum(x^2))



mytruncate<-function(x,s){
	xtruncate = rep(0, length(x))
	xtruncate[ which(rank(-abs(x))<= s)] = x[which(rank(-abs(x))<= s)]
	xtruncate
}

#######################################################################
## This function estimates the low rank parameter only ################
## gamma_ar is the step size and it needs tuning ######################
## if gamma_ar is too large, then the algorithm will not converge #####
## gamma_ar may set to decrease with iter number for faster convergence  
#######################################################################

mySymGLM0<-function(Z, X, R=null, sparsity=null, niter=50){
	n<-dim(Z)[1]
	N<-dim(Z)[3]
	p<-dim(X)[2]
    ### initialize ####		
    Abar<-apply(Z@data,c(1,2),mean)
    Abar[Abar<=0]<-EPS
    Abar[Abar>=1]<-1-EPS	
	c<-rep(1,N)
    initial<-svd(logit(Abar))
    w<-initial$d[1:R]*(sign(initial$u[1,1:R])*sign(initial$v[1,1:R]))
    A<-initial$u[,1:R] 

    iter<-1; cond<-FALSE; 
    loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum(log(1+exp((A%*%diag(w)%*%t(A)%o%c))))

    while(!cond){
    	loglikelihood_old<-loglikelihood
    	A_tmp<-A
    	w_tmp<-w
    	vec_a<-ttm(ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c),t(A),m = 2),t(c),m=3)[,,1]@data
    	direction<-sweep(vec_a,2,w,"*")
    	gamma_ar<-1/apply(abs(A_tmp)%o%c,2,sum)/abs(w_tmp)/2
    	A_tmp<-A_tmp+sweep(direction,2,gamma_ar,"*")
    	w<-apply(A_tmp,2,norm_vec)*apply(A_tmp,2,norm_vec)*w_tmp
    	A<-apply(A_tmp,2,scale_vec)
    	Theta<-(A%*%diag(w)%*%t(A))
    	iter<-iter+1
    	loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum(log(1+exp((A%*%diag(w)%*%t(A)%o%c))))
    	cond<-(iter > niter)| (loglikelihood-loglikelihood_old <= EPS)
    }
    BIC<--loglikelihood+log(N*(n+1)*n/2)*(R*n)+2*(R*n)*log(n*n/2)
    return(list(loglikelihood=loglikelihood,BIC=BIC,w=w,A=A))   
}

mySymGLM0_normal<-function(Z, X, R=null, sparsity=null, niter=50){
  n<-dim(Z)[1]
  N<-dim(Z)[3]
  p<-dim(X)[2]
  ### initialize ####		
  Abar<-apply(Z@data,c(1,2),mean)
  # Abar[Abar<=0]<-EPS
  # Abar[Abar>=1]<-1-EPS	
  c<-rep(1,N)
  #initial<-svd(logit(Abar))
  initial<- svd(Abar)
  #w<-initial$d[1:R]*(sign(initial$u[1,1:R])*sign(initial$v[1,1:R]))
  w<-initial$d[1:R]
  A<-initial$u[,1:R] 
  
  iter<-1; cond<-FALSE; 
  # loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum(log(1+exp((A%*%diag(w)%*%t(A)%o%c))))
  loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum((A%*%diag(w)%*%t(A)%o%c)^2)/2
  
  while(!cond){
    loglikelihood_old<-loglikelihood
    A_tmp<-A
    w_tmp<-w
    #vec_a<-ttm(ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c),t(A),m = 2),t(c),m=3)[,,1]@data
    vec_a<-ttm(ttm(Z-((A%*%diag(w)%*%t(A))%o%c),t(A),m = 2),t(c),m=3)[,,1]@data
    direction<-sweep(vec_a,2,w,"*")
    gamma_ar<-1/apply(abs(A_tmp)%o%c,2,sum)/abs(w_tmp)/2
    A_tmp<-A_tmp+sweep(direction,2,gamma_ar,"*")
    w<-apply(A_tmp,2,norm_vec)*apply(A_tmp,2,norm_vec)*w_tmp
    A<-apply(A_tmp,2,scale_vec)
    Theta<-(A%*%diag(w)%*%t(A))
    iter<-iter+1
    #loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum(log(1+exp((A%*%diag(w)%*%t(A)%o%c))))
    loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum((A%*%diag(w)%*%t(A)%o%c)^2)/2
    cond<-(iter > niter)| (loglikelihood-loglikelihood_old <= EPS)
  }
  BIC<--loglikelihood+log(N*(n+1)*n/2)*(R*n)+2*(R*n)*log(n*n/2)
  return(list(loglikelihood=loglikelihood,BIC=BIC,w=w,A=A))   
}

##################################################################################
## This estimates both the low rank and the sparse parameters ####################
## gamma_ar and gamma_BB are the step sizes and both require some tuning #########
## they may decrease with iter number for faster convergence  ####################
## use eBIC when dimension is larger in real applications ########################
##################################################################################

mySymGLM<-function(Z, X,  R=null, sparsity=null, niter=50){
  
    n<-dim(Z)[1]
    N<-dim(Z)[3]
    p<-dim(X)[2]

    ### initialize ####     
    Abar<-apply(Z@data,c(1,2),mean)
    Abar[Abar<=0]<-EPS
    Abar[Abar>=1]<-1-EPS    
    c<-rep(1,N)
    initial<-svd(logit(Abar))
    w<-initial$d[1:R]*(sign(initial$u[1,1:R])*sign(initial$v[1,1:R]))
    A<-initial$u[,1:R] 
    BB_truncate<-new("Tensor",3L,c(n,n,p),data=0)

    iter<-1; cond<-FALSE; 
    loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(log(1+exp(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data))) 
    
    while(!cond){
        loglikelihood_old<-loglikelihood
        A_tmp<-A
        w_tmp<-w
        vec_a<-ttm(ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3)@data),t(A),m = 2),t(c),m=3)[,,1]@data*2
        direction<-sweep(vec_a,2,w,"*")
        gamma_ar<-1/apply(abs(A_tmp)%o%c,2,sum)/abs(w_tmp)/4
        A_tmp<-A_tmp+sweep(direction,2,gamma_ar,"*")
        w<-apply(A_tmp,2,norm_vec)*apply(A_tmp,2,norm_vec)*w_tmp
        A<-apply(A_tmp,2,scale_vec)
        Theta<-(A%*%diag(w)%*%t(A))
        gamma_BB<-1/20/(iter+2)
        tsnr_b<-ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate, X, m = 3)@data),t(X),m=3)@data*2
        BB<-BB_truncate + gamma_BB * tsnr_b
        BB_truncate<-new("Tensor",3L,c(n,n,p),data=mytruncate(vec(BB),sparsity))
        iter<-iter+1
        loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(log(1+exp(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data))) 
        cond<-(iter > niter) | (loglikelihood-loglikelihood_old <= EPS) 
    }
    BIC<--loglikelihood+log(N*(n+1)*n/2)*(R*n+sparsity/2)#+2*(R*n+sparsity/2)*log(n*n/2+n*n*p/2)
    return(list(loglikelihood=loglikelihood,BIC=BIC,w=w,A=A,BB=BB_truncate))  
}


mySymGLM_normal<-function(Z, X,  R=null, sparsity=null, niter=50){ #normal version
  
  n<-dim(Z)[1]
  N<-dim(Z)[3]
  p<-dim(X)[2]
  
  ### initialize ####     
  Abar<-apply(Z@data,c(1,2),mean)
  #Abar[Abar<=0]<-EPS
  #Abar[Abar>=1]<-1-EPS    
  c<-rep(1,N)
  #initial<-svd(logit(Abar))
  initial<- svd(Abar)
  #w<-initial$d[1:R]*(sign(initial$u[1,1:R])*sign(initial$v[1,1:R]))
  w<-initial$d[1:R]
  A<-initial$u[,1:R] 
  BB_truncate<-new("Tensor",3L,c(n,n,p),data=0)
  
  iter<-1; cond<-FALSE; 
  #loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(log(1+exp(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data))) 
  loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data^2)/2 
  while(!cond){
    loglikelihood_old<-loglikelihood
    A_tmp<-A
    w_tmp<-w
    #vec_a<-ttm(ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3)@data),t(A),m = 2),t(c),m=3)[,,1]@data*2
    vec_a<-ttm(ttm(Z-((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3)@data),t(A),m = 2),t(c),m=3)[,,1]@data*2
    direction<-sweep(vec_a,2,w,"*")
    gamma_ar<-1/apply(abs(A_tmp)%o%c,2,sum)/abs(w_tmp)/4
    A_tmp<-A_tmp+sweep(direction,2,gamma_ar,"*")
    w<-apply(A_tmp,2,norm_vec)*apply(A_tmp,2,norm_vec)*w_tmp
    A<-apply(A_tmp,2,scale_vec)
    Theta<-(A%*%diag(w)%*%t(A))
    gamma_BB<-1/20/(iter+2)
    #tsnr_b<-ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate, X, m = 3)@data),t(X),m=3)@data*2
    tsnr_b<-ttm(Z-((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate, X, m = 3)@data),t(X),m=3)@data*2
    BB<-BB_truncate + gamma_BB * tsnr_b
    BB_truncate<-new("Tensor",3L,c(n,n,p),data=mytruncate(vec(BB),sparsity))
    iter<-iter+1
    #loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(log(1+exp(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data))) 
    loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data^2)/2 
    cond<-(iter > niter) | (loglikelihood-loglikelihood_old <= EPS) 
  }
  BIC<--loglikelihood+log(N*(n+1)*n/2)*(R*n+sparsity/2)#+2*(R*n+sparsity/2)*log(n*n/2+n*n*p/2)
  return(list(loglikelihood=loglikelihood,BIC=BIC,w=w,A=A,BB=BB_truncate))  
}

########################################################
###################### tuning  #########################
####### range for R and s can be adjusted   ############
########################################################

best_SymGLM<-function(Z,X){
    n<-dim(Z)[1]
    p<-dim(X)[2]
    result<-list()
    BIC<-c()
    for(R in 2:15){
        output<-mySymGLM0(Z=Z,X=X,R=R,sparsity=0,niter=100)
        BIC[R]<-output$BIC
        plot(BIC,main="BIC for rank selection")
    }
    R<-which.min(BIC)
    BIC2<-c()
    for(i in c(1:21)){
        s<-c(10^seq(-2,0,0.1))[i]
        output<-mySymGLM(Z=Z,X=X,R=R,sparsity=n*n*p*s,niter=100)
        result[[i]]<-output
        BIC2[i]<-output$BIC
        plot(BIC2,main="BIC for sparsity selection")
    }
    
    # original line is
    # result[[which.min(BIC)]]

    s = c(10^seq(-2,0,0.1))[which.min(BIC2)]
    # result[[which.min(BIC2)]] # typo here
    return(list(result = result[[which.min(BIC2)]], best_R = R, best_s = s))
}

best_SymGLM_normal<-function(Z,X){
  n<-dim(Z)[1]
  p<-dim(X)[2]
  result<-list()
  BIC<-c()
  for(R in 2:15){
    output<-mySymGLM0_normal(Z=Z,X=X,R=R,sparsity=0,niter=100)
    BIC[R]<-output$BIC
    plot(BIC,main="BIC for rank selection")
  }
  R<-which.min(BIC)
  BIC2<-c()
  for(i in c(1:21)){
    s<-c(10^seq(-2,0,0.1))[i]
    output<-mySymGLM_normal(Z=Z,X=X,R=R,sparsity=n*n*p*s,niter=100)
    result[[i]]<-output
    BIC2[i]<-output$BIC
    plot(BIC2,main="BIC for sparsity selection")
  }
  
  # original line is
  # result[[which.min(BIC)]]
  
  s = c(10^seq(-2,0,0.1))[which.min(BIC2)]
  # result[[which.min(BIC2)]] # typo here
  return(list(result = result[[which.min(BIC2)]], best_R = R, best_s = s))
}




####################################################
################# Simulation #######################
####################################################

generateTensor<-function(N,n,X,Theta,BB,symmetric=TRUE){
    l<-length(vec(Theta))
    Z<-new("Tensor",3L,c(n,n,N),data=rbinom(n=l,size=1,prob=logistic(vec(Theta+ttm(BB,X,m = 3)))))
    if(symmetric==TRUE){
        for(i in 1:N){
            slice<-matrix(0,n,n)
            slice[upper.tri(slice)]<-Z[,,i]@data[upper.tri(Z[,,i]@data)]
            slice<-slice+t(slice) 
            diag(slice)<-diag(Z[,,i]@data)      
            Z[,,i]@data<-slice
        }
        }else{
            Z<-new("Tensor",3L,c(n,n,N),data=rbinom(n=l,lambda=exp(vec(Theta+ttm(BB,X,m = 3)))))
        }
    Z
}


generateTensor_normal<-function(N,n,X,Theta,BB,symmetric=TRUE){
  l<-length(vec(Theta))
  #Z<-new("Tensor",3L,c(n,n,N),data=rbinom(n=l,size=1,prob=logistic(vec(Theta+ttm(BB,X,m = 3)))))
  Z<-new("Tensor",3L,c(n,n,N),data=rnorm(n=l,vec(Theta+ttm(BB,X,m = 3))))
  if(symmetric==TRUE){
    for(i in 1:N){
      slice<-matrix(0,n,n)
      slice[upper.tri(slice)]<-Z[,,i]@data[upper.tri(Z[,,i]@data)]
      slice<-slice+t(slice) 
      diag(slice)<-diag(Z[,,i]@data)      
      Z[,,i]@data<-slice
    }
  }else{
    #Z<-new("Tensor",3L,c(n,n,N),data=rbinom(n=l,lambda=exp(vec(Theta+ttm(BB,X,m = 3)))))
    Z<-new("Tensor",3L,c(n,n,N),data=rnorm(n=l,vec(Theta+ttm(BB,X,m = 3))))
  }
  Z
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

