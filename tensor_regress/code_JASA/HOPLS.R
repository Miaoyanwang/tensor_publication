## This part follows algorithm 1+2 --> 3 in Q. Zhao et al. 2012
## HOPLS: for N-order tensor, N \geq 3
library(rTensor)
## X is a matrix
HOPLS = function(X,tsr, R, Kn, tol){
  # R: number of latent vectors
  # Kn: vector, contains number of loadings in tsr (see fig 2 in Q. Zhao et al. 2012)
  # tol: tolerance, epsilon in algorithm
  E = list(); Fr = list() ## denote residuals of X, tsr
  t = list() ## latent vectors
  Q = list() ## loadings
  D = list() ## tensor of tsr (see fig 2 in Q. Zhao et al. 2012)
  E[[1]] = X; Fr[[1]] = as.tensor(tsr)
  for(r in 1:R){
    if(sum(E[[r]]^2) > tol & sum(Fr[[r]]@data^2) > tol){
      Cr = ttm(Fr[[r]],t(E[[r]]), m = 1)
      tckr = tucker(Cr, ranks = c(1,Kn))
      t_inte = E[[r]]%*%tckr$U[[1]]  ## pr = tckr$U[[1]]  vector
      t[[r]] = t_inte/sqrt(sum(t_inte^2))
      Q[[r]] = tckr$U ; Q[[r]][[1]] = t[[r]]
      Qrt = lapply(1:(1+length(Kn)), function(x) t(Q[[r]][[x]]))  ## get t(Q[[r]])
      #D_inte = ttm(Fr[[r]],t(t[[r]]), m = 1)
      #D[[r]] = ttl(D_inte,Qrt,ms = seq(2,length(Kn)+1))
      D[[r]] = ttl(Fr[[r]],Qrt,ms = seq(1,length(Kn)+1))
      E[[r+1]] = E[[r]] - t[[r]]%*%t(tckr$U[[1]])
      #inte = ttm(D[[r]], t[[r]], m = 1)
      #Fr[[r+1]] = Fr[[r]] - ttl(inte, Q[[r]], m = seq(2,length(Kn)+1))
      Fr[[r+1]] = Fr[[r]] - ttl(D[[r]], Q[[r]], m = seq(1,length(Kn)+1))
    }else break
  }
  pre = 0
  for(r in 1:R){
      inte = ttl(D[[r]], Q[[r]], m = seq(1,length(Kn)+1))
      pre = pre + inte
  }
  return(list(Q = Q, D = D, t = t,pre=pre))
}

### X is a tensor
HOPLS2 = function(X,tsr, R, Ln,Kn, tol){
    # R: number of latent vectors
    # Kn: vector, contains number of loadings in tsr (see fig 2 in Q. Zhao et al. 2012)
    # tol: tolerance, epsilon in algorithm
    E = list(); Fr = list() ## denote residuals of X, tsr
    t = list() ## latent vectors
    G = list() ## loadings
    D = list() ## tensor of tsr (see fig 2 in Q. Zhao et al. 2012)
    
    E[[1]] = as.tensor(X); Fr[[1]] = as.tensor(tsr)
    for(r in 1:R){
        if(sum(E[[r]]@data^2) > tol & sum(Fr[[r]]@data^2) > tol){
            Cr = inner(E[[r]],Fr[[r]])
            tckr = tucker(Cr, ranks = c(Ln,Kn))
            
            pre=preG=preD=list()
            for(i in 1:length(Ln)){
                preG[[i+1]]=pre[[i]]=t(tckr$U[[i]])
            }
            temp=ttl(E[[r]],pre,ms=2:(1+length(Ln)))
            preG[[1]]=preD[[1]]=t[[r]]=rbind(svd(k_unfold(temp,1)@data)$u[,1])
            G[[r]]=ttl(E[[r]],preG,ms=(1:(1+length(Ln))))
              
             for(i in 1:length(Kn)){
                 preD[[i+1]]=t(tckr$U[[(i+length(Ln))]])
             }
             D[[r]]=ttl(Fr[[r]],preD,ms=(1:(1+length(Kn))))
            
            preE=preF=list()
            for(i in 1:(1+length(Ln))){
             preE[[i]]=t(preG[[i]])  
            }
            for(i in 1:(1+length(Kn))){
                preF[[i]]=t(preD[[i]])  
            }
            
            
            E[[r+1]] = E[[r]] - ttl(G[[r]],preE,ms=(1:(1+length(Ln))))
            Fr[[r+1]] = Fr[[r]] - ttl(D[[r]],preF,ms=(1:(1+length(Ln))))
        }else break
    }
    
    return(list(P=preE, Q=preF, G = G, D = D, t = t,pre=pre))
}

inner=function(X,Y){
    dim1=dim(X)
    dim2=dim(Y)
    d=c(dim1[-1],dim2[-1])
    r=dim1[1]
    tensor=array(t(k_unfold(X,1)@data)%*%k_unfold(Y,1)@data,dim=d)

    return(as.tensor(tensor))
}

#result2=ttm(as.tensor(tsr),t(X_covar1),1)
#result=inner(as.tensor(tsr),as.tensor(X_covar1))
#HOPLS = HOPLS2(tsr,tsr,10,c(2,2),c(3,3),0.01)


pre_err = function(HOPLS,tsr,Kn){
  pre = 0
  for(r in 1:R){
    inte = ttl(HOPLS$D[[r]], HOPLS$Q[[r]], m = seq(1,length(Kn)+1))
    pre = pre + inte
  }
  err = sum((pre@data - tsr)^2)
  return(err)
}
#pre_err(HOPLS,tsr,c(2,2))













