

##------------------------   TPG + ITP
library(rTensor)
## This part follows algorithm 2 in Rose Yu et al. 2016
## ITP: for N-order tensor, N = 3
ITP = function(W, R, U){
    # R: scalar constraining tucker rank
    # W: Coefficient tensor: d1, d2, T
    # U: list containing factor matrices
    W = as.tensor(W)
    U1 = U[[1]]; U2 = U[[2]]; U3 = U[[3]]
    i = 1
    while(i <= R){
        u10 = U1[,i]; u20 = U2[,i]; u30 = U3[,i]
        for(j in 1:100){  ##  power iteration
            u11 = ttl(W,list(t(as.matrix(u20)), t(as.matrix(u30))), ms = c(2,3))@data
            u11 = as.vector(u11); u11 = u11/sum(u11^2)
            u21 = ttl(W,list(t(as.matrix(u11)), t(as.matrix(u30))), ms = c(1,3))@data
            u21 = as.vector(u21); u21 = u21/sum(u21^2)
            u31 = ttl(W,list(t(as.matrix(u11)), t(as.matrix(u21))), ms = c(1,2))@data
            u31 = as.vector(u31); u31 = u31/sum(u31^2)
            u10 = u11; u20 = u21; u30 = u31
            if((sum((u11 - u10)^2) + sum((u21 - u20)^2) + sum((u31 - u30)^2)) <= 0.05) break
        }
        U1[,i] = u11;  U2[,i] = u21;  U3[,i] = u31 
        i = i + 1
    }
    U = list(U1,U2,U3)
    UUT = lapply(1:3, function(x) U[[x]]%*%t(U[[x]]))
    W = ttl(W, UUT, m = c(1,2,3))@data
    return(list(W = W, U = U))
}

sketch = function(x){
    i = sample(1:length(x),1)
    x[i] = sample(c(1,-1),1,prob = c(0.5,0.5))
    return(x)
}

## This part follows algorithm 1 in Rose Yu et al. 2016
## TPG: for N-order tensor, N = 3
TPG = function(X, tsr, R, eta, iter, tol){
    # tsr: d3, d2, T
    # X: d3, d1
    # W: Coefficient tensor: d1, d2, T
    # eta: step size
    # iter: iteration times
    # tol: tolerance
    X_input=X
    d3 = dim(X)[1]; Tt = dim(tsr)[3]; d2 = dim(tsr)[2]; d1 = dim(X)[2]
    #l = sample(seq(d3-1))  # sketch to lower dim
    #S = matrix(0,l, d3)
    #S = as.matrix(sapply(as.data.frame(S), sketch))
    l=d3; S=diag(d3) ## no sketch
    tsr = as.tensor(tsr)
    tsr = ttm(tsr,S,m = 1)
    X = S %*% X
    W = array(0, c(d1,d2,Tt));  W = as.tensor(W)
    gr = ttm(tsr - ttm(W,X,m = 1),t(X),m=1)
    W = W - eta*gr
    U = tucker(W, ranks = rep(R,3))$U
    W = W@data
    loss = 0
    for(i in 1:iter) {
        loss0 = loss
        gr = -ttm(tsr - ttm(as.tensor(W),X,m = 1),t(X),m=1)@data
        W = W - eta*gr
        ITPres = ITP(W, R, U)
        W = ITPres$W
        U = ITPres$U
        loss = sum((tsr@data - ttm(as.tensor(W),X,m = 1)@data)^2)
        cat("step", i, ":", loss, "\n")
        #print(W)
        if(loss <= tol | sum((loss - loss0)^2)<= 0.01 ){
            break
        }
    }
    U=ttm(as.tensor(W),X_input,m = 1)@data
    return(list(W=W,U=U))
}
