######### massive GLM to estimate the covariate effects ###
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

#### valid rank range for core tensor ## the maximum rank should be no smaller than the product of other.
## Input: rank_range = expand.grid(...)
valid_rank=function(rank_range){
    N=dim(rank_range)[1]
    output=NULL
    for(i in 1:N){
        rank=rank_range[i,]
        rank=sort(rank)
        if(rank[1]*rank[2]>=rank[3])
        output=rbind(output,rank_range[i,])
    }
    return(output)
}

### compute loglikelihood under different models
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

#### wrapped simulation function to estimate covariate effects using tensor methods or massive glm methods
conv_rate = function(seed,signal=10,Nsim=20,cons="no",lambda = 1,alpha=10,solver ="GC",c_range,dist="binary",dup=10,d_range,p_range,naive=TRUE){
    #cons can be "non","vanilla","penalty"
    n=m=nrow(d_range)
    s=nrow(c_range)
    error=error_naive=array(0,dim=c(n,s,dup))
    for(i in 1:n){
        for(k in 1:s){
            data = sim_data(seed, d_range[i,], c_range[k,],p_range[i,],dist, dup, signal)
            X_covar1=data$X_covar1
            X_covar2=data$X_covar2
            X_covar3=data$X_covar3
            for(l in 1:dup){
                result=tensor_regress(data$tsr[[l]],X_covar1,X_covar2,X_covar3,c_range[k,],Nsim,cons,lambda,alpha,solver,dist)
                if(naive==TRUE) {
                    naive_C=massive_glm(data$tsr[[l]],X_covar3,dist)
                    error_naive[i,k,l]=mean((naive_C-data$C_ts)^2)
                }
                error[i,k,l]=mean((result$C_ts-data$C_ts)^2) ## use mean not the total F-norm because X has been rescaled. See sim_data(...)
                print(paste(l,"-th replicate---- when dimension is ",d_range[i,1],d_range[i,2],d_range[i,3],"-- covariate is",p_range[i,1],p_range[i,2],p_range[i,3]," --------core is ",c_range[k,1],c_range[k,2],c_range[k,3]))
            }
        }
    }
    return(list(error,error_naive))
}


#### This function is to select the lambda in the penalty constrain version
sele_lambda = function(seed, lambda, ...){
    re = lapply(lambda, FUN = conv_rate, seed = seed, ...)
    re = lapply(seq(length(re)), function(x) re[[x]]$RMSE)
    return(re)
}

