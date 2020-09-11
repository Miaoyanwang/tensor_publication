install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
library(tensorregress)
source("HOLRR.R")
source("HOPLS.R")
source("Yu.R")
library(ggplot2)


############# covariate ######
signal=2.5 ## low noise (sd = 1/2.5)
dup=10;
core_shape=c(3,3,3)
drange=1;
d=40
final=array(0,dim=c(3,4))
finalsd=array(0,dim=c(3,4))
whole_shape = rep(d,3)
err_res=err_holrr=err_pls=err_tpg=matrix(0,nrow=3,ncol=dup)
dist="normal"
#for(s in 1:2){
#for(i in 1:drange){


for(i in 1:3){
    for(n in 1:dup){
        
        if(i==1){ ## one mode covaraite
            data=sim_data_new(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0,0),dist=dist, dup=1, signal=2.5)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }else if(i==2){ ## two mode covaraite
            data=sim_data_new(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0.4*d,0),dist=dist, dup=1, signal=2.5)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,X_covar2 = data$X_covar2,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }else if (i==3){## three model covariate
            data=sim_data_new(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0.4*d,0.4*d),dist=dist, dup=1, signal=2.5)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }
        err_res[i,n]= mean((res$U - data$U)^2)
        
        ### holrr
        holrr = HOLRR(X = data$X_covar1, data$tsr[[1]], core_shape = core_shape)
        err_holrr[i,n] = mean((holrr$pre@data -data$U)^2)
        
        
        ## pls
        pls = HOPLS(X=data$X_covar1,data$tsr[[1]],2,c(4,6),0.01) ## rank
        err_pls[i,n]=mean((pls$pre@data - data$U)^2)
        
        ## Yu
        test=NULL
        for(R in 1:8){
            tpg=TPG(X=data$X_covar1,data$tsr[[1]],R=R,eta=.1,iter=15,tol=10-6)
            test=c(test,mean((tpg$U - data$U)^2))
        }
        err_tpg[i,n]=min(test)
    }
    final[i,1]=mean(err_res[i,])
    final[i,2]=mean(err_holrr[i,])
    final[i,3]=mean(err_pls[i,])
    final[i,4]=mean(err_tpg[i,])
    
    
    finalsd[i,1]=sd(err_res[i,])
    finalsd[i,2]=sd(err_holrr[i,])
    finalsd[i,3]=sd(err_pls[i,])
    finalsd[i,4]=sd(err_tpg[i,])
    
}
#}
