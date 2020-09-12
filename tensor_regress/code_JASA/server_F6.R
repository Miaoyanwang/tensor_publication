source("simulation.R")
source("comparison.R")
library(tensorregress)

set.seed(0)
seed=0
## combinations of signal and rank
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))
d=c(20,40,60,80,100)

##
dup=30;
final=array(0,dim=c(2,2,5,4))
finalsd=array(0,dim=c(2,2,5,4))
err_res=err_holrr=err_pls=err_tpg=array(0,dim=c(2,2,5,dup))
dist="normal";


for(s in 1:2){ ## signal level
    for(r in 1:2){ ## rank
        for(i in 1:5){ ## dimension with informative mode
            whole_shape = c(d[i],20,20)
            
            for(n in 1:dup){ ## simulation replicates
                core_shape=core_range[r,]
                signal=signal_range[s]
                
                data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d[i],0,0),dist=dist, dup=dup, signal=signal)
                res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
                
                err_res[s,r,i,n]= mean((res$U - data$U)^2)
                
                ### holrr
                holrr = HOLRR(X = data$X_covar1, data$tsr[[n]], core_shape = core_shape) ## same rank as ours
                err_holrr[s,r,i,n] = mean((holrr$pre@data -data$U)^2)
                
                ## pls ## choose hyperparameter using the best one
                test=NULL
                for(R in 1:8){
                    pls = HOPLS(X=data$X_covar1,data$tsr[[n]],R=R,Kn=c(4,6),tol=0.01) ## rank
                    test=c(test,mean((pls$pre@data - data$U)^2))
                }
                err_pls[s,r,i,n]=min(test)
                
                ## Yu ## choose hyperparameter using the best one
                test=NULL
                for(R in 1:8){
                    tpg=TPG(X=data$X_covar1,data$tsr[[n]],R=R,eta=.1,iter=15,tol=10-6)
                    test=c(test,mean((tpg$U - data$U)^2))
                }
                err_tpg[s,r,i,n]=min(test)
            }
            final[s,r,i,1]=mean(err_res[s,r,i,])
            final[s,r,i,2]=mean(err_holrr[s,r,i,])
            final[s,r,i,3]=mean(err_pls[s,r,i,])
            final[s,r,i,4]=mean(err_tpg[s,r,i,])
            
            
            finalsd[s,r,i,1]=sd(err_res[s,r,i,])
            finalsd[s,r,i,2]=sd(err_holrr[s,r,i,])
            finalsd[s,r,i,3]=sd(err_pls[s,r,i,])
            finalsd[s,r,i,4]=sd(err_tpg[s,r,i,])
        }
    }
}

save(final,finalsd,file="Figure6.RData")
