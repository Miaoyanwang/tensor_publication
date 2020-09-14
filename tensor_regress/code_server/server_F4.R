source("simulation.R")
library(tensorregress)

set.seed(0)
seed=0
## one-side covariates
c_range=c(5,5,5) ## selected rank
c_range=rbind(c_range,c(10,10,5),c(15,15,5))
optimal_rank=c(5,5,5)
optimal_rank=rbind(optimal_rank,c(10,10,5),c(14,14,5)) ## to be determined
dup=30

### simulate data from normal model
table=table_naive=matrix(0,nrow=3,ncol=dup)
##############################
for(s in 1:3){
    data=sim_data(whole_shape=c(20,20,50),core_shape=c_range[s,],p=c(0,0,5),dist="normal",dup=dup,signal=5,block=c(TRUE,TRUE,FALSE))
    for(i in 1:dup){
        result=tensor_regress(tsr=data$tsr[[i]],X_covar1 =data$X_covar1, X_covar2 = data$X_covar2,X_covar3= data$X_covar3, core_shape=optimal_rank[s,], Nsim=10, cons="no", lambda = 0.1, alpha = 1, solver ="CG",dist="normal")
        table[s,i]=mean((result$C_ts-data$C_ts)^2) ## tensor method
        naive_C=massive_glm(data$tsr[[i]],data$X_covar3,dist="normal")
        table_naive[s,i]=mean((naive_C-data$C_ts)^2) ## GLM method
    }
}
save(table,table_naive,file="comparison_normal.RData")

### simulate data from poisson model
table=table_naive=matrix(0,nrow=3,ncol=dup)
##############################
for(s in 1:3){
    data=sim_data(whole_shape=c(20,20,50),core_shape=c_range[s,],p=c(0,0,5),dist="poisson",dup=dup,signal=5,block=c(TRUE,TRUE,FALSE))
    for(i in 1:dup){
        result=tensor_regress(tsr=data$tsr[[i]],X_covar1 =data$X_covar1, X_covar2 = data$X_covar2,X_covar3= data$X_covar3, core_shape=optimal_rank[s,], Nsim=10, cons="no", lambda = 0.1, alpha = 1, solver ="CG",dist="poisson") ## tensor method
        table[s,i]=mean((result$C_ts-data$C_ts)^2)
        naive_C=massive_glm(data$tsr[[i]],data$X_covar3,dist="poisson") ## glm method
        table_naive[s,i]=mean((naive_C-data$C_ts)^2)
    }
}
save(table,table_naive,file="comparison_poisson.RData")


### simulate data from binary model
table=table_naive=matrix(0,nrow=3,ncol=dup)
##############################
for(s in 1:3){
    data=sim_data(whole_shape=c(20,20,50),core_shape=c_range[s,],p=c(0,0,5),dist="binary",dup=dup,signal=5,block=c(TRUE,TRUE,FALSE))
    for(i in 1:dup){
        result=tensor_regress(tsr=data$tsr[[i]],X_covar1 =data$X_covar1, X_covar2 = data$X_covar2,X_covar3= data$X_covar3, core_shape=optimal_rank[s,], Nsim=10, cons="no", lambda = 0.1, alpha = 1, solver ="CG",dist="binary")
        table[s,i]=mean((result$C_ts-data$C_ts)^2) ## tensor method
        naive_C=massive_glm(data$tsr[[i]],data$X_covar3,dist="binary")
        table_naive[s,i]=mean((naive_C-data$C_ts)^2) ## GLM method
    }
}

save(table,table_naive,file="comparison_binary.RData")

