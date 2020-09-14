source("simulation.R")
library(tensorregress)

set.seed(0) ## set seed for this simulation
seed=0

### all covariates included ##
d=c(20,40)
rank_range=array(0,dim=c(3,10,3))
rank_range[1,,]=rbind(c(3,3,3),c(3,3,2),c(3,2,2),c(2,2,2),c(3,3,4),c(3,4,4),c(4,4,4),c(2,2,4),c(3,2,3),c(2,3,2))
rank_range[2,,]=rbind(c(4,4,6),c(4,4,5),c(3,3,6),c(4,5,7),c(3,4,6),c(4,5,6),c(5,5,5),c(3,3,5),c(3,3,4),c(3,4,7))
rank_range[3,,]=rbind(c(6,8,8),c(6,7,8),c(5,8,8),c(7,8,8),c(5,7,7),c(5,7,8),c(6,7,7),c(5,6,7),c(5,6,6),c(5,5,5))
core_shape=rbind(c(3,3,3),c(4,4,6),c(6,8,8))

########################## normal model ##########################
dist="normal"
dup=30
final=array(0,dim=c(3,2,dup,3))
for(s in 1:3){
    for(i in 1:2){
        whole_shape = c(d[i],d[i],d[i])
        p=0.4*whole_shape
        data=sim_data(seed, whole_shape = whole_shape, core_shape=core_shape[s,],p=p,dist=dist, dup=dup, signal=4)
        table=NULL
        
        for(j in 1:dup){
            res=sele_rank(data$tsr[[j]],X_covar1 = data$X_covar1, X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,rank_range=rank_range[s,,],Nsim=10,cons = 'non',dist=dist)
            table=rbind(table,res$rank)
        }
        
        final[s,i,,]=table
    }
}

save(final,file="normal_rank.RData")


########################## poisson model ##########################
dist="poisson"
final=array(0,dim=c(3,2,dup,3))
for(s in 1:3){
    for(i in 1:2){
        whole_shape = c(d[i],d[i],d[i])
        p=0.4*whole_shape
        data=sim_data(seed, whole_shape = whole_shape, core_shape=core_shape[s,],p=p,dist=dist, dup=dup, signal=4)
        table=NULL
        
        for(j in 1:dup){
            res=sele_rank(data$tsr[[j]],X_covar1 = data$X_covar1, X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,rank_range=rank_range[s,,],Nsim=10,cons = 'non',dist=dist)
            table=rbind(table,res$rank)
        }
        final[s,i,,]=table
    }
}

save(final,file="poisson_rank.RData")

## d = 20
apply(final[1,1,,],2,mean) ## ground truth (3,3,3)
apply(final[2,1,,],2,mean) ## ground truth (4,4,6)
apply(final[3,1,,],2,mean) ## ground truth (6,8,8)
apply(final[1,1,,],2,sd) ## ground truth (3,3,3)
apply(final[2,1,,],2,sd) ## ground truth (4,4,6)
apply(final[3,1,,],2,sd) ## ground truth (6,8,8)

## d = 40
apply(final[1,2,,],2,mean) ## ground truth (3,3,3)
apply(final[2,2,,],2,mean) ## ground truth (4,4,6)
apply(final[3,2,,],2,mean) ## ground truth (6,8,8)
apply(final[1,2,,],2,sd) ## ground truth (3,3,3)
apply(final[2,2,,],2,sd) ## ground truth (4,4,6)
apply(final[3,2,,],2,sd) ## ground truth (6,8,8)
