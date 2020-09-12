##Code for figure 4. MSE comparison between tensor method vs. massive GLMs.
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
library(tensorregress)
library(ggplot2)

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line #60 to run the code ##########

#################### plot MSE for normal model ###############
load("presaved/Figure4_normal.RData")
mean=c(apply(table,1,mean),apply(table_naive,1,mean))
sd=c(apply(table,1,sd),apply(table_naive,1,sd))
data=cbind(round(mean,2),sd,c(1,1,1,2,2,2),c(5,10,15,5,10,15))
data=data.frame(data)
data[,3]=as.factor(data[,3])
levels(data[,3])=c("Ours","GLM")
names(data)=c("MSE","sd","method","rank")


pdf("Figure4_normal.pdf",width=6,height=4)
p=ggplot(data=data, aes(x=as.factor(rank), y=MSE, fill=as.factor(method))) +
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=.2,
position=position_dodge(.9))+scale_fill_manual(values=c('#069AA0','#C0C0C0'))+labs(x="number of blocks",y="MSE")+coord_cartesian(ylim = c(0, 1.3))+ theme(text = element_text(size = 15))
p
dev.off()


#################### plot MSE for poisson model ###############
load("presaved/Figure4_poisson.RData")
mean=c(apply(table,1,mean),apply(table_naive,1,mean))
sd=c(apply(table,1,sd),apply(table_naive,1,sd))
data=cbind(round(mean,2),sd,c(1,1,1,2,2,2),c(5,10,15,5,10,15))
data=data.frame(data)
data[,3]=as.factor(data[,3])
levels(data[,3])=c("Ours","GLM")
names(data)=c("MSE","sd","method","rank")

pdf("Figure4_poisson.pdf",width=6,height=4)
p=ggplot(data=data, aes(x=as.factor(rank), y=MSE, fill=as.factor(method))) +
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=.2,
position=position_dodge(.9))+scale_fill_manual(values=c('#069AA0','#C0C0C0'))+labs(x="number of blocks",y="MSE")+coord_cartesian(ylim = c(0, 0.9))+ theme(text = element_text(size = 15))
p
dev.off()

#################### plot MSE for binary model ###############
load("presaved/Figure4_binary.RData")
### remove non-convergent solution from GLM; also remove corresponding points from our method for fair comparison
cutoff=quantile(table_naive)[3]+quantile(table_naive)[4]
index=which(table_naive>=cutoff,arr.ind=T)
table[index]=NA
table_naive[index]=NA

mean=c(apply(table,1,function(x) mean(x,na.rm=TRUE)),apply(table_naive,1,function(x) mean(x,na.rm=TRUE)))
sd=c(apply(table,1,function(x)sd(x,na.rm=TRUE)),apply(table_naive,1,function(x) sd(x,na.rm=TRUE)))
data=cbind(round(mean,2),sd,c(1,1,1,2,2,2),c(5,10,15,5,10,15))
data=data.frame(data)
data[,3]=as.factor(data[,3])
levels(data[,3])=c("Ours","GLM")
names(data)=c("MSE","sd","method","rank")

pdf("Figure4_binary.pdf",width=6,height=4)
p=ggplot(data=data, aes(x=as.factor(rank), y=MSE, fill=as.factor(method))) +
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=.2,
position=position_dodge(.9))+scale_fill_manual(values=c('#069AA0','#C0C0C0'))+labs(x="number of blocks",y="MSE")+coord_cartesian(ylim = c(0, 15))+ theme(text = element_text(size = 15))
p
dev.off()

##################### end of plotting ##############################

### If a new run of simulation is desired, please run the code from here and save the results as .RData. Then run the above code to generate figures. ###
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
save(table,table_naive,file="Figure4_normal.RData")

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
save(table,table_naive,file="Figure4_poisson.RData")


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

save(table,table_naive,file="Figure4_binary.RData")

