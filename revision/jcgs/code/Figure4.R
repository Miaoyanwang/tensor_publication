###### Mreg vs STD comparison ######

## R package rTensor has been removed from the R CRAN. 
## Install the old version here.
library(remotes)
install_version("rTensor", "1.3")
library("rTensor")

## Mreg code
source("netglm_code_adjust.R") # original data has a typo in the output of best_SymGLM() 
library(tensorregress)


###### run simulation. PMSE vs. number of informative modes
set.seed(0)
seed=0
## combinations of signal and rank
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))
##
dup=10;
d=20;
whole_shape = rep(d,3);
final_mode=array(0,dim=c(2,2,3,2)) # ours vs mreg
final_sd_mode=array(0,dim=c(2,2,3,2))
err_res=err_mreg=array(0,dim=c(2,2,3,dup))
dist="binary";


for(s in 1:2){ ## signal level
  for(r in 1:2){ ## rank
    for(i in 1:3){ ## number of informative covariates
      for(n in 1:dup){ ## simulation replicates
        core_shape=core_range[r,]
        signal=signal_range[s]
        cat("signal = ",signal, ", rank = ", core_shape, ", info number = ",i,", dup = ", n,"\n" )
        if(i==1){ ## one mode covaraite
          data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0,0,0.4*d),dist=dist, dup=dup, signal=signal)
          res = tensor_regress(data$tsr[[n]],X_covar3 = data$X_covar3,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }else if(i==2){ ## two mode covaraite
          data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0,0.4*d,0.4*d),dist=dist, dup=dup, signal=signal)
          res = tensor_regress(data$tsr[[n]],X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }else if (i==3){## three model covariate
          data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0.4*d,0.4*d),dist=dist, dup=dup, signal=signal)
          res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }
        err_res[s,r,i,n]= mean((res$U - data$U)^2)
        
        mreg_res = best_SymGLM(as.tensor(data$tsr[[n]]),data$X_covar3)$result # result with best rank and sparsity
        # calculate the prediction
        Theta = mreg_res$A%*%diag(mreg_res$w)%*%t(mreg_res$A)
        mreg_B = ttm(mreg_res$BB,data$X_covar3,m = 3)@data
        mreg_U = mreg_B
        for (k in 1:dim(mreg_B)[3]) {
          mreg_U[,,k] = mreg_B[,,k] + Theta
        }
        err_mreg[s,r,i,n] = mean((mreg_U - data$U)^2)
      }
      
      final_mode[s,r,i,1]=mean(err_res[s,r,i,])
      final_mode[s,r,i,2]=mean(err_mreg[s,r,i,])
      
      final_sd_mode[s,r,i,1]=sd(err_res[s,r,i,])
      final_sd_mode[s,r,i,2]=sd(err_mreg[s,r,i,])
      
    }
  }
}


### plot figure 4 ####

final = final_mode
finalsd = final_sd_mode

library(ggplot2)
library(patchwork)
new_color = c("#069AA0","#CCC591")

signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))


s=1;r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,])/signal_range[s],sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/signal_range[s],Method=rep(c('STD (Our method)',"Mreg"),3),Category=c(rep(1,2),rep(2,2),rep(3,2)))
data[,3]=factor(data[,3],levels=c("STD (Our method)","Mreg"))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE")+coord_cartesian(ylim = c(0, 0.75)) +  labs(title = "Low Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))


s=1;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,])/signal_range[s],sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/signal_range[s],Method=rep(c('STD (Our method)',"Mreg"),3),Category=c(rep(1,2),rep(2,2),rep(3,2)))
data[,3]=factor(data[,3],levels=c("STD (Our method)","Mreg"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE")+coord_cartesian(ylim = c(0, 0.75)) +  labs(title = "Low Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))

s=2;r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,])/signal_range[s],sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/signal_range[s],Method=rep(c('STD (Our method)',"Mreg"),3),Category=c(rep(1,2),rep(2,2),rep(3,2)))
data[,3]=factor(data[,3],levels=c("STD (Our method)","Mreg"))
p3=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE")+coord_cartesian(ylim = c(0, 0.75)) +  labs(title = "High Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))

s=2;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,])/signal_range[s],sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/signal_range[s],Method=rep(c('STD (Our method)',"Mreg"),3),Category=c(rep(1,2),rep(2,2),rep(3,2)))
data[,3]=factor(data[,3],levels=c("STD (Our method)","Mreg"))
p4=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE")+coord_cartesian(ylim = c(0, 0.75)) +  labs(title = "High Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))

pdf("Mreg_mode.pdf",width=10,height=6)
(p1|p2)/
  (p3|p4)
dev.off()


