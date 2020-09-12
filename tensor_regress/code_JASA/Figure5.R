## Code for Figure 5. Comparison of mean sqaured prediction error between our method and the other three previous methods
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
source("comparison.R")
library(tensorregress)
library(ggplot2)

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line #40 to run the code ##########

######## plot MSPE for three methods under different (signal, rank) setting #######
load("Figure5.RData")
library(ggplot2)
library(wesanderson)
library(patchwork)

signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))
pdf("Figure5.pdf",width=20,height=3)
s=r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,]),Method=rep(c('GTD (Our method)','HOLRR (Low-rank regression)','HOPLS (Partial least sqaure)','TPG (Projection gradient)'),3),Category=c(rep(1,4),rep(2,4),rep(3,4)))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE"+coord_cartesian(ylim = c(0, max(final[s,r,,])+0.1))) +  labs(title = "Low Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 15))+scale_fill_manual(values=c(rev(wes_palette(n=4, name="Moonrise2"))[1:3],'#C0C0C0'))+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

s=1;r=2;+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,]),Method=rep(c('GTD (Our method)','HOLRR (Low-rank regression)','HOPLS (Partial least sqaure)','TPG (Projection gradient)'),3),Category=c(rep(1,4),rep(2,4),rep(3,4)))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE"+coord_cartesian(ylim = c(0, max(final[s,r,,])+0.1))) +  labs(title = "Low Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 15))+scale_fill_manual(values=c(rev(wes_palette(n=4, name="Moonrise2"))[1:3],'#C0C0C0'))+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

s=2;r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,]),Method=rep(c('GTD (Our method)','HOLRR (Low-rank regression)','HOPLS (Partial least sqaure)','TPG (Projection gradient)'),3),Category=c(rep(1,4),rep(2,4),rep(3,4)))
p3=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE"+coord_cartesian(ylim = c(0, max(final[s,r,,])+0.1))) +  labs(title = "High Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 15))+scale_fill_manual(values=c(rev(wes_palette(n=4, name="Moonrise2"))[1:3],'#C0C0C0'))+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

s=2;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,]),Method=rep(c('GTD (Our method)','HOLRR (Low-rank regression)','HOPLS (Partial least sqaure)','TPG (Projection gradient)'),3),Category=c(rep(1,4),rep(2,4),rep(3,4)))
p4=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE"+coord_cartesian(ylim = c(0, max(final[s,r,,])+0.1))) +  labs(title = "High Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 15))+scale_fill_manual(values=c(rev(wes_palette(n=4, name="Moonrise2"))[1:3],'#C0C0C0'))+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

(p1|p2|p3|p4)
##################### end of plotting ##############################

### If a new run of simulation is desired, please run the code from here and save the results as .RData. Then run the above code to generate figures. ###
set.seed(0)
seed=0
## combinations of signal and rank
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))
##
dup=30;
d=20;
whole_shape = rep(d,3);
final=array(0,dim=c(2,2,3,4))
finalsd=array(0,dim=c(2,2,3,4))
err_res=err_holrr=err_pls=err_tpg=array(0,dim=c(2,2,3,dup))
dist="normal";


for(s in 1:2){ ## signal level
    for(r in 1:2){ ## rank
        for(i in 1:3){ ## number of informative covariates
            for(n in 1:dup){ ## simulation replicates
        core_shape=core_range[r,]
        signal=signal_range[s]
        if(i==1){ ## one mode covaraite
            data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0,0),dist=dist, dup=dup, signal=signal)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }else if(i==2){ ## two mode covaraite
            data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0.4*d,0),dist=dist, dup=dup, signal=signal)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,X_covar2 = data$X_covar2,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }else if (i==3){## three model covariate
            data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d,0.4*d,0.4*d),dist=dist, dup=dup, signal=signal)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,X_covar2 = data$X_covar2,X_covar3 = data$X_covar3,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
        }
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
save(final,finalsd,file="Figure5.RData")
