#### test new pack #####

library(tensorregress)
library(MASS)
library(pracma)
library(tictoc)
library(ggplot2)
library(patchwork)

### data 
set.seed(0054)
dup = 10; d=20; whole_shape = rep(d,3); dist = "normal"
r_range = rbind(c(3,3,3), c(4,5,6))
signal_range = c(3,6)

err0 = err1 = err2 = err3 = array(0,dim = c(2,2,dup)) # signal, rank, metric, dup. metric: mse, 1 - cor, time 
final_err = final_sd_err = array(0, dim = c(2,2,4)) # signal, rank, method

cor0 = cor1 = cor2 = cor3 = array(0,dim = c(2,2,dup)) # signal, rank, metric, dup. metric: mse, 1 - cor, time 
final_cor = final_sd_cor = array(0, dim = c(2,2,4)) # signal, rank, method

times0 = times1 = times2 =times3 = array(0,dim = c(2,2,dup)) # signal, rank, metric, dup. metric: mse, 1 - cor, time 
final_time = final_sd_time = array(0, dim = c(2,2,4)) # signal, rank, method

for (i in 1:2) {  # signal
  for (j in 1:2) {  # rank
    for (n in 1:dup) {
      
      #i  = 2; j = 2; n = 3
      
      signal = signal_range[i]
      core_shape = r_range[j,]
      
      cat("signal = ",signal, ", rank = ", core_shape,", dup = ", n,"\n" )
      
      data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0,0,0.4*d),dist=dist, dup=dup, signal=signal)
      
      ptm = proc.time()
      res = tensor_regress(data$tsr[[n]],X_covar3 = data$X_covar3,
                           core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
      ptime = proc.time() - ptm
      times0[i,j,n] = ptime[2]
      err0[i,j,n] = mean((res$U - data$U)^2)
      cor0[i,j,n] = cor(as.vector(res$U), as.vector(data$U))
      
      
      ptm = proc.time()
      res1 = tensor_regress1(data$tsr[[n]],X_covar3 = data$X_covar3,
                             core_shape=core_shape,Nsim=10, cons = 'non', dist = dist,
                             initial = "de_tucker", alg = "alter")
      ptime = proc.time() - ptm
      times1[i,j,n] = ptime[2]
      err1[i,j,n] = mean((res1$U - data$U)^2)
      cor1[i,j,n] = cor(as.vector(res1$U), as.vector(data$U))
      
      ptm = proc.time()
      res2 = tensor_regress1(data$tsr[[n]],X_covar3 = data$X_covar3,
                             core_shape=core_shape,Nsim=10, cons = 'non', dist = dist,
                             initial = "de_tucker", alg = "unsup")
      ptime = proc.time() - ptm
      times2[i,j,n] = ptime[2]
      err2[i,j,n] = mean((res2$U - data$U)^2)
      cor2[i,j,n] = cor(as.vector(res2$U), as.vector(data$U))
      
      ptm = proc.time()
      res3 = tensor_regress1(data$tsr[[n]],X_covar3 = data$X_covar3,
                             core_shape=core_shape,Nsim=10, cons = 'non', dist = dist,
                             initial = "tucker", alg = "unsup")
      ptime = proc.time() - ptm
      times3[i,j,n] = ptime[2]
      err3[i,j,n] = mean((res3$U - data$U)^2)
      cor3[i,j,n] = cor(as.vector(res3$U), as.vector(data$U))
      
    }
    
    final_err[i,j,1] = mean(err0[i,j,]); final_err[i,j,2] = mean(err1[i,j,]); final_err[i,j,3] = mean(err2[i,j,]); final_err[i,j,4] = mean(err3[i,j,])
    final_sd_err[i,j,1] = sd(err0[i,j,]); final_sd_err[i,j,2] = sd(err1[i,j,]); final_sd_err[i,j,3] = sd(err2[i,j,]); final_sd_err[i,j,4] = sd(err3[i,j,])
    
    final_cor[i,j,1] = mean(cor0[i,j,]); final_cor[i,j,2] = mean(cor1[i,j,]); final_cor[i,j,3] = mean(cor2[i,j,]);final_cor[i,j,4] = mean(cor3[i,j,])
    final_sd_cor[i,j,1] = sd(cor0[i,j,]); final_sd_cor[i,j,2] = sd(cor1[i,j,]); final_sd_cor[i,j,3] = sd(cor2[i,j,]); final_sd_cor[i,j,4] = sd(cor3[i,j,])
    
    final_time[i,j,1] = mean(times0[i,j,]); final_time[i,j,2] = mean(times1[i,j,]); final_time[i,j,3] = mean(times2[i,j,]); final_time[i,j,4] = mean(times3[i,j,])
    final_sd_time[i,j,1] = sd(times0[i,j,]); final_sd_time[i,j,2] = sd(times1[i,j,]); final_sd_time[i,j,3] = sd(times2[i,j,]); final_sd_time[i,j,4] = sd(times3[i,j,])
  }
}


save(final_err,final_sd_err, final_cor, final_sd_cor, final_time, final_sd_time, file = "new_pack_check.RData")


new_color = c("#069AA0","#CCC591","#BCA455","#D6CFC4")

#### rank vs mse
final = final_err
finalsd = final_sd_err

s=1;
data=data.frame(PMSE=c(final[s,1,],final[s,2,]),sd=c(finalsd[s,1,],finalsd[s,2,]),Method=rep(c('Origin',"de_Tucker","de_Tucker_un", "Tucker"),2),Category=c(rep("low",4),rep("high",4)))
data[,3]=factor(data[,3],levels=c('Origin',"de_Tucker","de_Tucker_un", "Tucker"))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Rank",y="MSE")+coord_cartesian(ylim = c(0, 0.06)) +  labs(title = "Low Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p1

s=2;
data=data.frame(PMSE=c(final[s,1,],final[s,2,]),sd=c(finalsd[s,1,],finalsd[s,2,]),Method=rep(c('Origin',"de_Tucker","de_Tucker_un", "Tucker"),2),Category=c(rep("low",4),rep("high",4)))
data[,3]=factor(data[,3],levels=c('Origin',"de_Tucker","de_Tucker_un", "Tucker"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Rank",y="MSE")+coord_cartesian(ylim = c(0, 0.06)) +  labs(title = "High Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p2


pdf("newpack_err.pdf",height = 3, width = 8)
(p1|p2)
dev.off()

#### rank vs 1 - correlation
final = 1- final_cor
finalsd = final_sd_cor

s=1;
data=data.frame(PMSE=c(final[s,1,],final[s,2,]),sd=c(finalsd[s,1,],finalsd[s,2,]),Method=rep(c('Origin',"de_Tucker","de_Tucker_un", "Tucker"),2),Category=c(rep("low",4),rep("high",4)))
data[,3]=factor(data[,3],levels=c('Origin',"de_Tucker","de_Tucker_un", "Tucker"))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Rank",y="1 - Correlation")+coord_cartesian(ylim = c(0, 0.3)) +  labs(title = "Low Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p1

s=2;
data=data.frame(PMSE=c(final[s,1,],final[s,2,]),sd=c(finalsd[s,1,],finalsd[s,2,]),Method=rep(c('Origin',"de_Tucker","de_Tucker_un", "Tucker"),2),Category=c(rep("low",4),rep("high",4)))
data[,3]=factor(data[,3],levels=c('Origin',"de_Tucker","de_Tucker_un", "Tucker"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Rank",y="1 - Correlation")+coord_cartesian(ylim = c(0, 0.1)) +  labs(title = "High Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p2


pdf("newpack_cor.pdf",height = 3, width = 8)
(p1|p2)
dev.off()

#### rank vs time
final = final_time
finalsd = final_sd_time

s=1;
data=data.frame(PMSE=c(final[s,1,],final[s,2,]),sd=c(finalsd[s,1,],finalsd[s,2,]),Method=rep(c('Origin',"de_Tucker","de_Tucker_un", "Tucker"),2),Category=c(rep("low",4),rep("high",4)))
data[,3]=factor(data[,3],levels=c('Origin',"de_Tucker","de_Tucker_un", "Tucker"))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Rank",y="time")+coord_cartesian(ylim = c(0, 0.4)) +  labs(title = "Low Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p1

s=2;
data=data.frame(PMSE=c(final[s,1,],final[s,2,]),sd=c(finalsd[s,1,],finalsd[s,2,]),Method=rep(c('Origin',"de_Tucker","de_Tucker_un", "Tucker"),2),Category=c(rep("low",4),rep("high",4)))
data[,3]=factor(data[,3],levels=c('Origin',"de_Tucker","de_Tucker_un", "Tucker"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Rank",y="time")+coord_cartesian(ylim = c(0, 0.35)) +  labs(title = "High Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p2


pdf("newpack_time.pdf",height = 3, width = 8)
(p1|p2)
dev.off()


######### test
# original results
dat=sim_data(whole_shape = whole_shape, core_shape=c(4,5,6),p=c(0,0,0.4*d),dist=dist, dup=dup, signal=6)

ptm = proc.time()
res = tensor_regress(dat$tsr[[1]],X_covar1 = dat$X_covar1,X_covar2 = dat$X_covar2,X_covar3 = dat$X_covar3,
                     core_shape=c(4,5,6),Nsim=10, cons = 'non', dist="normal")### tensor_regress
ptime = proc.time() - ptm
ptime[2]

ptm = proc.time()
res1 = tensor_regress1(dat$tsr[[1]],X_covar1 = dat$X_covar1,X_covar2 = dat$X_covar2,X_covar3 = dat$X_covar3,
                       core_shape=c(3,3,3),Nsim=10, cons = 'non', dist="normal",
                       initial = "tucker", alg = "alter")### tensor_regress
ptime = proc.time() - ptm
ptime[2]

ptm = proc.time()
res2 = tensor_regress1(dat$tsr[[1]],X_covar1 = dat$X_covar1,X_covar2 = dat$X_covar2,X_covar3 = dat$X_covar3,
                       core_shape=c(3,3,3),Nsim=10, cons = 'non', dist="normal",
                       initial = "tucker", alg = "unsup")### tensor_regress
ptime = proc.time() - ptm
ptime[2]

class(ptime[2]) 

cor(as.vector(res$U), as.vector(dat$U))

mean((res$U - dat$U)^2)
mean((res1$U - dat$U)^2)
mean((res2$U - dat$U)^2)

cor(as.vector(res$U), as.vector(dat$U))
cor(as.vector(res1$U), as.vector(dat$U))
cor(as.vector(res2$U), as.vector(dat$U))
