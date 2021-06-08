####### Code for Figure 7: Robustness to sparsity ######

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line  to run the code ##########
library(ggplot2)
library(patchwork)

new_color = c("#069AA0","#D6CFC4")

load("presaved/fig7.RData")

#### cor_level vs 1 - Cor
final = 1 - final_s_cor
finalsd = final_sd_s_cor

s=1;r=1; # signal, rank
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"GLSNet"),3),Category=c(rep(0,2),rep(0.1,2),rep(0.9,2)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"GLSNet"))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Sparsity level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.1)) +  labs(title = "Low Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))+theme(legend.position = "none")
p1

s=1;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"GLSNet"),3),Category=c(rep(0,2),rep(0.1,2),rep(0.9,2)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"GLSNet"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Sparsity level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.15)) +  labs(title = "Low Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10)) +theme(legend.position = "none")
p2

s=2;r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"GLSNet"),3),Category=c(rep(0,2),rep(0.1,2),rep(0.9,2)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"GLSNet"))
p3=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Sparsity level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.1)) +  labs(title = "High Signal",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p3

s=2;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"GLSNet"),3),Category=c(rep(0,2),rep(0.1,2),rep(0.9,2)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"GLSNet"))
p4=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Sparsity level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.15)) +  labs(title = "High Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p4

pdf("Figures/Figure7.pdf", width = 8, height = 3)
(p1|p3)
dev.off()

############# loading packages for a new run of simulation #############
install.packages("function/tensorregress_4.0.tar.gz", repos = NULL, type="source")
library(tensorregress)

source("function/netglm_code_adjust.R") # GLSNet
source('function/funcs.R')

########### run simulation. Robustness to sparsity #######

d = 20; d3 = 20; p = 8; dup = 30
signal_range = c(2,4); r_range = c(2,4); s_range = c(0.1,0.5,0.9)

final_s=final_sd_s=array(0,dim=c(2,2,3,2)) # signal, rank, sparse_level, method
err_res=err_mreg=array(0,dim=c(2,2,5,dup))

final_s_cor=final_sd_s_cor=array(0,dim=c(2,2,3,2)) # ours vs glsnet
cor_res=cor_mreg=array(0,dim=c(2,2,5,dup))

for (i in 1:2) { # signal
  for (j in 1:2) { # rank
    for (k in 1:3) { # sparsity
      
      signal = signal_range[i]; R = r3 = r_range[j]; s = s_range[k]
      
      for (m in 1:dup) {
        
        seed = 10000*i + 1000*j + 100*k + m
        set.seed(seed)
        
        dat = sparse_sim(d, d3, p, R, r3, s, signal, dup)
        
        #### fit tensorregress
        std_fit = tensor_regress(dat$tsr[[m]], X_covar3 = dat$X, core_shape = c(R,R,r3), niter = 20, cons = "non", dist = "binary", initial = "QR_tucker")
        err_res = mean((std_fit$U - dat$U)^2)
        cor_res = cor(as.vector(std_fit$U), as.vector(dat$U))
        
        ##### fit GLSNet
        mreg_res = mySymGLM(as.tensor(dat$tsr[[m]]), dat$X, R, sparsity=d*d*p*s, niter=50)
        
        Theta = mreg_res$A%*%diag(mreg_res$w)%*%t(mreg_res$A)
        mreg_B = ttm(mreg_res$BB,dat$X,m = 3)@data
        mreg_U = mreg_B
        for (l in 1:dim(mreg_B)[3]) {
          mreg_U[,,l] = mreg_B[,,l] + Theta
        }
        err_mreg =  mean((mreg_U - dat$U)^2)
        cor_mreg = cor(as.vector(mreg_U), as.vector(dat$U))
      } # end dup
      
      final_s[i,j,k,1] = mean(err_res[i,j,k,])
      final_s[i,j,k,2] = mean(err_mreg[i,j,k,])
      
      final_sd_s[i,j,k,1] = sd(err_res[i,j,k,])
      final_sd_s[i,j,k,2] = sd(err_mreg[i,j,k,])
      
      final_s_cor[i,j,k,1] = mean(cor_res[i,j,k,])
      final_s_cor[i,j,k,2] = mean(cor_mreg[i,j,k,])
      
      final_sd_s_cor[i,j,k,1] = sd(cor_res[i,j,k,])
      final_sd_s_cor[i,j,k,2] = sd(cor_mreg[i,j,k,])
      
    }
    
  }
  
}

save(final_s, final_sd_s,final_s_cor,final_sd_s_cor,file = "presaved/fig7.RData")