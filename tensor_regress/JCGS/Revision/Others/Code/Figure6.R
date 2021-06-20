####### Code for Figure 6: Robustness to non-i.i.d. noise ######

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line  to run the code ##########
library(ggplot2)
library(patchwork)

new_color = c("#069AA0","#BCA455","#7B533E")

load("presave/fig6.RData")

#### cor_level vs 1 - Cor
final = 1 - final_c_cor
finalsd = final_sd_c_cor

s=1;r=1; # rank, signal
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"Envelope","SupCP"),3),Category=c(rep(0,3),rep(0.3,3),rep(0.5,3)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"Envelope","SupCP"))
p1=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Correlation level",y="Misallignment error")+coord_cartesian(ylim = c(0, 1)) +  labs(title = "Low Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p1

s=1;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"Envelope","SupCP"),3),Category=c(rep(0,3),rep(0.3,3),rep(0.5,3)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"Envelope","SupCP"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Correlation level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.15)) +  labs(title = "Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10)) +theme(legend.position = "none")
p2

s=2;r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"Envelope","SupCP"),3),Category=c(rep(0,3),rep(0.3,3),rep(0.5,3)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"Envelope","SupCP"))
p3=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Correlation level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.6)) +  labs(title = "Low Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p3

s=2;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,]),sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/sqrt(30),Method=rep(c('STD (Our method)',"Envelope","SupCP"),3),Category=c(rep(0,3),rep(0.3,3),rep(0.5,3)))
data[,3]=factor(data[,3],levels=c('STD (Our method)',"Envelope","SupCP"))
p4=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Correlation level",y="Misallignment error")+coord_cartesian(ylim = c(0, 0.15)) +  labs(title = "High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
  scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))
p4

pdf("Figures/Figure6.pdf", width = 8, height = 3)
(p2|p4)
dev.off()

############# end of plotting ###########

############# loading packages for a new run of simulation #############
install.packages("function/tensorregress_4.0.tar.gz", repos = NULL, type="source")
library(tensorregress)

if(!require(TRES)){ # Envelope
  install.packages(TRES)
  library(TRES)
} 

if(!require(rmatio)){ # write .mat file for SupCP
  install.packages(rmatio)
  library(ramtio)
}

source("function/netglm_code_adjust.R") # GLSNet
source('function/funcs.R')

########### run simulation. Robustness to non-i.i.d. noise #######

D = c(20,20); n = 20; p = 5 ;dup = 30
u_range = rbind(c(3,3), c(4,5))
signal_level_range = c(3,6)
cor_level_range = c(0, 0.3 ,0.5) # 0 = i.i.d.

final_c=final_sd_c=array(0,dim=c(2,2,3,3)) # u, signal, cor_level, method
err_res=err_env=array(0,dim=c(2,2,5,dup))

final_c_cor=final_sd_c_cor=array(0,dim=c(2,2,3,3)) # ours vs envelope
cor_res=cor_env=array(0,dim=c(2,2,5,dup))

tsr_list_i = list(); x_list_i = list(); u_list_i = list();para_list_i = list()
for (i in 1:2) { # u
  tsr_list_j = list(); x_list_j = list(); u_list_j = list();para_list_j = list()
  for (j in 1:2) { # signal
    tsr_list_k = list(); x_list_k = list(); u_list_k = list();  para_list_k = list()
    for (k in 1:3 ) { # cor_level
      tsr_list = list(); x_list = list(); u_list = list();  para_list = list()
      for (m in 1:dup) { # dup
        
        seed = 10000*i + 1000*j + 100*k + m
        set.seed(seed)
        
        u = u_range[i,]; signal_level = signal_level_range[j]; cor_level = cor_level_range[k]
        
        # generate data
        dat = env_sim(D, u, p, n, signal_level, cor_level, dup)
        
        ### fit envelope 
        env_fit = TRR.fit(t(dat$X), dat$tsr[[m]]@data, u, method = "FG")
        err_env = mean((env_fit$fitted.values@data - dat$U@data)^2)
        cor_env = cor(as.vector(env_fit$fitted.values@data), as.vector(dat$U@data))
        
        #### fit tensorregress
        std_fit = tensor_regress(dat$tsr[[m]]@data, X_covar3 = dat$X, core_shape = c(u,p), niter = 10, cons = "non", dist = "normal", initial = "QR_tucker")
        err_res = mean((std_fit$U - dat$U@data)^2)
        cor_res = cor(as.vector(std_fit$U), as.vector(dat$U@data))
        
        # save data to matlab
        # need to transpose the data
        tsr_list[[m]] = tensor_trans(dat$tsr[[m]]@data)
        x_list[[m]] = dat$X
        u_list[[m]] = tensor_trans(dat$U@data)
        para_list[[m]] = c("u", u, "signal", signal_level, "cor_level", cor_level, "dup",m)
        
      }# end dup
      
      final_c[i,j,k,1] = mean(err_res[i,j,k,])
      final_c[i,j,k,2] = mean(err_env[i,j,k,])
      
      final_sd_c[i,j,k,1] = sd(err_res[i,j,k,])
      final_sd_c[i,j,k,2] = sd(err_env[i,j,k,])
      
      final_c_cor[i,j,k,1] = mean(cor_res[i,j,k,])
      final_c_cor[i,j,k,2] = mean(cor_env[i,j,k,])
      
      final_sd_c_cor[i,j,k,1] = sd(cor_res[i,j,k,])
      final_sd_c_cor[i,j,k,2] = sd(cor_env[i,j,k,])
      
      # save data
      tsr_list_k[[k]] = tsr_list; x_list_k[[k]] = x_list; u_list_k[[k]] = u_list; para_list_k[[k]] = para_list 
    }
    # save data
    tsr_list_j[[j]] = tsr_list_k; x_list_j[[j]] = x_list_k; u_list_j[[j]] = u_list_k; para_list_j[[j]] = para_list_k
    
  }
  # save data
  tsr_list_i[[i]] = tsr_list_j; x_list_i[[i]] = x_list_j; u_list_i[[i]] = u_list_j; para_list_i[[i]] = para_list_j
}

sample_list = list(tsr = tsr_list_i, y = x_list_i, x_true = u_list_i, para = para_list_i)
write.mat(sample_list,"mat_data/fig6.mat")

#load("for_chtc/chtc_results/fig6.RData")
# read matlab results
supcp_result = read.mat("mat_output/fig6.mat")

final_c[,,,3] = supcp_result$final_c_sup
final_sd_c[,,,3] = supcp_result$final_sd_c_sup

final_c_cor[,,,3] = supcp_result$final_c_cor_sup
final_sd_c_cor[,,,3] = supcp_result$final_sd_c_cor_sup

save(final_c, final_sd_c,final_c_cor,final_sd_c_cor,file = "presaved/fig6.RData")

