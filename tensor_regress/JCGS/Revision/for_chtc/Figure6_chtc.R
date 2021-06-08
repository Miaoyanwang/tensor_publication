######## Figure 6 (CHTC version) #######
rm(list = ls())

install.packages("tensorregress_4.0.tar.gz", repos = NULL, type="source")
library(tensorregress)

if(!require(TRES)){ # Envelope
  install.packages("TRES")
  library(TRES)
} 

if(!require(pracma)){ # Envelope
  install.packages("pracma")
  library(pracma)
} 

if(!require(rmatio)){ # write .mat file for SupCP
  install.packages("rmatio")
  library(ramtio)
}

source("funcs.R") # useful functions


args <- commandArgs(T) # env_dim, signal_level, cor_level, dup
i = as.numeric(args[1]); j = as.numeric(args[2]); 
k = as.numeric(args[3]); m = as.numeric(args[4])

seed = 10000*i + 1000*j + 100*k + m
set.seed(seed)

D = c(20,20); n = 20; p = 5 ;dup = 30
u_range = rbind(c(3,3), c(4,5))
signal_level_range = c(3,6)
cor_level_range = c(0, 0.3 ,0.5) # 0 = i.i.d.

err_res=err_env=c()
cor_res=cor_env=c()

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

result = data.frame(err = c(err_res,err_env), cor = c(cor_res,cor_env))
save(result, file = paste0("i_",i,"_j_",j,"_k_",k,"_m_",m,"_fig6.RData") )

# save data to matlab
# need to transpose the data
tsr = tensor_trans(dat$tsr[[m]]@data)
x = dat$X
u = tensor_trans(dat$U@data)
para = c("u", u, "signal", signal_level, "cor_level", cor_level, "dup",m)

sample_list = list(tsr = tsr, y = x, x_true = u, para = para)
write.mat(sample_list,paste0("i_",i,"_j_",j,"_k_",k,"_m_",m,"_fig6.mat"))