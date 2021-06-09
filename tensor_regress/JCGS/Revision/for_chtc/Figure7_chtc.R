######## Figure 7 (CHTC version) #######
rm(list = ls())

install.packages("tensorregress_4.0.tar.gz", repos = NULL, type="source")
library(tensorregress)

source("netglm_code_adjust.R") # GLSNet

source("funcs.R") # useful functions


args <- commandArgs(T) # signal, rank, sparsity, dup
i = as.numeric(args[1]); j = as.numeric(args[2]); 
k = as.numeric(args[3]); m = as.numeric(args[4])

#i = 1; j = 2; k = 3; m = 30
seed = 10000*i + 1000*j + 100*k + m
set.seed(seed)

d = 20; d3 = 20; p = 8; dup = 30
signal_range = c(2,4); r_range = c(2,4); s_range = c(0.1,0.5,0.9)

err_res=err_mreg=c()
cor_res=cor_mreg=c()

signal = signal_range[i]; R = r3 = r_range[j]; s = s_range[k]
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

result = data.frame(err = c(err_res,err_mreg), cor = c(cor_res,cor_mreg))
save(result, file = paste0("i_",i,"_j_",j,"_k_",k,"_m_",m,"_fig7.RData") )