######### Figure 4 (CHTC version) ########
rm(list = ls())

install.packages("tensorregress_4.0.tar.gz", repos = NULL, type="source")
library(tensorregress)

if(!require(TRES)){ # Envelope
  install.packages(TRES)
  library(TRES)
} 

if(!require(rrpack)){ # mRRR
  install.packages(rrpack)
  library(rrpack)
} 

if(!require(rmatio)){ # write .mat file for SupCP
  install.packages(rmatio)
  library(ramtio)
}

source("netglm_code_adjust.R") # GLSNet

source("funcs.R") # useful functions

########### run simulation. 1 - Correlation vs effective sample size #######

####### Gaussian #######
args <- commandArgs(T) # signal, rank, dim, dup
s = as.numeric(args[1]); r = as.numeric(args[2]); 
i = as.numeric(args[3]); n = as.numeric(args[4])

seed = 10000*s + 1000*r + 100*i + n
set.seed(seed)

## combinations of signal and rank
d=c(20,40,60,80,100); dup = 30
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))

err_res=err_env=err_rrr=err_mreg1=array(0,dim=c(2,2,5,dup))

cor_res=cor_env=cor_rrr=cor_mreg1 =array(0,dim=c(2,2,5,dup))

dist="normal"

whole_shape = c(20,20,d[i]); 
core_shape=core_range[r,]
signal=signal_range[s]

data=sim_data(seed,whole_shape = whole_shape, core_shape=core_shape,p=c(0,0,0.4*d[i]),dist=dist, dup=dup, signal=signal, ortho = T)
res = tensor_regress(data$tsr[[n]],X_covar3 = data$X_covar3,core_shape=core_shape,niter=10, cons = 'non', dist = dist, initial = "QR_tucker")### tensor_regress
        
# tensorregress
err_res= mean((res$U - data$U)^2)
cor_res= cor(as.vector(res$U), as.vector(data$U))

# envelope
env_fit = TRR.fit(t(data$X_covar3), data$tsr[[n]], core_shape[1:2], method = "FG")

err_env= mean((env_fit$fitted.values@data - data$U)^2)
cor_env= cor(as.vector(env_fit$fitted.values@data), as.vector(data$U))


###### grid search range for GLSNet
rank_range = c(2,3,4,5,6); sparsity_range = (d^3)*0.4*c(0.1,0.3,0.5,0.7,0.9)

# Mreg1
mreg_res = search_GLM(as.tensor(data$tsr[[n]]), data$X_covar3, true_U = data$U, rank_range,sparsity_range, dist = "normal")

err_mreg1 = mreg_res$err_opt
cor_mreg1 = mreg_res$cor_opt

# rrr
# unfold the tensor on mode 3
data0 = unfold(as.tensor(data$tsr[[n]]), row_idx = 3, col_idx = c(1,2))@data
U0 = unfold(as.tensor(data$U), row_idx =  3, col_idx = c(1,2))@data

# rrreg fit
#rrr_fit = rrr.fit(data0, X = data$X_covar3, nrank = core_shape[3])
rrr_fit = mrrr(data0, X = data$X_covar3, family = list(gaussian()),penstr = list(penaltySVD = "rankCon", lambdaSVD = core_shape[3]))

err_rrr= mean((rrr_fit$mu - U0)^2)
cor_rrr= cor(as.vector(rrr_fit$mu), as.vector(U0))

result = data.frame(err = c(err_res,err_env,err_mreg1,err_rrr), cor = c(cor_res,cor_env,cor_mreg1,cor_rrr))
save(result, file = paste0("s_",s,"_r_",r,"_i_",i,"_n_",n,"_fig4_sample.RData") )

# save data to matlab
# need to transpose the data
tsr = tensor_trans(data$tsr[[n]])
x = data$X_covar3
u= tensor_trans(data$U)
para = c("signal", signal, "rank",core_shape,"dim",d[i],"dup",n)

sample_list = list(tsr = tsr, y = x, x_true = u, para = para)
write.mat(sample_list,paste0("s_",s,"_r_",r,"_i_",i,"_n_",n,"_fig4_sample.mat"))
