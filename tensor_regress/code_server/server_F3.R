source("simulation.R")
library(tensorregress)
set.seed(0) ## set seed for this simulation
seed=0

### all covariates
dup=30
d_pre=c(30,35,40,45,50,55,60)
p_pre=round(0.4*d_pre)
c_pre=c(2,4,6)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(p_pre,p_pre,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)


F1_normal=conv_rate(seed=NA,signal=10,cons="non",c_range=c_range,dist="normal",alpha=10,dup=dup,d_range=d_range,p_range=p_range,naive=FALSE)
save(F1_normal,file="F1_normal.RData")

F1_binary=conv_rate(seed=NA,signal=10,cons="non",c_range=c_range,dist="binary",alpha=10,dup=dup,d_range=d_range,p_range=p_range,naive=FALSE)
save(F1_binary,file="F1_binary.RData")

F1_poisson=conv_rate(seed=NA,signal=10,cons="non",c_range=c_range,dist="poisson",alpha=10,dup=dup,d_range=d_range,p_range=p_range,naive=FALSE)
save(F1_poisson,file="F1_poisson.RData")


