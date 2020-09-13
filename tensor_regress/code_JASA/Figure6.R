### Code for Figure 7. HCP data analysis
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
library(tensorregress)

data=load("presaved/HCP.RData")
tsr=tensor
X=attr[,4:5]
levels(X[,2])=c("22-25","26-30","31+","31+") ## three groups
contrasts(X[,1]) <- contr.sum
contrasts(X[,2]) <- contr.sum
X_covar3=model.matrix(~as.factor(X[,1])+as.factor(X[,2]))

core_shape=c(10,10,4)
result=tensor_regress(tsr,X_covar1=NULL,X_covar2=NULL,X_covar3,core_shape,Nsim=50,cons="non",lambda=0.1,alpha=10,solver="CG",dist="binary")

save(result,file="output_HCP.RData")
