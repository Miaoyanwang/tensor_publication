### Code for Table 3. Nations data analysis
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
library(tensorregress)
set.seed(0)
seed=0
data=load("presaved/nations.RData")
tsr=R
tsr[is.na(tsr)]=0
X=cov

#rank_range=valid_rank(expand.grid(c(3:5),c(3:5),c(3:5)))
#res=sele_rank(tsr,X_covar1=cov,X_covar2=cov,X_covar3=NULL,rank_range=rank_range,Nsim=10,cons = 'non',dist="binary")


core_shape=c(4,4,4)
result=tensor_regress(tsr,X_covar1=cov,X_covar2=cov,X_covar3=NULL,core_shape,Nsim=10,cons="penalty",lambda=0.1,alpha=10,solver="CG",dist="binary")
save(result,file="output_nations.RData")
clustering=kmeans(result$W[[3]],4,nstart=5)
relnames=names(tsr[1,1,])
relnames[clustering$cluster==1]
relnames[clustering$cluster==2]
relnames[clustering$cluster==3]
relnames[clustering$cluster==4]

