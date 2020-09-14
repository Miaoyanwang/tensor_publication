### Code for Figure 8. Nations data analysis
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
library(tensorregress)
library(lattice)
library(RColorBrewer)
my.palette <- brewer.pal(n = 12, name = "RdBu")
source("simulation.R")
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
library(lattice)
clustering=kmeans(result$W[[3]],4,nstart=5)
relnames=names(tsr[1,1,])
relnames[clustering$cluster==1]
relnames[clustering$cluster==2]
relnames[clustering$cluster==3]
relnames[clustering$cluster==4]


myPanel <- function(x, y, z,...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  round(M[cbind(x,y)],2)) ## use handy matrix indexing
}

pdf("Figure8.pdf",width=6,height=6)
M=result$C_ts[,,which(relnames=="warning")]
levelplot(M,panel=myPanel,col.regions =my.palette ,at=seq(-1.6, 1.6, length.out=12),main="Relation: warning")


M=result$C_ts[,,which(relnames=="violentactions")]
levelplot(M,panel=myPanel,col.regions =my.palette ,at=seq(-1.6, 1.6, length.out=12),main="Relation: violentactions")


M=result$C_ts[,,which(relnames=="treaties")]
levelplot(M,panel=myPanel,col.regions =my.palette ,at=seq(-1.6, 1.6, length.out=12),main="Relation: treaties")



M=result$C_ts[,,which(relnames=="aidenemy")]
levelplot(M,panel=myPanel,col.regions =my.palette ,at=seq(-1.6, 1.6, length.out=12),main="Relation: aidenemy")


dev.off()
