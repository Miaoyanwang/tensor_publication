#### code for Figure 3. Convergence of MSE with respect to dimension for three models under different ranks
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
library(tensorregress)
library(ggplot2)

### covariates on all three modes
d_pre=c(30,35,40,45,50,55,60)
p_pre=round(0.4*d_pre)
c_pre=c(2,4,6)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(p_pre,p_pre,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line #83 to run the code ##########

#################### plot MSE for normal model ###############
load("presaved/Figure3_normal.RData")
table=F1_normal

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(dl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

fun.1 <- function(x) max(1.2*res[,4]*res[,1]^2)/x

pdf("Figure3_normal.pdf",width=5,height=4)
figure=ggplot(res, aes(x =d^3/d, y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)+geom_point(size=1)+theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+labs(x="effective sample size",y="mean sqaured error (MSE)")

figure = figure + geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))+coord_cartesian(ylim = c(0,0.6))+stat_function(fun = fun.1,linetype = 2)
figure
dev.off()

################## plot MSE for binary model ###############
load("presaved/Figure3_binary.RData")
table=F1_binary

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(dl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

fun.1 <- function(x) max(1.2*res[,4]*res[,1]^2)/x

pdf("Figure3_binary.pdf",width=5,height=4)
figure=ggplot(res, aes(x =d^3/d, y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)+geom_point(size=1)+theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+labs(x="effective sample size",y="mean sqaured error (MSE)")

figure = figure + geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))+coord_cartesian(ylim = c(0,3.5))+stat_function(fun = fun.1,linetype = 2)
figure
dev.off()

################## plot MSE for poisson model ###############
load("presaved/Figure3_poisson.RData")
table=F1_poisson

MSE_matrix=apply(table[[1]],c(1,2),mean)
sd_matrix=apply(table[[1]],c(1,2),sd)
dl=length(d_pre)
cl=length(c_pre)
res=cbind(rep(d_pre,cl),rep(p_pre,cl),rep(c_pre,rep(dl,cl)),c(MSE_matrix),c(sd_matrix))
res=data.frame(res)
colnames(res)=c("d","p","r","MSE","sd")

fun.1 <- function(x) max(1.2*res[,4]*res[,1]^2)/x

pdf("Figure3_poisson.pdf",width=5,height=4)
figure=ggplot(res, aes(x =d^3/d, y = MSE)) + geom_line(aes(color = as.factor(r)),size = 1)+geom_point(size=1)+theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+labs(x="effective sample size",y="mean sqaured error (MSE)")

figure = figure + geom_errorbar(aes(ymin=MSE-sd,ymax=MSE+sd),width=0.001,position=position_dodge(0.01))+coord_cartesian(ylim = c(0,0.22))+stat_function(fun = fun.1,linetype = 2)
figure
dev.off()

############################## end of plotting ##############################

### If a new run of simulation is desired, please run the code from here and save the results as .RData. Then run the above code to generate figures. ###

set.seed(0)
### covariates on all three modes
dup=30
d_pre=c(30,35,40,45,50,55,60)
p_pre=round(0.4*d_pre)
c_pre=c(2,4,6)
d_range=cbind(d_pre,d_pre,d_pre)
p_range=cbind(p_pre,p_pre,p_pre)
c_range=cbind(c_pre,c_pre,c_pre)

Figure3_normal=conv_rate(seed=NA,signal=10,cons="non",c_range=c_range,dist="normal",alpha=10,dup=dup,d_range=d_range,p_range=p_range,naive=FALSE)
save(Figure3_normal,file="presaved/Figure3_normal.RData")

Figure3_binary=conv_rate(seed=NA,signal=10,cons="non",c_range=c_range,dist="binary",alpha=10,dup=dup,d_range=d_range,p_range=p_range,naive=FALSE)
save(Figure3_binary,file="presaved/Figure3_binary.RData")

Figure3_poisson=conv_rate(seed=NA,signal=10,cons="non",c_range=c_range,dist="poisson",alpha=10,dup=dup,d_range=d_range,p_range=p_range,naive=FALSE)
save(Figure3_poisson,file="presaved/Figure3_poisson.RData")

