### Code for Figure 8. Comparison between classical GLM and tensor methods in the HCP data anaysis
source("simulation.R")
data=load("presaved/HCP.RData")
############ Method 1: tensor method  ##############
load("presaved/output_HCP.RData")
table(attr[,5])
## age group: 22-25 26-30 31+
## number of individuals: 35  58    43
X=attr[,4:5]
levels(X[,2])=c("22-25","26-30","31+","31+") ## three age groups
contrasts(X[,1]) <- contr.sum
contrasts(X[,2]) <- contr.sum
### intercept, F 1, M,-1; 22-25 [1,0]; 26-30 [0 1]; 30+
coef=array(0,dim=c(68,68,5))
coef[,,1:4]=result$C_ts[,,1:4]##intercept, gender, age 22-25, age 26-30
coef[,,5]=-coef[,,3]-coef[,,4] ## age 30+


############## Method 2: Classical GLM  ##############
tsr=tensor
X_covar3=model.matrix(~as.factor(X[,1])+as.factor(X[,2])) ## baseline female and age 22-25
massive=massive_glm(tsr,X_covar3,"binary")
mass=array(0,dim=c(68,68,5))
mass[,,1]=massive[,,1]
mass[,,2]=massive[,,2]
mass[,,3]=massive[,,3]
mass[,,4]=massive[,,4]
mass[,,5]=-massive[,,4]-massive[,,3]


### Plot histogram to compare two methods###
pdf("compare.pdf",width=15,height=6)
par(mfrow=c(2,4))
hist(c(mass[,,1]),nclass=40,xlab="Effect size",main="Intercept (Classical GLM)")
hist(c(mass[,,2]),nclass=40,xlab="Effect size",main="Gender (Classical GLM)")
hist(c(mass[,,4]),nclass=40,xlab="Effect size",main="Age 26-30 (Classical GLM)")
hist(c(mass[,,5]),nclass=40,xlab="Effect size",main="Age 31+ (Classical GLM)")

hist(c(coef[,,1]),nclass=40,xlab="Effect size",main="Intercept (Tensor regression)")
hist(c(coef[,,2]),nclass=40,xlab="Effect size",main="Gender (Tensor regression)")
hist(c(coef[,,4]),nclass=40,xlab="Effect size",main="Age 26-30 (Tensor regression)")
hist(c(coef[,,5]),xlab="Effect size",nclass=40,main="Age 31+ (Tensor regression)")
dev.off()
