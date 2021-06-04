############# tests ##############

install.packages("tensorregress_4.0.tar.gz", repos = NULL, type="source")
library(tensorregress)

library(plot.matrix)
library(RColorBrewer)
library(pracma)


angle_mat = function(A,B_null) {
  # sin(A,B) = || A^T B^{\perp} || (max singular value) 
  mat = t(A) %*% B_null
  sin_theta = svd(mat)$d[1]
  return(sin_theta) 
}

##### Ex1: whether the supervision on third mode affects the factor on first two modes ######
set.seed(0)
d = 10

Core = array(0,dim = c(4,4,4))
for(i in 1:4){
  Core[i,i,i]=2
}
Core[1,1,1] = 1

A = diag(10)[,1:4]
B = diag(10)[,1:4]
C = diag(10)[,1:4]

U = ttl(as.tensor(Core), list(A,B,C), c(1,2,3))@data
X = as.matrix(A[,1])
Y = U + array(rnorm(d^3, sd = 0), dim = c(d,d,d))

### unsupervised
un_res = tensor_regress(Y, X_covar1=NULL,X_covar2=NULL,X_covar3=NULL,c(1,1,1),cons="non",dist="normal", initial = "QR_tucker")

angle_mat(un_res$W$W1, diag(10)[,-c(1,2,3,4)])
angle_mat(un_res$W$W1, diag(10)[,-c(2,3,4)])

angle_mat(un_res$W$W2, diag(10)[,-c(2,3,4)])
angle_mat(un_res$W$W3, diag(10)[,-c(2,3,4)])

### supervised
sup_res = tensor_regress(Y, X_covar1=X,X_covar2=NULL,X_covar3=NULL,c(1,1,1),niter=10,cons="non",dist="normal", initial = "QR_tucker")

angle_mat(X%*%sup_res$W$W1, diag(10)[,-c(1)])
angle_mat(sup_res$W$W2, diag(10)[,-c(2,3,4)])
angle_mat(sup_res$W$W3, diag(10)[,-c(2,3,4)])

angle_mat(sup_res$W$W2, diag(10)[,-c(1)])
angle_mat(sup_res$W$W3, diag(10)[,-c(1)])


###### end Ex 1 ######



##### Ex2: whether our model is robust to overparameterization ######

set.seed(0)
d = 10; alpha = 3

core = array(runif(2^3,-alpha,alpha), dim = c(2,2,2))
W1 = randortho(d)[,1:2]
W2 = randortho(d)[,1:2]
W3 = randortho(2)[,1:2]


X0 = randortho(d)
X = X0[,1:2]

B = ttl(as.tensor(core), list(W1,W2,W3), c(1,2,3))
U = ttm(B, X, 3)@data
Y = U + array(rnorm(d^3, sd = 0.1), dim = c(d,d,d))

# regular
res = tensor_regress(Y, X_covar1 = NULL,X_covar2=NULL,X_covar3=X,c(2,2,2),niter = 10,cons="non",dist="normal", initial = "QR_tucker")
plot(res$C_ts[,,1], col=brewer.pal(n = 11, name = "RdBu"), breaks = seq(-2, 2, length.out=12),border = NA)
# overparameterizaiton
X1 = X0[,1:3]
res1 = tensor_regress(Y, X_covar1=NULL,X_covar2=NULL,X_covar3=X1,c(2,2,3),niter = 10,cons="non",dist="normal", initial = "QR_tucker")
plot(res1$C_ts[,,1], col=brewer.pal(n = 11, name = "RdBu"), breaks = seq(-2, 2, length.out=12),border = NA)

mean((res$C_ts[,,1] - res1$C_ts[,,1])^2, na.rm=FALSE)

