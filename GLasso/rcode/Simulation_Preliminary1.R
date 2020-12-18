# Simulation Preliminary

library(MASS)
library(ggplot2)
#### Source main functions ####
source("falsePN.R")
source("graphicalLasso.R")
source("GLasso.R")

#### Bricks ####
emp_cov_matrix = function(X,mu) {
  
  # this function calculate the empirical covariance matrix
  
  p = dim(X)[2]
  emp_cov = matrix(0,nrow = p,ncol = p)
  for (i in 1:dim(X)[1]) {
    X_v = as.matrix(X[i,])
    emp_cov = emp_cov + (X_v-mu)%*%t(X_v-mu)
  }
  
  return(emp_cov/dim(X)[1])
  
}

cal_mean_PN = function(Theta,Sigma,duplicate,n,p){
  # calculate the mean FN or FP with duplication given Sigam
  dup_FN = rep(0,duplicate)
  dup_FP = rep(0,duplicate)
  dup_tune = rep(0,duplicate)
  
  for (dup in 1:duplicate) {
    X = mvrnorm(n, mu = rep(0,p), Sigma)
    S = emp_cov_matrix(X, mu = rep(0,p))
    
    glasso = GLasso(S,n,5)

    dup_tune[dup] = glasso$rho_opt
    PN_rate = falsePN(Theta,glasso$Omega)
    dup_FN[dup] = PN_rate$FN
    dup_FP[dup] = PN_rate$FP
  }
  
  return(list(tune = mean(dup_tune), FN = mean(dup_FN), FP = mean(dup_FP)))
}


#### Designed simulation experiments ####

# We focus on the false positive/negative rate under different network settings
# We consider 5 settings: sparse, dense (in Tibshirani's paper), 
# chain, nearst-neighbourhood, and scale-free (in Guo's paper)

# basic settings
#p_seq = seq(30,150,30)
n = 200 # number of sample for one parameter
duplicate = 3

p_seq = seq(20,120,20)


# initialization
FN_matrix = matrix(0, nrow = 5,ncol = length(p_seq))
FP_matrix = matrix(0, nrow = 5,ncol = length(p_seq))
tune_para = matrix(0, nrow = 5,ncol = length(p_seq))


# simulation
set.seed(2020)
for (p_iter in 1:length(p_seq)) {
  
  # for each p, generate the true parameter and data
  p = p_seq[p_iter]
  
  ## Sparse setting ##
  # Theta_ii = 1, Theta_{i,i-1} = Theta_{i-1,i} = 0.5, others are 0.
  
  Theta_sparse = matrix(0,nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i == j){ # Theta_ii
        Theta_sparse[i,j] = 1
      }
      if(i == (j-1)){# Theta_i-1,i
        Theta_sparse[i,j] = 0.5
      }
      if(i == (j+1)){# Theta_i,i-1
        Theta_sparse[i,j] = 0.5
      }
    }# end for j
  }# end for i
  
  Sigma_sparse = solve(Theta_sparse)
  glasso = cal_mean_PN(Theta_sparse,Sigma_sparse,duplicate,n,p)
  
  tune_para[1,p_iter] = glasso$tune
  FN_matrix[1,p_iter] = glasso$FN
  FP_matrix[1,p_iter] = glasso$FP
  
  cat("Finish sparse setting, p = ",p,"\n")
  ## Dense setting ##
  # Theta_ii = 2, Theta_ii' = 1
  
  Theta_dense = matrix(0,nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i == j){ # Theta_ii
        Theta_dense[i,j] = 2
      }else{
        Theta_dense[i,j] = 1
      }
    }# end for j
  }# end for i
  
  Sigma_dense = solve(Theta_dense)
  glasso = cal_mean_PN(Theta_dense,Sigma_dense,duplicate,n,p)
  
  tune_para[2,p_iter] = glasso$tune
  FN_matrix[2,p_iter] = glasso$FN
  FP_matrix[2,p_iter] = glasso$FP
  
  cat("Finish dense setting, p = ",p,"\n")
  
  ## Chain network setting ##
  # sigma_jj' = exp(-|sj-sj'|/2), s1<...<sp, sj-sj-1 \sim U(0.5,1)
  dup_FN = rep(0,duplicate)
  dup_FP = rep(0,duplicate)
  dup_tune = rep(0,duplicate)
  
  for (dup in 1:duplicate) {
    Sigma_chain = matrix(0,nrow = p, ncol = p)
    s_seq = cumsum(runif(p,0.5,1))
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma_chain[i,j] = exp(-abs(s_seq[i] - s_seq[j])/2)
      }# end for j
    }# end for i
    
    Theta_chain = solve(Sigma_chain)
    Theta_chain[abs(Theta_chain) <= 1e-16] = 0
    
    glasso_chain = cal_mean_PN(Theta_chain,Sigma_chain,1,n,p)
    
    dup_tune[dup] = glasso_chain$tune
    dup_FN[dup] = glasso_chain$FN
    dup_FP[dup] = glasso_chain$FP
  } # end for dup
  
  tune_para[3,p_iter] = mean(dup_tune)
  FN_matrix[3,p_iter] = mean(dup_FN)
  FP_matrix[3,p_iter] = mean(dup_FP)
  
  cat("Finish chain network setting, p = ",p,"\n")
  
  
  ## Nearest-neighbour network setting ##
  # generate p points (xi,yi) in the square [0,1]x[0,1], 
  # calculate distance between each pair of points
  # m-nearest points are connected with omega \sim [-1,-0.5]cup[0.5, 1] uniformly
  
  dup_FN = rep(0,duplicate)
  dup_FP = rep(0,duplicate)
  dup_tune = rep(0,duplicate)
  
  for (dup in 1:duplicate) {
  
    x = runif(p)
    y = runif(p)
    distance = matrix(Inf, nrow = p,ncol = p)
    Theta_mn = matrix(0,nrow = p, ncol = p)
    
    # calculate distance
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        distance[i,j] = sqrt( (x[i]-x[j])^2 + (y[i]-y[j])^2 )
      }# end for j
    }# end for i
    
    distance[lower.tri(distance)] = t(distance)[lower.tri(distance)] 
    
    # find m nearest
    for (k in 1:p) {
      ind = order(distance[k,])[1:3] # 3-nearset neighbour
      ind = ind[ind > k] # to avoid repeatence
      
      for(l in ind){
        Theta_mn[k,l] = sign(runif(1,-1,1))*runif(1,0.5,1)
      }# end for l
      
    }#end for k
    
    Theta_mn[lower.tri(Theta_mn)] = t(Theta_mn)[lower.tri(Theta_mn)]
    #Theta_mn = Theta_mn + diag(rep(1,p))
    Theta_mn = Theta_mn + diag( colSums(abs(Theta_mn)) + rep(1,p) ) # ensure Theta_mn is invertible and Sigma is pos-definite
    Theta_mn[abs(Theta_mn) <= 1e-16] = 0
    
    Sigma_mn = solve(Theta_mn)
    glasso_mn = cal_mean_PN(Theta_mn,Sigma_mn,1,n,p)
    
    dup_tune[dup] = glasso_mn$tune
    dup_FN[dup] = glasso_mn$FN
    dup_FP[dup] = glasso_mn$FP
  }# end for dup
  
  tune_para[4,p_iter] = mean(dup_tune)
  FN_matrix[4,p_iter] = mean(dup_FN)
  FP_matrix[4,p_iter] = mean(dup_FP)
  
  cat("Finish m-nearest network setting, p = ",p,"\n")
  
  ## Scale-free network setting
  # 20 % edges are connected
  
  dup_FN = rep(0,duplicate)
  dup_FP = rep(0,duplicate)
  dup_tune = rep(0,duplicate)
  
  for (dup in 1:duplicate) {
    Theta_rand = matrix(0,nrow = p, ncol = p)
    
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        if(runif(1)<0.2)
          Theta_rand[i,j] = sign(runif(1,-1,1))*runif(1,0.5,1)
      }# end for j
    }# end for i
    
    Theta_rand[lower.tri(Theta_rand)] = t(Theta_rand)[lower.tri(Theta_rand)]
    Theta_rand = Theta_rand + diag( colSums(abs(Theta_rand)) + rep(1,p) ) 
    Theta_rand[abs(Theta_rand) <= 1e-16] = 0
    
    Sigma_rand = solve(Theta_rand)
    glasso_rand = cal_mean_PN(Theta_rand,Sigma_rand,1,n,p)
    
    dup_tune[dup] = glasso_rand$tune
    dup_FN[dup] = glasso_rand$FN
    dup_FP[dup] = glasso_rand$FP
  
  }#end for dup
  
  tune_para[5,p_iter] = mean(dup_tune)
  FN_matrix[5,p_iter] = mean(dup_FN)
  FP_matrix[5,p_iter] = mean(dup_FP)
  
  cat("Finish scale-free network setting, p = ",p,"\n")
  
}# end for p_iter

#tune_para
#FN_matrix
#FP_matrix

###### ggplot ########

glasso_df = data.frame(FN = as.vector(FN_matrix),
                       FP = as.vector(FP_matrix),
                       p = rep(seq(20,120,20),each = 5), 
                       setting = rep(c("Sparse","Dense","Chain","M-nearest","Random"),times = 6))
glasso_df$setting = as.factor(glasso_df$setting)

levels(glasso_df$setting) = c("Sparse","Dense","Chain","M-nearest","Random")
new_color = c("#293757","#568D4B","#D5BB56","#D26A1B","#A41D1A")
simplot = ggplot(glasso_df,aes(x = p, y = FN, fill = setting)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=new_color) +
  labs(title = "False Negative rate vs p, under different settings") 
simplot 
