###### test simulations for Gaussian tensor matching #####
##### Feb 28 22 ########

source("functions.R")

library(ggplot2)
library(patchwork)

### show seeded algorithm works well #####

d = 50
rho_vec = c(0.2,0.5,0.8)
seed_n_vec = c(2,4,6)
dup = 5

err_collect = array(0, dim = c(length(rho_vec), length(seed_n_vec)))
err_collect_sd = array(0, dim = c(length(rho_vec), length(seed_n_vec)))

for (r in 1:length(rho_vec)) { # 
  for (s in 1:length(seed_n_vec)) {
    
    #r = 2; s = 3
    
    rho = rho_vec[r]
    seed_n = seed_n_vec[s]
    
    err_vec = rep(0,dup)
    
    for (n in 1:dup) {
      #n = 3
    
      cat("rho = ",rho," seed_n = ", seed_n, " dup = ", n, "\n")
      
      dat = sim_corr_tensor(d = d, rho = rho)
      
      seed_set = dat$pi[sample(1:d, seed_n),]
      mt = seeded_bipartite_matching(dat$A,dat$B,seed_set)
      
      err_vec[n] = match_error(mt,dat$pi)
    }# n
    
    err_collect[r, s] = mean(err_vec)
    err_collect_sd[r, s] = sd(err_vec)
  }# s
}# r


data = data.frame( err = c(err_collect[,1], err_collect[,2],err_collect[,3]),
                   seed_n = rep(seed_n_vec, each = 3),
                   rho = rep(rho_vec, times = 3))
data$seed_n = as.factor(data$seed_n)
data$rho = as.factor(data$rho)

g0 = ggplot(data, aes(x = seed_n, y = err , fill = rho))+
  geom_bar(position="dodge",stat="identity")+
  labs(title = "p = 50")+
  theme_light()

pdf("figures/err_vs_seedn.pdf", width = 6, height = 4)
g0
dev.off()

#### show emp method does not work 

d_vec = c(10,30,50)
rho_vec = c(0.8,0.9,0.95,1)
dup = 5

err_collect = array(0, dim = c(length(rho_vec), length(d_vec), 2)) # rho, d, p
err_collect_sd = array(0, dim = c(length(rho_vec), length(d_vec), 2))

for (r in 1:length(rho_vec)) {
  for (i in 1:length(d_vec)) {
    #r = 2; i = 2
    
    rho = rho_vec[r]
    d = d_vec[i]
    
    err_vec_p1 = rep(0,dup)
    err_vec_p2 = rep(0,dup)
    
    for (n in 1:dup) {
      
      cat(" rho = ", rho,"d = ",d, " dup = ", n, "\n")
      
      dat = sim_corr_tensor(d = d, rho = rho)
      
      mt1 = tensor_matching_unseed(dat$A, dat$B, p =1)
      mt2 = tensor_matching_unseed(dat$A, dat$B, p =2)
      
      err_vec_p1[n] = match_error(mt1,dat$pi)
      err_vec_p2[n] = match_error(mt1,dat$pi)
    }#n
    
    err_collect[r,i,1] = mean(err_vec_p1)
    err_collect[r,i,2] = mean(err_vec_p2)
    
    err_collect_sd[r,i,1] = sd(err_vec_p1)
    err_collect_sd[r,i,2] = sd(err_vec_p2)
  }# i
}#r

# err vs rho
p = 1
data = data.frame( err = c(err_collect[,1,p], err_collect[,2,p],err_collect[,3,p]),
                   d = rep(d_vec, each = 4),
                   rho = rep(rho_vec, times = 3))
data$rho = as.factor(data$rho)
data$d = as.factor(data$d)

g11 = ggplot(data, aes(x = rho, y = err , fill = d))+
  geom_bar(position="dodge",stat="identity")+
  labs(title = "unseeded, err vs rho, d1")+
  theme_light()
g11

p = 2
data = data.frame( err = c(err_collect[,1,p], err_collect[,2,p],err_collect[,3,p]),
                   d = rep(d_vec, each = 4),
                   rho = rep(rho_vec, times = 3))
data$rho = as.factor(data$rho)
data$d = as.factor(data$d)

g12 = ggplot(data, aes(x = rho, y = err , fill = d))+
  geom_bar(position="dodge",stat="identity")+
  labs(title = "unseeded, err vs rho, d2")+
  theme_light()
g12

pdf("figures/err_vs_rho_unseed.pdf", width = 10, height = 4)
(g11|g12)
dev.off()


# # specific example
# d = 50
# rho = 0.99
# dat = sim_corr_tensor(seed = 1, d = d, rho = rho)
# 
# t1 = tensor_matching_unseed(dat$A, dat$B, p =1)
# t1
# match_error(t1,dat$pi)
# 
# 
# t3 = tensor_matching_unseed(dat$A, dat$B, p =1, test = TRUE)
# #t3
# match_error(t3,dat$pi)
# 
# 
# 
# # check distance matrix
# A = dat$A
# B = dat$B
# d_collect = matrix(0, nrow = d, ncol = d)
# for (i in 1:d) {
#   for (k in 1:d) {
#     # i = 1; k = 15
#     d_collect[i,k] = dist_p(as.vector(A[i,,]), as.vector(B[k,,]),p)
#   }
# }
# 
# d_collect
# 
# # (8,1) true pair
# d_collect[8,1]
# 
# plot(sort(as.vector(A[8,,])), sort(as.vector(B[1,,])))
# 
# # (10,7) fake pair
# d_collect[10,7]
# 
# plot(sort(as.vector(A[10,,])), sort(as.vector(B[7,,])))
# 
# match_error(HungarianSolver(d_collect)$pairs,dat$pi) 
# 
# 
