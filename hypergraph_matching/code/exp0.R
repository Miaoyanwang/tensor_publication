######## Experiment 0: acc vs sigma with unseeded and seeded algorithm ########

### Mar 8 22, Jiaxin Hu ###

source("functions.R")

library(ggplot2)
library(patchwork)


### start the simulation -----------------------------

d_vec = c(50, 80)
sigma_vec = c(0, 0.1, 0.25, 0.4)
rho_vec = sqrt(1 - sigma_vec^2)

dup = 5

# consider method: 
# 1. unseed; 2. seed without clean up; 3. seed with clean up
err_collect = array(0, dim = c(length(d_vec), length(sigma_vec), 3)) # d, sigma, method 
err_collect_sd = array(0, dim = c(length(d_vec), length(sigma_vec), 3))

for (i in 1:length(d_vec)) { # dim
  for (j in 1:length(sigma_vec)) { # sigma
    
    #i = 2; j = 3
    
    d = d_vec[i]
    rho = sqrt(1 - sigma_vec[j]^2)
    
    xi = 5*sqrt((log(d))^(1/2))
    zeta = 3*sqrt(sigma_vec[j]/d^2)
    
    err_vec1 = rep(NA,dup) # unseed
    err_vec2 = rep(NA,dup) # seed without clean up
    err_vec3 = rep(NA,dup) # seed with clean up
    
    for (n in 1:dup) {
      
      #n = 1
      
      cat("d = ", d ,"rho = ",rho," xi = ", xi," zeta = ", zeta , " dup = ", n, "\n")
      
      seed = i*1000 + j*100 + n
      dat = sim_corr_tensor(seed = seed, d = d, rho = rho)
      
      # unseed
      res_un = tensor_matching_unseed(dat$A, dat$B, p = 1)
      err_vec1[n] = match_error(res_un, dat$pi)
      
      # seed without clean up
      res_seed = tensor_matching_seed(dat$A, dat$B, p = 1, xi = xi, zeta = zeta, clean_up = F)
      if(! is.null(res_seed )){
        err_vec2[n] = match_error(res_seed, dat$pi)
      }
      
      
      # seed with clean up
      res_seed_c = tensor_matching_seed(dat$A, dat$B, p = 1, xi = xi, zeta = zeta, clean_up = T)
      if(! is.null(res_seed_c)){
        err_vec3[n] = match_error(res_seed_c, dat$pi)
      }
      
    } # n
    
    err_collect[i, j, 1] = mean(err_vec1, na.rm = T)
    err_collect[i, j, 2] = mean(err_vec2, na.rm = T)
    err_collect[i, j, 3] = mean(err_vec3, na.rm = T)
    
    err_collect_sd[i, j, 1] = sd(err_vec1, na.rm = T)
    err_collect_sd[i, j, 2] = sd(err_vec2, na.rm = T)
    err_collect_sd[i, j, 3] = sd(err_vec3, na.rm = T)
    
  } # j
} # i 


err_collect

save(err_collect,err_collect_sd, file = "results/exp0.RData")


### visualization ---------------------

load("results/exp0.RData")

d_vec = c(50, 80)
sigma_vec = c(0, 0.1, 0.25, 0.4)

data = data.frame( err = c(err_collect[1,,1], err_collect[1,,2], err_collect[1,,3],
                           err_collect[2,,1], err_collect[2,,2], err_collect[2,,3]),
                   method = rep(c("unseed", "seed", "seed+clean"), each = length(sigma_vec), times = length(d_vec)),
                   d = rep(d_vec, each = 3*length(sigma_vec)),
                   sigma = rep(sigma_vec, times = 3*length(d_vec))
                   )
data$method = factor(data$method, levels = c("unseed", "seed", "seed+clean"))
data$d = as.character(data$d)
#data$sigma = as.factor(data$sigma)

g0 = ggplot(data, aes(x = sigma, y = err,  linetype=d))+
  geom_line(aes(color = method))+
  geom_point(aes(shape=method))+
  scale_color_manual(values=c("#999999", "#E69F00","#56B4E9"))+
  scale_linetype_manual(values=c("solid","dashed" )) +
  labs(title = "err vs sigma") +
  coord_cartesian(ylim = c(0, 1))+
  theme_light()
g0

pdf("figures/exp0.pdf", height = 6, width = 8)
g0
dev.off()
