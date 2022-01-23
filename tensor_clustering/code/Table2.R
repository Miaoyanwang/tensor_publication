# Code for Table 2. Verification BIC criterion for cluster number selection

# load dependencies ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)


# report simulation statistics from prior simulation -------------------------------------------------------------
# if a new run of simulation is desired, please go to line #58 to run the code ----------------------------------

load("presaved/Table2.RData")

# start a new run of simulation ------------------------------------------------------------------------------------------------

r_vec = c(2,4)
sigma_vec = c(0.5,1)
p_vec = c(50,80)

mis = mis_sd = array(0,dim = c(2,2,2))
dup = 30

s_min <- 0.5
s_max <- 20
delta = 0.5

for (i in 1:length(r_vec)) { # r
  for (j in 1:length(sigma_vec)) { # sigma
    for (k in 1:length(p_vec)) { # p
      
      r_res_vec = rep(0,dup)
      for (n in 1:dup) { # dup
        
        r = r_vec[i]
        sigma = sigma_vec[j]
        p = p_vec[k]
        
        seed = 10000*i+ 1000*j + 100*k + n
        set.seed(seed)
        
        dat = sim_dTBM(seed = seed, imat = F, asymm = F, p = rep(p,3), r = rep(r,3),
                       core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
                       dist = "normal", sigma = sigma,
                       theta_dist = "abs_normal")
        
        r_range = cbind(2:6, 2:6, 2:6)
        r_res_vec[n] = select_r(dat$Y, r_range = r_range, asymm = F)$r[1]
        
      } # end dup
      
      mis[i,j,k] = mean(r_res_vec)
      mis_sd[i,j,k] = sd(r_res_vec)
      
    } # end p
  } # end sigma
} # end r


save(mis, mis_sd, file = "presaved/Table2.RData")