############### Experiment 1: gap between stats and comp limits ###########
library(RSKC)
library(LICORS)
source("functions.R")
source("HOLloyd.R")


########## simulation #########
p <- 100
r <- 5

delta_lower = 0.7
delta_upper = 1.8

delta.candidate = seq(delta_lower,delta_upper,0.1)

dist <- "normal"
theta_dist <- "non"
mis <- array(NA, dim = c(2, length(delta.candidate)))
dup <- 5

# gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5)
# alpha_range =  1/(1 - 1700*p^gamma_range) 

#alpha_range <- 38000 * p^gamma_range
# s_min <- 0.5

for (g in 1:length(delta.candidate)) {
  
  delta =  200*p^(-delta.candidate[g])
  comp_err <- c()
  stat_err <- c()
  Lloyd_err <- c()
  oracle_error <- c()
  for (n in 1:dup) {
    #g = 1; n = 1
    cat("delta = ", delta, "dup = ", n, "\n")
    
    dat = sim_hDCBM_network(seed = NA, p, r, core_conrtol = "signal_m",delta = delta, dist =  dist, sigma = 1, theta_dist = theta_dist)
    
    # comp
    ini_result <- wkmedian(dat$Y, r)
    re_result <- hALloyd(dat$Y, ini_result$z0, max_iter = 20)
    comp_err <- c(comp_err, CER(re_result, dat$z))
    
    # stats
    or_result <- hALloyd(dat$Y, dat$z, max_iter = 20)
    stat_err <- c(stat_err, CER(or_result, dat$z))
    
    # han HSC + HLloyd
    HSC = HO.SC(as.tensor(dat$Y), rep(r,3))
    Lloyd = HO.Lloyd(as.tensor(dat$Y), HSC)
    
    Lloyd_err = c(Lloyd_err, CER(Lloyd[[1]],dat$z))
    
    # oracle
    
    oracle = HO.Lloyd(as.tensor(dat$Y), list(dat$z,dat$z,dat$z))
    oracle_error = c(oracle_error, CER(oracle[[1]],dat$z))
  }
  
  mis[1, g] <- mean(comp_err)
  mis[2, g] <- mean(stat_err)
  mis[3, g] <- mean(Lloyd_err)
  mis[4, g] <- mean(oracle_error)
}


mis

plot(delta.candidate, mis[1, ] * 2, type = "b", xlab = "-gamma", ylab = "MCR")
lines(delta.candidate, mis[2, ] * 2, type = "b", col = 2)
lines(delta.candidate, mis[3, ] * 2, type = "b", col = 1, lty = 3)
lines(delta.candidate, mis[4, ] * 2, type = "b", col = 2, lty = 3)
legend("topleft", c("wkmedian+hALloyd", "true+hALloyd", "HSC+hLloyd", "true+hLloyd"), 
       col = c(1,2, 1,2), lty = c(1,1,3,3))

