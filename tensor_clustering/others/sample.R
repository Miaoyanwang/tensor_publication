####### comparison #########

library(RSKC)
source("functions.R")
source("HOLloyd.R")
source("tensor_score.R")

set.seed(1)

############## CER vs gamma
p <- 100
r <- 5

#gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)
##gamma_range <- c(-2.1, -1.9, -1.7, -1.5, -1.4)
gamma_range <- c(-1.5,-1.4,-1.3,-1.2,-1.1)


alpha_range = 1 + 100*p^(gamma_range/2)

s_min <- 0.5
delta = 0.5

dist <- "normal"
theta_dist <- "abs_normal"


mis <- array(NA, dim = c(10, length(gamma_range)))
mis_sd <- array(NA, dim = c(10, length(gamma_range)))
dup <- 5


for (g in 1:length(gamma_range)) {
  
  #g = 1; n = 1
  s_max <- s_min * alpha_range[g]
  # our_err <- c()
  our_err2 <- c()
  #han_err <- c()
  #hosvd_err1 <- c()
  #hosvd_err2 <- c()
  ke_err <- c()
  # ke_err2 <- c()
  # ke_err3 <- c()
  # ke_err4 <- c()
  # ke_err5 <- c()
  
  for (n in 1:dup) {
    
    cat("gamma = ", gamma_range[g], "s_max = ", s_max, "dup = ", n, "\n")
    
    seed = 100*g + n
    set.seed(seed)
    
    dat <- sim_hDCBM_network(seed = seed, p, r, core_conrtol = "signal_a1",
                             delta = delta, s_min = s_min, s_max = s_max,
                             dist =  dist, sigma = 0.75, theta_dist = theta_dist)
    
    # # ours
    # ini_result <- wkmedian(dat$Y, r)
    # re_result <- hALloyd(dat$Y, ini_result$z0, max_iter = 50)
    # our_err <- c(our_err, CER(re_result, dat$z))
    
    # ours2
    ini_result <- wkmeans(dat$Y, r)
    re_result <- hALloyd(dat$Y, ini_result$z0, max_iter = 20)
    our_err2 <- c(our_err2, CER(re_result, dat$z))
    
    # # han
    #  z.HOSC = HO.SC(as.tensor(dat$Y), r) # high-order initialization
    # z.Lloyd.HOSC = HO.Lloyd(as.tensor(dat$Y), z.HOSC) # HLloyd iteration
    # han_err = c(han_err,min(CER(z.Lloyd.HOSC[[1]], dat$z) ,CER(z.Lloyd.HOSC[[2]], dat$z) ,CER(z.Lloyd.HOSC[[3]], dat$z) ))
    # 
    # # hosvd1
    #hosvd1 = HOSVD_kmeans(dat$Y, r)
    # hosvd_err1 = c(hosvd_err1, CER(hosvd1, dat$z))
    # 
    # # hosvd2
    #  hosvd2 = HOSVD_score(dat$Y, r)
    #hosvd_err2 = c(hosvd_err2, CER(hosvd2, dat$z))
    
    # tensor score
    ke_result =  tensor_score_adj(dat$Y, r, rm_diag = T, hooi = T, 
                                  reg_hooi = T, score_reg = T, max_iter =  20)
    ke_err = c(ke_err, CER(ke_result,dat$z))
    
    # # no thresholding and diag rm + regularized hooi
    # ke_result2 = tensor_score_adj(dat$Y, r, rm_diag = F, hooi = T,
    #                               reg_hooi = T, score_reg = F, max_iter =  50)
    # ke_err2 = c(ke_err2, CER(ke_result2,dat$z))
    # 
    # # no thresholding and diag rm + no regularized hooi
    # ke_result3 = tensor_score_adj(dat$Y, r, rm_diag = F, hooi = T,
    #                               reg_hooi = F, score_reg = F, max_iter =  50)
    # ke_err3 = c(ke_err3, CER(ke_result3,dat$z))
    # 
    # # no thresholding and diag rm + no hooi
    # ke_result4 = tensor_score_adj(dat$Y, r, rm_diag = F, hooi = F,
    #                               reg_hooi = F, score_reg = F, max_iter =  50)
    # ke_err4 = c(ke_err4, CER(ke_result4,dat$z))
    # 
    # # no hooi
    # ke_result5 = tensor_score_adj(dat$Y, r, rm_diag = T, hooi = F,
    #                               reg_hooi = F, score_reg = T, max_iter =  50)
    # ke_err5 = c(ke_err5, CER(ke_result5,dat$z))
    
  }
  
  #mis[1, g] <- mean(our_err)
  mis[2, g] <- mean(our_err2)
  # mis[3, g] <- mean(han_err)
  #mis[4, g] <- mean(hosvd_err1)
  #mis[5, g] <- mean(hosvd_err2)
  mis[6, g] <- mean(ke_err)
  # mis[7, g] <- mean(ke_err2)
  # mis[8, g] <- mean(ke_err3)
  # mis[9, g] <- mean(ke_err4)
  #mis[10, g] <- mean(ke_err5)
}
library(ggplot2)
data=mis[2:6,]
rownames(data)=c("ours","han","hosvd","hosvd+","Ke")

require(reshape2)
dataL=melt(data)
dataL[,2]=rep(gamma_range,rep(5,5))
p  <- ggplot(dataL, aes_string(x="Var2", y="value", colour="Var1", group="Var1"))
p <- p + geom_line()+geom_point()

plot(gamma_range,mis[2,],type="l",ylim=c(0,0.25))
points(gamma_range,mis[3,],type="l",col="red")
points(gamma_range,mis[4,],type="l",col="blue")
points(gamma_range,mis[5,],type="l",col="green")
points(gamma_range,mis[6,],type="l",col="purple")
dev.copy(pdf,"compare.pdf")




############## CER vs a (shape parameter)

p <- 100
r <- 5

shape_range = seq(5.5,1.5,-0.5)

s_min <- 0.5
s_max <- 1
delta = 0.25

#dist <- "normal"
dist <- "binary"
theta_dist <- "pareto"


mis <- array(NA, dim = c(10, length(shape_range)))
mis_sd <- array(NA, dim = c(10, length(shape_range)))
dup <- 5


for (g in 1:length(gamma_range)) {
  #g = 1; n = 1
  
  shape_para = shape_range[g]
  scale_para = (shape_para-1)/shape_para
  
  #our_err <- c()
  our_err2 <- c()
  # han_err <- c()
  # hosvd_err1 <- c()
  # hosvd_err2 <- c()
  ke_err <- c()
  # ke_err2 <- c()
  # ke_err3 <- c()
  # ke_err4 <- c()
  # ke_err5 <- c()
  
  for (n in 1:dup) {
    
    cat("shape = ", shape_para, "scale = ", scale_para, "dup = ", n, "\n")
    
    seed = 100*g + n
    set.seed(seed)
    
    dat = sim_hDCBM_network(seed = seed, p, r, core_conrtol = "signal_a1",
                            delta = delta, s_min = s_min, s_max = s_max,
                            dist =  dist, sigma = 1, theta_dist = theta_dist, 
                            alpha = shape_para, beta = scale_para)
    
    # # ours
    # ini_result <- wkmedian(dat$Y, r)
    # re_result <- hALloyd(dat$Y, ini_result$z0, max_iter = 50)
    # our_err <- c(our_err, CER(re_result, dat$z))
    
    # ours2
    ini_result <- wkmeans(dat$Y, r)
    re_result <- hALloyd(dat$Y, ini_result$z0, max_iter = 20)
    our_err2 <- c(our_err2, CER(re_result, dat$z))
    
    # # han
    # z.HOSC = HO.SC(as.tensor(dat$Y), r) # high-order initialization
    # z.Lloyd.HOSC = HO.Lloyd(as.tensor(dat$Y), z.HOSC) # HLloyd iteration
    # han_err = c(han_err,min(CER(z.Lloyd.HOSC[[1]], dat$z) ,CER(z.Lloyd.HOSC[[2]], dat$z) ,CER(z.Lloyd.HOSC[[3]], dat$z) ))
    # 
    # # hosvd1
    # hosvd1 = HOSVD_kmeans(dat$Y, r)
    # hosvd_err1 = c(hosvd_err1, CER(hosvd1, dat$z))
    # 
    # # hosvd2
    # hosvd2 = HOSVD_score(dat$Y, r)
    # hosvd_err2 = c(hosvd_err2, CER(hosvd2, dat$z))
    
    # tensor score
    ke_result =  tensor_score_adj(dat$Y, r, rm_diag = T, hooi = T, 
                                  reg_hooi = T, score_reg = T, max_iter =  20)
    ke_err = c(ke_err, CER(ke_result,dat$z))
    
    # # no thresholding and diag rm + regularized hooi
    # ke_result2 = tensor_score_adj(dat$Y, r, rm_diag = F, hooi = T,
    #                               reg_hooi = T, score_reg = F, max_iter =  50)
    # ke_err2 = c(ke_err2, CER(ke_result2,dat$z))
    # 
    # # no thresholding and diag rm + no regularized hooi
    # ke_result3 = tensor_score_adj(dat$Y, r, rm_diag = F, hooi = T,
    #                               reg_hooi = F, score_reg = F, max_iter =  50)
    # ke_err3 = c(ke_err3, CER(ke_result3,dat$z))
    # 
    # # no thresholding and diag rm + no hooi
    # ke_result4 = tensor_score_adj(dat$Y, r, rm_diag = F, hooi = F,
    #                               reg_hooi = F, score_reg = F, max_iter =  50)
    # ke_err4 = c(ke_err4, CER(ke_result4,dat$z))
    # 
    # # no hooi
    # ke_result5 = tensor_score_adj(dat$Y, r, rm_diag = T, hooi = F,
    #                               reg_hooi = F, score_reg = T, max_iter =  50)
    # ke_err5 = c(ke_err5, CER(ke_result5,dat$z))
    
  }
  
  #mis[1, g] <- mean(our_err)
  mis[2, g] <- mean(our_err2)
  # mis[3, g] <- mean(han_err)
  # mis[4, g] <- mean(hosvd_err1)
  # mis[5, g] <- mean(hosvd_err2)
  mis[6, g] <- mean(ke_err)
  # mis[7, g] <- mean(ke_err2)
  # mis[8, g] <- mean(ke_err3)
  # mis[9, g] <- mean(ke_err4)
  #mis[10, g] <- mean(ke_err5)
}

mis

