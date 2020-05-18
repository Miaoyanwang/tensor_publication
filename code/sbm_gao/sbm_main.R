#### Algorithm for SBM in Gao ####

library(plot.matrix)

source("/Users/March/Desktop/myresearch/code/sbm_gao/bricks_sbm.R")
source("/Users/March/Desktop/myresearch/code/sbm_gao/get_sbm_data.R")
source("/Users/March/Desktop/myresearch/code/sbm_gao/greedy_initial_clustering.R")
source("/Users/March/Desktop/myresearch/code/sbm_gao/refinement_sbm.R")


## generate data
set.seed(1001)
test_data = get_sbm_data(50,3,20,10)
plot(test_data$x)
plot(test_data$truthX)
plot(test_data$prob)
test_data$truthC


## initialization clustering \sigma^0
test_ini = greedy_initial_clustering(test_data$x,3,mu = 0.5)
test_ini

svd(test_data$truthX)$d
# when lambda k of truth x is larger, the test ini is more accurate.

## refinement
test_refine = refinement_sbm(test_data$x,3,mu = 0.5)
test_refine$sigma_est
test_data$truthC

est_x = matrix(0,nrow = 50, ncol = 50)
for (i in 1:3) {
  for (j in 1:3) {
    est_x[test_refine$sigma_est == i,test_refine$sigma_est == j] = test_refine$B_est[i,j]
  }
}
plot(est_x)
plot(test_refine$B_est)

# plot x
axis = 1:50
data1 = expand.grid(X=axis, Y=axis)
gridx = test_data$x
dim(gridx) = c(2500,1)
data1$gridx = gridx

gridtruth = test_data$truthX
dim(gridtruth) = c(2500,1)
data1$gridtruth = gridtruth

gridjudge = est_x
dim(gridjudge) = c(2500,1)
data1$gridjudge = gridjudge

sbmx = ggplot(data1, aes(X, Y, fill= gridx)) + 
  geom_tile() +
  guides( fill = guide_legend(title = "x")) 

sbmt = ggplot(data1, aes(X, Y, fill= gridtruth)) + 
  geom_tile() +
  guides( fill = guide_legend(title = "truth")) 

sbme = ggplot(data1, aes(X, Y, fill= gridjudge)) + 
  geom_tile()+
  guides( fill = guide_legend(title = "judge")) 

sbm_show = plot_grid(sbmx,sbmt,sbme,labels = c('x', 'truth',"estimate"),nrow = 1,ncol = 3)
#sbm_show
# cal loss
gene_conf_mat_one(3,test_data$truthC,test_refine$sigma_est)
cal_loss_sbm(3,test_data$truthC,test_refine$sigma_est)

## simulation
# need to fix when cluster size is 1

# k = 3
d = seq(40,60,10)
mcr_seq_1 = rep(0,length(d))

for (iter in 1:length(d)) {
  
  mcr_dup = rep(0,10)
  for (dup in 1:10) {
    print(dup)
    sbm_data = get_sbm_data(d[iter],3,d[iter]*3/4,d[iter]*1/3)
    sbm_refine = refinement_sbm(sbm_data$x,3,mu = 0.5)
    mcr_dup[dup] = cal_loss_sbm(3,sbm_data$truthC,sbm_refine)
  }
  
  mcr_seq_1[iter] = mean(mcr_dup)
  
}


d = seq(70,90,10)
mcr_seq_2 = rep(0,length(d))

for (iter in 1:length(d)) {
  
  mcr_dup = rep(0,10)
  for (dup in 1:10) {
    print(dup)
    sbm_data = get_sbm_data(d[iter],3,d[iter]*3/4,d[iter]*1/3)
    sbm_refine = refinement_sbm(sbm_data$x,3,mu = 0.5)
    mcr_dup[dup] = cal_loss_sbm(3,sbm_data$truthC,sbm_refine)
  }
  
  mcr_seq_2[iter] = mean(mcr_dup)
  
}

mcr_seq = c(mcr_seq_1,mcr_seq_2)
dchange = seq(40,90,10)
sbm_sim_result = data.frame(mcr = mcr_seq, n = dchange)
sbm_sim_result$expn = 10^20*exp(-sbm_sim_result$n)
library(ggplot2)
library(cowplot)

sbm_plot = ggplot(data = sbm_sim_result,aes(x = n,y = mcr_seq)) +
  geom_line() +
  labs(x = "n", y = "mcr_sbm",title = "MCR vs Dimension")
sbm_plot

sbm_plot2 = ggplot(data = sbm_sim_result,aes(x = expn ,y = mcr_seq)) +
  geom_point() +
  labs(x = "exp(-n)", y = "mcr_sbm",title = "MCR vs exp(-Dimension)") 
sbm_plot2

plot_grid(sbm_plot,sbm_plot2,nrow = 1,ncol = 2)

