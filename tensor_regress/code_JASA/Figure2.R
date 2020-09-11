#### code for Figure 2. Log-likelihood trajectory for three models under different (d,r) combinations
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
library(tensorregress)
library(ggplot2)
library(patchwork)


set.seed(0)## set seed for the entire simulation
model = c("normal","binary","poisson")
dimen = c(30,25)
rank = c(6,3)

################### simulation ##################

lglk_sum = list()
true_lglk_sum = list()
time_sum = list()


for (i in 1:3) {## three data types
  lglk = list()
  true_lglk = c()
  time = c()
  k = 1
  for (j in 1:2) {## two possbiel dimensions
      for (r in 1:2) {## two possible ranks
        if(i==1) signal=3
        else if(i==2) signal=10
        else if(i==3) signal=6
      my_data = sim_data(seed = NA, rep(dimen[j],3), rep(rank[r],3), p = rep(floor(0.4*dimen[j]) ,3),dist = model[i], dup = 1,signal = signal)
      ptm = proc.time()
      my_reg =  tensor_regress(tsr = my_data$tsr[[1]], X_covar1 = my_data$X_covar1, X_covar2 = my_data$X_covar2,X_covar3 = my_data$X_covar3,rep(rank[r],3), Nsim = 20, cons = "non",dist = model[i],traj_long=TRUE)
      time1 = proc.time() - ptm
      time[k] = time1[3] ## running time
      true_lglk[k] = loglike(my_data$tsr[[1]], my_data$U,dist = model[i])
      lglk[[k]] = my_reg$lglk
      k = k+1
    }
  }
  lglk_sum[[i]] = lglk
  true_lglk_sum[[i]] = true_lglk
  time_sum[[i]] = time
  
}

############ ggplot ##################

#### normal model with d = 30, r = 6,3
lglk_normal = c(lglk_sum[[1]][[1]] ,lglk_sum[[1]][[2]])
case_normal=  c(rep("d3r6",33), rep("d3r3",33))
iter = rep(1:33,times = 2)
olglk1 = true_lglk_sum[[1]][1]
olglk2 = true_lglk_sum[[1]][2]
plotdata = data.frame(lglk = lglk_normal,case = case_normal,iter = iter)

p1 = ggplot(data = plotdata, aes(x = iter,y = lglk/10000)) + 
  geom_line(aes(color = case),lwd = 1) +
  geom_point(aes(shape = case),size = 2) +
  scale_color_manual(values=c("#6F9B3C","#F2AB1D"),labels = c("d = 30, r = 3 (0.7 sec/iter)", "d = 30, r = 6 (4.1 sec/iter)")) + 
  scale_shape_manual(values = c(20,17),labels = c("d = 30, r = 3 (0.7 sec/iter)","d = 30, r = 6 (4.1 sec/iter)")) +
  geom_hline(yintercept=olglk1/10000, color = "#F2AB1D" , lty = "dashed",lwd = 1)+
  geom_hline(yintercept=olglk2/10000, color = "#6F9B3C",lty = "dashed",lwd = 1)+
  guides( color = guide_legend(title = "Setting"), linetype = guide_legend("Setting"), shape = guide_legend("Setting"))+
  labs(y = "log-likelihood x 10^4",x = "iteration",size = 16) +
  theme(axis.text.y = element_text(size = 16 ) , 
        axis.text.x = element_text(size = 16 ),
        legend.text=element_text(size=16), 
        legend.title = element_text(size=16),
        axis.title=element_text(size=16),
        legend.position = "none")
p1

#### normal model with d = 25, r = 6,3
lglk_normal = c(lglk_sum[[1]][[3]] ,lglk_sum[[1]][[4]])
case_normal=  c(rep("d2r6",33), rep("d2r3",33))
iter = rep(1:33,times = 2)
olglk1 = true_lglk_sum[[1]][3]
olglk2 = true_lglk_sum[[1]][4]
plotdata = data.frame(lglk = lglk_normal,case = case_normal,iter = iter)

p2 = ggplot(data = plotdata, aes(x = iter,y = lglk/10000)) + 
  geom_line(aes(color = case),lwd = 1) +
  geom_point(aes(shape = case),size = 2) +
  scale_color_manual(values=c("#43ABE4","#C76AB2"),labels = c("d = 25, r = 3 (0.7 sec/iter)", "d = 25, r = 6 (1.8 sec/iter)")) + 
  scale_shape_manual(values = c(18,4),labels = c("d = 25, r = 3 (0.7 sec/iter)", "d = 25, r = 6 (1.8 sec/iter)")) +
  geom_hline(yintercept=olglk1/10000, color = "#C76AB2" , lty = "dashed",lwd = 1)+
  geom_hline(yintercept=olglk2/10000, color = "#43ABE4",lty = "dashed",lwd = 1)+
  guides( color = guide_legend(title = "Setting"), shape = guide_legend("Setting"))+
  labs(title = "a. Normal" ,y = "log-likelihood x 10^4",size = 16) +
  theme(axis.text.y = element_text(size = 16 ) , 
        axis.text.x =  element_text(size = 16 ),
        legend.text =element_text(size=16), 
        legend.title = element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x = element_blank(),
        plot.title = element_text( size = 20, face="bold"),
        legend.position = "none")
p2

#### binary model with d = 30, r = 6,3
lglk_bin = c(lglk_sum[[2]][[1]] ,lglk_sum[[2]][[2]])
case_bin=  c(rep("d3r6",33), rep("d3r3",33))
iter = rep(1:33,times = 2)
olglk1 = true_lglk_sum[[2]][1]
olglk2 = true_lglk_sum[[2]][2]
plotdata = data.frame(lglk = lglk_bin,case = case_bin,iter = iter)

p3 = ggplot(data = plotdata, aes(x = iter,y = lglk/10000)) + 
  geom_line(aes(color = case),lwd = 1) +
  geom_point(aes(shape = case),size = 2) +
  scale_color_manual(values=c("#6F9B3C","#F2AB1D"),labels = c("d = 30, r = 3 (0.7 sec/iter)", "d = 30, r = 6 (4.1 sec/iter)")) + 
  scale_shape_manual(values = c(20,17),labels = c("d = 30, r = 3 (0.7 sec/iter)","d = 30, r = 6 (4.1 sec/iter)")) +
  geom_hline(yintercept=olglk1/10000, color = "#F2AB1D" , lty = "dashed",lwd = 1)+
  geom_hline(yintercept=olglk2/10000, color = "#6F9B3C",lty = "dashed",lwd = 1)+
  guides( color = guide_legend(title = "Setting"), linetype = guide_legend("Setting"), shape = guide_legend("Setting"))+
  labs(x = "iteration",size = 16) +
  theme(axis.text.y = element_text(size = 16 ) , 
        axis.text.x = element_text(size = 16 ),
        legend.text=element_text(size=16), 
        legend.title = element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x =element_text(size=16),
        legend.position = "none")

p3

#### binary model with d = 25, r = 6,3
lglk_bin = c(lglk_sum[[2]][[3]] ,lglk_sum[[2]][[4]])
case_bin=  c(rep("d2r6",33), rep("d2r3",33))
iter = rep(1:33,times = 2)
olglk1 = true_lglk_sum[[2]][3]
olglk2 = true_lglk_sum[[2]][4]
plotdata = data.frame(lglk = lglk_bin,case = case_bin,iter = iter)


p4 = ggplot(data = plotdata, aes(x = iter,y = lglk/10000)) + 
  geom_line(aes(color = case),lwd = 1) +
  geom_point(aes(shape = case),size = 2) +
  scale_color_manual(values=c("#43ABE4","#C76AB2"),labels = c("d = 25, r = 3 (0.7 sec/iter)", "d = 25, r = 6 (1.8 sec/iter)")) + 
  scale_shape_manual(values = c(18,4),labels = c("d = 25, r = 3 (0.7 sec/iter)", "d = 25, r = 6 (1.8 sec/iter)")) +
  geom_hline(yintercept=olglk1/10000, color = "#C76AB2" , lty = "dashed",lwd = 1)+
  geom_hline(yintercept=olglk2/10000, color = "#43ABE4",lty = "dashed",lwd = 1)+
  guides( color = guide_legend(title = "Setting"), shape = guide_legend("Setting"))+
  labs(title = "b. Binary" ,size = 16) +
  theme(axis.text.y = element_text(size = 16 ) , 
        axis.text.x =  element_text(size = 16 ),
        legend.text =element_text(size=16), 
        legend.title = element_text(size=16),
        axis.title = element_blank(),
        plot.title = element_text( size = 20, face="bold"),
        legend.position = "none")

p4

#### poisson model with d = 30, r = 6,3
lglk_pos = c(lglk_sum[[3]][[1]] ,lglk_sum[[3]][[2]])
case_pos=  c(rep("d3r6",33), rep("d3r3",33))
iter = rep(1:33,times = 2)
olglk1 = true_lglk_sum[[3]][1]
olglk2 = true_lglk_sum[[3]][2]
plotdata = data.frame(lglk = lglk_pos,case = case_pos,iter = iter)

p5 = ggplot(data = plotdata, aes(x = iter,y = lglk/10000)) + 
  geom_line(aes(color = case),lwd = 1) +
  geom_point(aes(shape = case),size = 2) +
  scale_color_manual(values=c("#6F9B3C","#F2AB1D"),labels = c("d = 30, r = 3", "d = 30, r = 6")) + 
  scale_shape_manual(values = c(20,17),labels = c("d = 30, r = 3","d = 30, r = 6")) +
  geom_hline(yintercept=olglk1/10000, color = "#F2AB1D" , lty = "dashed",lwd = 1)+
  geom_hline(yintercept=olglk2/10000, color = "#6F9B3C",lty = "dashed",lwd = 1)+
  guides( color = guide_legend(title = "Setting"), linetype = guide_legend("Setting"), shape = guide_legend("Setting"))+
  labs(x = "iteration",size = 16) +
  theme(axis.text.y = element_text(size = 16 ) , 
        axis.text.x = element_text(size = 16 ),
        legend.text=element_text(size=16), 
        legend.title = element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16)
  )

p5


###### poisson model with d = 25, r = 6,3
lglk_pos = c(lglk_sum[[3]][[3]] ,lglk_sum[[3]][[4]])
case_pos=  c(rep("d2r6",33), rep("d2r3",33))
iter = rep(1:33,times = 2)
olglk1 = true_lglk_sum[[3]][3]
olglk2 = true_lglk_sum[[3]][4]
plotdata = data.frame(lglk = lglk_pos,case = case_pos,iter = iter)

p6 = ggplot(data = plotdata, aes(x = iter,y = lglk/10000)) + 
  geom_line(aes(color = case),lwd = 1) +
  geom_point(aes(shape = case),size = 2) +
  scale_color_manual(values=c("#43ABE4","#C76AB2"),labels = c("d = 25, r = 3", "d = 25, r = 6")) + 
  scale_shape_manual(values = c(18,4),labels = c("d = 25, r = 3", "d = 25, r = 6")) +
  geom_hline(yintercept=olglk1/10000, color = "#C76AB2" , lty = "dashed",lwd = 1)+
  geom_hline(yintercept=olglk2/10000, color = "#43ABE4",lty = "dashed",lwd = 1)+
  guides( color = guide_legend(title = "Setting"), shape = guide_legend("Setting"))+
  labs(title = "c. Poisson" ,size = 16) +
  theme(axis.text.y = element_text(size = 16 ) , 
        axis.text.x =  element_text(size = 16 ),
        legend.text =element_text(size=16), 
        legend.title = element_text(size=16),
        axis.title = element_blank(),
        plot.title = element_text( size = 20, face="bold"))
p6

##### plot all figures
(p2|p4|p6)/
(p1|p3|p5)


