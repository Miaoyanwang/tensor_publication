######### Coauthor tensor matching analysis

######### Oct 20 Jiaxin Hu

# Dependenices ----------------------------------------------

library(data.table) # fread
library(tidyr) # separate
library(gtools) # permutation
library(ggplot2)
#library(patchwork)
library(cowplot)
source("functions_mm.R") # methods


# Load data --------------------------------------------------

# # node - author name
# node_label = fread("coauth-DBLP/coauth-DBLP-node-labels.txt", sep = ",", header = F)
# colnames(node_label) = "node_label"
# node_label = separate(node_label, "node_label", into = c("node", "label"), sep = "^\\S*\\K\\s+")

# number of authors in the pub
nverts = fread("real_data/coauth-DBLP/coauth-DBLP-nverts.txt", header = F)
nverts = nverts$V1

# timestamps
timestamps = fread("real_data/coauth-DBLP/coauth-DBLP-times.txt", header = F)
timestamps = timestamps$V1

# simplex - authors in one pub
simplices = fread("real_data/coauth-DBLP/coauth-DBLP-simplices.txt", header = F)
simplices = simplices$V1
simplex_list = split(simplices, as.factor( rep(1:length(nverts), times = nverts) ) )
remove(simplices)

# Sub-sampling -------------------------------------------------

pub_year = c(2009,2010,2011,2012)

# find publications in each year with >= 3 authors
pub_ind1 = timestamps %in% pub_year[1]
pub_ind2 = timestamps %in% pub_year[2]
pub_ind3 = timestamps %in% pub_year[3]
pub_ind4 = timestamps %in% pub_year[4]

simplex1 = simplex_list[pub_ind1 & nverts >= 3]
simplex2 = simplex_list[pub_ind2 & nverts >= 3]

# select highly active authors with number of publications >=3 in first two year
freq_author1 = as.data.frame(table(unlist(simplex1)))
freq_author2 = as.data.frame(table(unlist(simplex2)))

author1 = unique(as.numeric(freq_author1[freq_author1$Freq >= 3, ]$Var1))
author2 = unique(as.numeric(freq_author2[freq_author2$Freq >= 3, ]$Var1))

active_author = Reduce(intersect, list(author1, author2))

# select authors that have pubs with other two active authors in the first two years
rel_pub1_ind = unlist(lapply(simplex1, function(x) sum(x %in% active_author) >= 3))
sum(rel_pub1_ind == 1)

rel_pub2_ind = unlist(lapply(simplex2, function(x) sum(x %in% active_author) >= 3))
sum(rel_pub2_ind == 1)

sel_author1 = intersect(unique(unlist(simplex1[rel_pub1_ind])), active_author) 
sel_author2 = intersect(unique(unlist(simplex2[rel_pub2_ind])), active_author)

# selected authors
sel_author = intersect(sel_author1, sel_author2)

# find all pubs related to the selected authors in each year
sel_pub1_ind = unlist(lapply(simplex_list[pub_ind1], function(x) sum(x %in% sel_author) >= 1))
sel_pub2_ind = unlist(lapply(simplex_list[pub_ind2], function(x) sum(x %in% sel_author) >= 1))
sel_pub3_ind = unlist(lapply(simplex_list[pub_ind3], function(x) sum(x %in% sel_author) >= 1))
sel_pub4_ind = unlist(lapply(simplex_list[pub_ind4], function(x) sum(x %in% sel_author) >= 1))

sel_simplex1 = simplex_list[pub_ind1][sel_pub1_ind]
sel_simplex2 = simplex_list[pub_ind2][sel_pub2_ind]
sel_simplex3 = simplex_list[pub_ind3][sel_pub3_ind]
sel_simplex4 = simplex_list[pub_ind4][sel_pub4_ind]

sel_simplex_list = list(sel_simplex1, sel_simplex2,sel_simplex3,sel_simplex4)

# Construct coauthorship tensors ----------------------------------------

co_tensor_list = list()

for (j in 1:4) {
  
  co_tensor = array(0, dim = rep(length(sel_author),3))
  
  cat("co_tensor ", j, "\n")
  simplex_j = sel_simplex_list[[j]]
  
  # coauthorship tensor 1
  for (i in 1:length(simplex_j)) {
    
    aut_ind = which(sel_author %in% simplex_j[[i]])
    
    if(length(aut_ind) == 1){
      #next
      ind_arr = t(as.matrix(rep(aut_ind, 3)))
    }else if(length(aut_ind) == 2){
      #ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)[2:7,]
      ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)
    }else if(length(aut_ind) > 2){
      # aut_ind = c(1,2,3)
      # ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = T)[-c(1,14,27),]
      ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = T)
      #ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = F)
    }
    
    # for (j in 1:dim(ind_arr)[1]) {
    #   co_tensor1[ind_arr[j,][1], ind_arr[j,][2], ind_arr[j,][3]] = 1
    # }
    # 
    co_tensor[ind_arr] = co_tensor[ind_arr] + 1
    
  }
  co_tensor_list[[j]] = co_tensor
}

# Analysis and Matching ---------------------------------------------

comb = combinations(n = 4, r = 2, v = 1:4)
pi_star = list(1:length(sel_author), 1:length(sel_author), 1:length(sel_author))

# record test results 
cov_vec = rep(0, dim(comb)[1])
ame_mat = array(0, dim = c(6, dim(comb)[1])) # 6 methods, combs

for (i in 1:dim(comb)[1]) {
  
  cat("this is combination ", comb[i, ], " \n")
  num1 = comb[i,1]
  num2 = comb[i,2]
  
  # record correlation
  cov_vec[i] = cor(co_tensor_list[[num1]], co_tensor_list[[num2]])
  
  # our method: initialization
  mm_ini = W1_initial_mm(co_tensor_list[[num1]],co_tensor_list[[num2]])
  ame_mat[1, i] = mean(matching_error(mm_ini, pi_star))
  
  # our method: alternating
  mm_re = alter_cleanup_mm(co_tensor_list[[num1]],co_tensor_list[[num2]],mm_ini, 20)
  ame_mat[2, i] = mean(matching_error(mm_re, pi_star)) 
  
  cat("done ours \n")
  
  # ss
  pi_ss = slice_sum_mm(co_tensor_list[[num1]],co_tensor_list[[num2]])
  ame_mat[3, i] = mean(matching_error(pi_ss, pi_star)) 
  
  # ss+ 
  pi_ss_re = alter_cleanup_mm(co_tensor_list[[num1]],co_tensor_list[[num2]],pi_ss, 20)
  ame_mat[4, i] = mean(matching_error(pi_ss_re, pi_star)) 
  
  cat("done ss \n")
  
  # ts 
  pi_ts = spectral_mm(co_tensor_list[[num1]],co_tensor_list[[num2]])
  ame_mat[5, i] = mean(matching_error(pi_ts, pi_star)) 
  
  # ts+ 
  pi_ts_re = alter_cleanup_mm(co_tensor_list[[num1]],co_tensor_list[[num2]],pi_ts, 20)
  ame_mat[6, i] = mean(matching_error(pi_ts_re, pi_star)) 
  
  cat("done ts \n")
}

sqrt(1 - cov_vec^2)
ame_mat
ame_mat0 = ame_mat

# obtain slice-wise correlation
cov_mat = array(0, dim = c(length(sel_author), dim(comb)[1]))

for (i in  1:dim(comb)[1]) {
  cat("this is combination ", comb[i, ], " \n")
  num1 = comb[i,1]
  num2 = comb[i,2]
  
  # record slice-wise correlation
  cor_vec = c()
  for (j in 1:length(sel_author)) {
    cor_vec = c(cor_vec, cor(as.vector(co_tensor_list[[num1]][,,j]), as.vector(co_tensor_list[[num2]][,,j])))
  }
  cov_mat[,i] = cor_vec
}

# Plotting -------------------------------------------------------

#  plot the error
# dat = data.frame(AME = c(ame_mat[1,], ame_mat[2,], ame_mat[3,],ame_mat[4,],ame_mat[5,],ame_mat[6,]),
#                  Year = rep( c("09-10", "09-11", "09-12","10-11","10-12","11-12"), times = 6 ),
#                  Algorithm = rep( c("W1", "W1+", "SS", "SS+", "TS", "TS+"), each = 6 ),
#                  Noise = rep(cov_vec, times = 6))
# dat$Year = as.factor(dat$Year)
# dat$Algorithm = as.factor(dat$Algorithm)

# sel_color = c("#F3C937", "#5D8722", "#2677A5")
ame_mat = ame_mat0
order_ind = order(ame_mat[1,])

ame_mat = ame_mat[,order_ind]
year = c("2009vs2010", "2009vs2011", "2009vs2012","2010vs2011","2010vs2012","2011vs2012")
year = year[order_ind]
noise = sqrt(1 - cov_vec^2)[order_ind]
sel_color = c("#F3C937", "#90A74A", "#639BC1")

dat = data.frame(AME = c(ame_mat[1,], ame_mat[3,],ame_mat[5,]),
                 Year = rep( year, times = 3 ),
                 Algorithm = rep( c("W1", "SS", "TS"), each = 6 ),
                 Noise = rep(noise, times = 6))
dat$Year = factor(dat$Year, levels =  year)
dat$Algorithm = factor(dat$Algorithm, levels =  c("W1", "SS", "TS"))


co_plot = ggplot(dat) +
  geom_bar(aes(x = Year, y = AME, fill = Algorithm ), position=position_dodge(), size=.3, stat='identity')+
  coord_cartesian(ylim=c(0.65,1))+
  #geom_line(aes(x = Year, y = Noise + 0.3), group = 1, linetype = "dashed", color = "#C53211", lwd = 1) +
  #geom_point(size = 3,aes(x = Year, y = Noise + 0.3), shape = 16,  color = "#C53211")+
  scale_fill_manual(values = sel_color) +
  ylab("average matching error") + xlab('pair of years') +
  #scale_y_continuous( sec.axis = sec_axis(~.-0.3, name = "estimated noise level"))+
  theme_light()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16),
        plot.title = element_text(hjust = 0.5,size = 11),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        #axis.text.x = element_text(angle = -15, hjust = -0.1)
        )
co_plot

# slice-wise correlation boxplot
year = c("2009vs2010", "2009vs2011", "2009vs2012","2010vs2011","2010vs2012","2011vs2012")
year = year[order_ind]
box_dat = data.frame(cor = c(cov_mat[,order_ind[1]], cov_mat[,order_ind[2]], cov_mat[,order_ind[3]],
                             cov_mat[,order_ind[4]],cov_mat[,order_ind[5]],cov_mat[,order_ind[6]]),
                     Year = rep(year, each = length(sel_author)))
box_dat$Year = factor(box_dat$Year, levels =  year)
box_dat$noise =  sqrt(1 - box_dat$cor^2)

box_dat = box_dat[!is.na(box_dat$cor),]

co_box = ggplot(box_dat, aes(x = Year, y = noise))+
  geom_boxplot(outlier.colour="red", 
               outlier.shape=8,
               outlier.size=4,  
               fill="blue",
               alpha=0.2,
               # Notch?
               notch=TRUE,
               notchwidth = 0.8,)+
  xlab('pair of years')+ ylab("slice-wise estimated noise level")+
  theme_light()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=16),
        plot.title = element_text(hjust = 0.5,size = 11),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        #axis.text.x = element_text(angle = -15, hjust = -0.1)
  )
  # geom_violin(trim=FALSE) + 
  # stat_summary(fun.y=mean, geom="point", shape=20, size= 10, color="red", fill="red")
co_box

#summary(as.vector(sqrt(1 - box_dat$cor^2)))

#apply(cov_mat, 2, mean, na.rm = T)

pdf("figures/coauthor.pdf", width = 16, height = 5)
plot_grid(co_plot, co_box, rel_widths = c(1.5,1), labels = c('a', 'b'), label_size = 20)
dev.off()

# cor_vec = c()
# for (i in 1:dim(co_tensor_list[[1]])[3]) {
#   cor_vec = c(cor_vec, cor(as.vector(co_tensor_list[[3]][,,i]), as.vector(co_tensor_list[[4]][,,i])))
# }
# summary(cor_vec)

