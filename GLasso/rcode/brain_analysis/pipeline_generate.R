#### Brain data analysis (Pipeline) ####
#### Jiaxin Hu 03/07/21 ####

# Here is the pipeline to generate ss data, implement graphical check,
# check the Matlab results and visualize the results. 

##### loading the dependency ####

library(plot.matrix)
library(RColorBrewer)
library(rmatio)

# load full gtex data 
load("corrected_expression_list_noY.RData")
# load preprocessed data
load("brain_final.rdata")

# source the functions for generation 
source("generate_funcs.R")

##### pipeline to generate various kinds of ss data and figure #####

# We consider five ways to generate the sample covariance matrices and their figures
# We need five folders to store the covariance matrices and figures
# the list of 13 sample covariance matrices is stored as a single .mat file
# the 13 figures for the covariance matrices are stored in a folder named ss_figure

# Following are the path for the five ways:

pre_pro = "./pre_pro"# path for pre-processed data
gtex_pre = "./gtex_pre"# path for gtex data with pre-selected genes
gtex_top = "./gtex_top"# path for gtex data with top variance genes
gtex_varvar = "./gtex_varvar" # path for gtex data with top var-var genes
gtex_meanvar = "./gtex_meanvar"# path for gtex data with top mean-var genes


num01 = c(rep(0,9),rep("",4)) # the files should be named with index 01,02 ....,09,10,11
brain_tissue_name = colnames(sample[1,,]) 


### 1. pre-processed data
sample_list = list()
for (i in 1:13) {
  sample_list[[i]] = sample[,,i]
}

# generate sample covariance data
generate_ss_data(sample_list,rownames(sample[,,1]),pre_pro)

# graphical check
pre_pro_ss = read.mat(paste0(pre_pro, "/ss_data.mat"))[[1]]
names(pre_pro_ss) = brain_tissue_name
generate_ss_figure(pre_pro_ss, paste0(pre_pro,"/ss_figure"))


### 2. gtex data with pre-selected genes 
brain_tissue_index = c(14,33,39:49)

gtex_brain = list()
for (i in 1:13) {
  gtex_brain[[i]] = corrected_expression_list[[brain_tissue_index[i]]]
}

# generate sample covariance data
generate_ss_data(gtex_brain,rownames(sample[,,1]),gtex_pre)

# graphical check
gtex_pre_ss = read.mat(paste0(gtex_pre, "/ss_data.mat"))[[1]]
names(gtex_pre_ss) = brain_tissue_name
generate_ss_figure(gtex_pre_ss, paste0(gtex_pre,"/ss_figure"))


### 3. gtex data with top variance genes

# get top variance genes
gtex_brain_combine = gtex_brain[[1]]
for (i in 2:13) {
  gtex_brain_combine = cbind(gtex_brain_combine, gtex_brain[[i]])
}

gtex_brain_gene_var = apply(gtex_brain_combine, 1, var)
brain_top_var_index = order(gtex_brain_gene_var,decreasing = T)[1:300]

# generate sample covariance data
generate_ss_data(gtex_brain,brain_top_var_index,gtex_top)

# graphical check
gtex_top_ss = read.mat(paste0(gtex_top, "/ss_data.mat"))[[1]]
names(gtex_top_ss) = brain_tissue_name
generate_ss_figure(gtex_top_ss, paste0(gtex_top,"/ss_figure"))


### 4. gtex data with top var-var genes
brain_gene_tissue_var = matrix(0,ncol = 13,nrow = dim(corrected_expression_list[[1]])[1])
for (i in 1:13) {
  t_exp = gtex_brain[[i]]
  brain_gene_tissue_var[,i] = apply(t_exp, 1, var)
}

brain_gene_var_var = apply(brain_gene_tissue_var,1,var)
brain_top_var_var_index = order(brain_gene_var_var,decreasing = T)[1:300]

# generate sample covariance data
generate_ss_data(gtex_brain,brain_top_var_var_index,gtex_varvar)

# graphical check
gtex_varvar_ss = read.mat(paste0(gtex_varvar, "/ss_data.mat"))[[1]]
names(gtex_varvar_ss) = brain_tissue_name
generate_ss_figure(gtex_varvar_ss, paste0(gtex_varvar,"/ss_figure"))


### 5. gtex data with top mean-var genes
brain_gene_tissue_mean = matrix(0,ncol = 13,nrow = dim(corrected_expression_list[[1]])[1])
for (i in 1:13) {
  t_exp = gtex_brain[[i]]
  brain_gene_tissue_mean[,i] = apply(t_exp, 1, mean)
}

brain_gene_mean_var = apply(brain_gene_tissue_mean,1,var)
brain_top_mean_var_index = order(brain_gene_mean_var, decreasing = T)[1:300]

# generate sample covariance data
generate_ss_data(gtex_brain,brain_top_mean_var_index,gtex_meanvar)

# graphical check
gtex_meanvar_ss = read.mat(paste0(gtex_meanvar, "/ss_data.mat"))[[1]]
names(gtex_meanvar_ss) = brain_tissue_name
generate_ss_figure(gtex_meanvar_ss, paste0(gtex_meanvar,"/ss_figure"))



######## Check Matlab results ####### 

# first we need the result path
U_path = # Revised!
Theta_path = # Revised!
Theta_fpath = # Revised! The path to store the graphics for Theta
  
  
### membership results
U = read.csv(U_path,header = F)

r = 3
cluster_result =list()
for (i in 1:r) {
  cluster_result[[i]] = brain_tissue_name[ U[,i] != 0 ]
}

# group with global pattern
cluster_result[[(r+1)]] = brain_tissue_name[rowSums(U) == 0]

# print the results
for (i in 1:(r+1)) {
  cat("This is cluster ",i, "\n", cluster_result[[i]], "\n")
}


### Network(Theta) results
Theta = read.csv(Theta_path,header = F)
Theta = as.matrix(Theta) # data.frame -> matrix

# number of non-zero entries
# change != 0 -->  abs() > a to set the threshold 
# if there are too many non-zero entries 
sum(Theta - diag(diag(as.matrix(Theta))) != 0 )

# summary of all entries
summary(as.vector(Theta1))

# graphical check
filename = paste0(Theta_fpath,"Theta.jpeg")
jpeg(filename,width = 800, height = 800)
plot(Theta, col = c(col=brewer.pal(n = 3, name = "RdBu")), breaks = c(-4,-0.0001,0.0001,4),
     main = "Theta",border=NA)
dev.off()


