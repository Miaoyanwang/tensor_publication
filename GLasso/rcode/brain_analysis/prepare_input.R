#### Brain data analysis (Pipeline) ####
#### Jiaxin Hu 03/07/21 ####

# Latest update: Jiaxin Hu 03/21/21

# Here is the pipeline to generate ss data, implement graphical check,
# check the Matlab results and visualize the results. 

##### loading the dependency ####

library(plot.matrix)
library(RColorBrewer)
library(rmatio)

# load full gtex data 
load("corrected_expression_list_noY.RData")

# source the functions for generation 
source("functions.R")

### final choise. gtex data with top var-var genes

num01 = c(rep(0,9),rep("",4)) # the files should be named with index 01,02 ....,09,10,11
brain_tissue_index = c(14,33,39:49)
brain_tissue_name = names(corrected_expression_list)[brain_tissue_index]

gtex_brain = list()
for (i in 1:13) {
    gtex_brain[[i]] = corrected_expression_list[[brain_tissue_index[i]]]
}

names(gtex_brain)=brain_tissue_name
brain_gene_tissue_var = matrix(0,ncol = 13,nrow = dim(corrected_expression_list[[1]])[1])
for (i in 1:13) {
    t_exp = gtex_brain[[i]]
    brain_gene_tissue_var[,i] = apply(t_exp, 1, var)
}

brain_gene_var_var = apply(brain_gene_tissue_var,1,var)
brain_top_var_var_index = order(brain_gene_var_var,decreasing = T)[1:500]

# generate sample covariance data
nvector = c()# get sample size vector
for (i in c(14,33,39:49)) {
  nvector = c(nvector, dim(corrected_expression_list[[i]])[2])
}
generate_ss_data(gtex_brain,brain_top_var_var_index,nvector,"input/")

# graphical check
gtex_varvar_ss = read.mat(paste0("input/", "/ss_data.mat"))[[1]]
names(gtex_varvar_ss) = brain_tissue_name
generate_ss_figure(gtex_varvar_ss, paste0("","Figure/"))


