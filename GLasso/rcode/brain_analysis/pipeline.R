#### Brain tissue analysis ####
#### Jiaxin Hu 03/01/21 ####

# Here is the pipeline to generate ss data, implement graphical check,
# check the Matlab results and visualize the results. 

# Remember to revise these paths!

# under each of the following path, we need two folders: ss_data, and ss_figure
pre_pro = # path for pre-processed data
gtex_pre = # path for gtex data with pre-selected genes
gtex_top = # path for gtex data with top variance genes
gtex_varvar = # path for gtex data with top var-var genes
gtex_meanvar = # path for gtex data with top mean-var genes


library(plot.matrix)
library(RColorBrewer)


# load full gtex data 
load("corrected_expression_list_noY.RData")
# load preprocessed data
load("brain_final.rdata")


######## Generate different kinds of ss data and ss figures ####### 

num01 = c(rep(0,9),rep("",4))
# the files should be named with index 01,02 ....,09,10,11
# otherwise Matlab will read the files in a different order.

brain_tissue_name = colnames(sample[1,,]) 

generate_ss_figure = function(exp_data, select_gene, path){
  # exp_data should be a list of gene expression data
  # selec_genes is the index for the selected genes
  
  for (i in 1:13) {
    cut_exp = exp_data[[i]][select_gene,]
    ss = cov(t(cut_exp))
    
    # write ss_data
    filename = paste0(path,"/ss_data/ss_",num01[i],i,".csv")
    write.csv(ss,filename)
    
    # covariance matrix graph
    g_filename = paste0(path,"/ss_figure/ss_f",num01[i],i,".jpeg")
    jpeg(g_filename,width = 800, height = 800)
    plot(ss, col = c("darkorange4","darkorange3",col=brewer.pal(n = 4, name = "RdBu"),"royalblue3","royalblue4"), breaks = c(-10,seq(-1.5,1.5,0.5),10),
         main = paste0(brain_tissue_name[i]),border=NA)
    dev.off()
  }
}

### 1. pre-processed data
sample_list = list()
for (i in 1:13) {
  sample_list[[i]] = sample[,,i]
}

generate_ss_figure(sample_list,select_gene = 1:(dim(sample[,,1])[1]), pre_pro)

test = "/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer/test"
### 2. gtex data with pre-selected genes 
brain_tissue_index = c(14,33,39:49)
pre_brain_gene = rownames(sample[,,1])

gtex_brain = list()
for (i in 1:13) {
  gtex_brain[[i]] = corrected_expression_list[[brain_tissue_index[i]]]
}

generate_ss_figure(gtex_brain,pre_brain_gene,gtex_pre)

### 3. gtex data with top variance genes

# get top variance genes
gtex_brain_combine = gtex_brain[[1]]
for (i in 2:13) {
  gtex_brain_combine = cbind(gtex_brain_combine, gtex_brain[[i]])
}

gtex_brain_gene_var = apply(gtex_brain_combine, 1, var)
brain_top_var_index = order(gtex_brain_gene_var,decreasing = T)[1:300]

generate_ss_figure(gtex_brain,brain_top_var_index,gtex_top)


### 4. gtex data with top var-var genes
brain_gene_tissue_var = matrix(0,ncol = 13,nrow = dim(corrected_expression_list[[1]])[1])
for (i in 1:13) {
  t_exp = gtex_brain[[i]]
  brain_gene_tissue_var[,i] = apply(t_exp, 1, var)
}

brain_gene_var_var = apply(brain_gene_tissue_var,1,var)
brain_top_var_var_index = order(brain_gene_var_var)[1:300]

generate_ss_figure(gtex_brain,brain_top_var_var_index,gtex_varvar)

### 5. gtex data with top mean-var genes
brain_gene_tissue_mean = matrix(0,ncol = 13,nrow = dim(corrected_expression_list[[1]])[1])
for (i in 1:13) {
  t_exp = gtex_brain[[i]]
  brain_gene_tissue_mean[,i] = apply(t_exp, 1, mean)
}

brain_gene_mean_var = apply(brain_gene_tissue_mean,1,var)
brain_top_mean_var_index = order(brain_gene_mean_var)[1:300]

generate_ss_figure(gtex_brain,brain_top_mean_var_index,gtex_meanvar)


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


