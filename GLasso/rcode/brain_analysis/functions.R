#### Functions to generate ss data and ss figures ####

### Jiaxin Hu 03/07/21 ###

# lastest update: Miaoyan Wang 03/13/21


library(rmatio)
library(plot.matrix)
library(RColorBrewer)
library("EnsDb.Hsapiens.v79")

generate_ss_data = function(exp_data, select_gene, path){
  # exp_data should be a list of gene expression data with tissue name
  # selec_genes is the index for the selected genes
  # path is where we store the .mat file
  
  # covert gene name
  cut_exp = exp_data[[1]][select_gene,]
  geneid = rownames(cut_exp)
  geneid = gsub("\\..*","",geneid)
  
  covert = ensembldb::select(EnsDb.Hsapiens.v79, keys= geneid, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
  # some gene can not be encoded
  index = geneid %in% covert$GENEID
  gene_name = geneid
  gene_name[index] = covert$SYMBOL
  

  ss_list = list()
  for (i in 1:13) {
    cut_exp = exp_data[[i]][select_gene,]
    
    rownames(cut_exp) = gene_name
    ss_list[[i]] = cov(t(cut_exp))
  }
  
  # write ss_list as a single .mat file
  output_list = list(ss_data = ss_list, gene_name = gene_name, tissue_name = names(exp_data))
  #output_list[[1]] = ss_list # then matlab recognize the mat as a single cell
  #names(output_list) = "ss_data"
  
  filename = paste0(path,"/ss_data.mat" )
  write.mat(output_list,filename, compression = F)
}

generate_ss_figure = function(ss_data,path){
  # ss_data should be a list of covariance matrices and the list name should be the tissue name
  # path is the folder where we store the figures
  
  list_name = names(ss_data)
  for (i in 1:length(list_name)) {
    g_filename = paste0(path,"/ss_f",num01[i],i,".jpeg")
    jpeg(g_filename,width = 800, height = 800)
    plot(ss_data[[i]], col = c("darkorange4","darkorange3",col=brewer.pal(n = 4, name = "RdBu"),"royalblue3","royalblue4"), breaks = c(-10,seq(-1.5,1.5,0.5),10),
         main = paste0(list_name[i]),border=NA)
    dev.off()
  }
  
}


## 03/11/21 update

# Code to generate varvar data
# load("/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/gene_full/corrected_expression_list_noY.RData")
# brain_tissue_index = c(14,33,39:49)
# 
# gtex_brain = list()
# for (i in 1:13) {
#   gtex_brain[[i]] = corrected_expression_list[[brain_tissue_index[i]]]
# }
# brain_tissue_name = colnames(sample[1,,]) 
#names(gtex_brain) = brain_tissue_name
# 
# write.csv(brain_tissue_name,"/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer/test1/varvar_withname/tissue_name.csv")

# a = c("Brain - Cortex")
# gsub(" - ", "",a)
# 
# brain_gene_tissue_var = matrix(0,ncol = 13,nrow = dim(corrected_expression_list[[1]])[1])
# for (i in 1:13) {
#   t_exp = gtex_brain[[i]]
#   brain_gene_tissue_var[,i] = apply(t_exp, 1, var)
# }
# 
# brain_gene_var_var = apply(brain_gene_tissue_var,1,var)
# brain_top_var_var_index = order(brain_gene_var_var,decreasing = T)[1:300]
# 
#generate_ss_data(gtex_brain,brain_top_var_var_index,"/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer/test1/varvar_withname")
