#### Functions to generate ss data and ss figures ####

### Jiaxin Hu 03/07021 ###

library(rmatio)
library(plot.matrix)
library(RColorBrewer)


generate_ss_data = function(exp_data, select_gene, path){
  # exp_data should be a list of gene expression data
  # selec_genes is the index for the selected genes
  # path is where we store the .mat file
  
  ss_list = list()
  for (i in 1:13) {
    cut_exp = exp_data[[i]][select_gene,]
    ss_list[[i]] = cov(t(cut_exp))
  }
  
  # write ss_list as a single .mat file
  output_list = list()
  output_list[[1]] = ss_list # then matlab recognize the mat as a single cell
  names(output_list) = "ss_data"
  
  filename = paste0(path,"/ss_data.mat")
  write.mat(output_list,filename)
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
