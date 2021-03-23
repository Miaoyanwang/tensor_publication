##### vertical barplots & preparation for network #####

# Origin:  03/12/21, Lastest update: Jiaxin Hu 03/21/21

## Here is the example to make barplot for membership matrix U and networks for precision matrix

## The code for barplot is originally used to plot a gene enrichment plot for the final project of BMI 826
## The original figure is an error bar figure using  geom_errorbar()

library(ggplot2)
library(igraph)
library(rmatio)

# load result data
result = read.mat("output_r3rho1500/result.mat")
brain_tissue_name = unlist(result$tissue_name)

# check the clustering result quickly
U = result$U
rownames(U) = brain_tissue_name


r = 3
cluster_result =list()
for (i in 1:r) {
  cluster_result[[i]] = brain_tissue_name[ U[,i] != 0 ]
}
cluster_result[[(r+1)]] = brain_tissue_name[rowSums(U) == 0]
for (i in 1:(r+1)) {
  cat("This is cluster ",i, "\n", cluster_result[[i]], "\n")
}


## bar plot visualization

tissue_name = brain_tissue_name
tissue_name = gsub("Brain - ", "",tissue_name)
tissue_name = gsub("\\s*\\([^\\)]+\\)","",tissue_name)

U_df = data.frame(tissue_name = tissue_name, U1 = U[,1], U2 = U[,2],U3 = U[,3],U4=rep(1,13)) ## modify based on outputs
U_df$Tissues = rep(0,13)

#### modify based on outputs
U_df$Tissues[U_df$U1 > 0] = "Other Tissues"
U_df$Tissues[U_df$U1 < 0] = "Cortex"
U_df$Tissues[U_df$U2 < 0] = "Cerebellum"
U_df$Tissues[U_df$U2 > 0] = "Basal Ganglia"
U_df$Tissues[U_df$U3 > 0] = "Substantia nigra"
U_df$Tissues[U_df$U3 < 0] = "Amygdala"


U_df$Tissues = as.factor(U_df$Tissues)
U_df = U_df[order(U_df$Tissues),]

selected_color = c("#6F9B3C","#F2AB1D","#D26A1B","#43ABE4","#214EB8","#C76AB2")
# check the color
barplot(1:length(selected_color), col=as.character(selected_color), names.arg = selected_color, las=2) 
U_df$color = rep(0,13)

for (i in 1:length(levels(U_df$Tissues))) {
  U_df$color[ U_df$Tissues %in% levels(U_df$Tissues)[i]  ] = selected_color[i]
}

# for group 0
pdf("Figure/bar3.pdf",width=6,height=10)
gg1 = ggplot(U_df, aes(x = tissue_name, y = U3, fill = Tissues)) +
  geom_bar(stat="identity", width=0.75)+
  scale_fill_manual(values = selected_color) +
  #scale_colour_manual(values = c("#6F9B3C","#F2AB1D","#214EB8","#43ABE4","#7776BC","#C76AB2")) +
  scale_x_discrete(limits=c(U_df$tissue_name))+
  theme_light()+
  theme(
    aspect.ratio = 3/1,
    axis.text.y = element_text(size = 20, color = as.character(U_df$color),
                               angle = 30),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    plot.title = element_text(hjust = 0.5,size = 20),
    legend.position = "none"
  )+
  ylim(-1,1) +
  ylab("Membership Loadings")+
  geom_hline(yintercept = 0, lty = "dotted", alpha = 0.8)+
  ggtitle(label = "Group 3") +
  coord_flip()

gg1
dev.off()



########### network plotting using igraph
# read data
gene_name = unlist(result$gene_name)

Theta0 = result$Theta0
rownames(Theta0) = gene_name; colnames(Theta0) = gene_name

Theta1 = result$T[[1]][[1]]
rownames(Theta1) = gene_name; colnames(Theta1) = gene_name

Theta2 = result$T[[2]][[1]]
rownames(Theta2) = gene_name; colnames(Theta2) = gene_name

Theta3 = result$T[[3]][[1]]
rownames(Theta3) = gene_name; colnames(Theta3) = gene_name

# get strong genes
index1=which((Theta1 - diag(diag(Theta1)))!=0,arr.ind=T)
index2=which((Theta2 - diag(diag(Theta2)))!=0,arr.ind=T)
index3=which((Theta3 - diag(diag(Theta3)))!=0,arr.ind=T)
index3=index3[order(abs(Theta3[index3]),decreasing = T)[1:20],] # first 10 strong connections in Theta3

index = unique(c(index1[,1], index1[,2], index2[,1], index2[,2],index3[,1], index3[,2]))

# get vertices
gene_name= gsub("\\_.*","",gene_name)
gene = as.data.frame(sort(gene_name[index]))

# get node-edge file
# for group 2
Theta = Theta2
network = NULL
for (i in 2:length(index)) {
  for (j in 1:(i-1)) {
    if( abs(Theta[index[i],index[j]])!=0){
      network = rbind(network, c(gene_name[index[i]], gene_name[index[j]], 
                                 abs(Theta[index[i],index[j]]), sign(Theta[index[i],index[j]])))
    }
  }
}
colnames(network) = c("from","to","value","sign")
network = as.data.frame(network)

# For Theta 0, cut off weak connecitons for better visualization
# For Theta_l, l = 1,2,3, ignore this step
# network = network[abs(as.numeric(network$value)) > 0.05,]



# plot
g1 = graph_from_data_frame(network, directed = F, vertices = gene)
width = as.numeric(network$value)*20 # For Theta0, change to *5
color = rep("#D4613E",length(network$sign)) # negative red
color[network$sign == 1] = "#466CA6" #positive blue

pdf("Figure/Theta0.pdf",width=10, height=10)
plot(g1, layout = layout.circle(g1),
     edge.width = width,
     edge.color = color,
     vertex.color = 'white',
     vertex.frame.color = "black")
dev.off()









