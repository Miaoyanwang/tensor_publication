##### vertical barplots & preparation for network #####

# Origin:  03/12/21, Lastest update: Jiaxin Hu 03/14/21

## Here is the example to make barplot for membership matrix U
## The code is originally used to plot a gene enrichment plot for the final project of BMI 826
## The original figure is an error bar figure using  geom_errorbar()

library(ggplot2)


# check the clustering result quickly
U = read.csv("output_r3rho1500/U_r3_rho1500.csv",header = TRUE)
rownames(U)=U[,1]
U=U[,-1]

brain_tissue_name=rownames(as.matrix(U))
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

# for group 1
pdf("Figure/bar1.pdf",width=6,height=10)
gg1 = ggplot(U_df, aes(x = tissue_name, y = U4, fill = Tissues)) +
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
  ggtitle(label = "Global") +
  coord_flip()

gg1
dev.off()



########### prepare node file for network plotting
Theta1=read.csv("output_r3rho1500/Theta_1_r3_rho1500.csv")
rownames(Theta1)=Theta1[,1]
Theta1=as.matrix(Theta1[,-1])
Theta2=read.csv("output_r3rho1500/Theta_2_r3_rho1500.csv")
rownames(Theta2)=Theta2[,1]
Theta2=as.matrix(Theta2[,-1])
Theta3=read.csv("output_r3rho1500/Theta_3_r3_rho1500.csv")
rownames(Theta3)=Theta3[,1]
Theta3=as.matrix(Theta3[,-1])

Theta0=read.csv("output_r3rho1500/Theta0_r3_rho1500.csv")
rownames(Theta0)=Theta0[,1]
Theta0=as.matrix(Theta0[,-1])

index1=which((Theta1 - diag(diag(Theta1)))!=0,arr.ind=T)
index2=which((Theta2 - diag(diag(Theta2)))!=0,arr.ind=T)
#index3=which((Theta3 - diag(diag(Theta3)))!=0,arr.ind=T)

# index=rbind(index1,index2)
# index=index[index[,1]-index[,2]!=0,]

index = unique(c(index1[,1], index1[,2], index2[,1], index2[,2]))
gene_name=row.names(Theta1) 

Theta = Theta0
network = NULL
for (i in 2:length(index)) {
  for (j in 1:(i-1)) {
    network = rbind(network, c(gene_name[index[i]], gene_name[index[j]], abs(Theta[index[i],index[j]]), 
                               sign(Theta[index[i],index[j]]), Theta[index[i],index[j]]!= 0 ))
  }
}

colnames(network)=c("source","target","value","sign","exists")  
write.table(network,"Figure/network_input/network0.txt",row.names=F,quote=F)

# network=NULL
# for(i in 1:28){
#     for(j in 1:i){
#         if(abs(Theta1[index[i,1],index[j,1]])!=0)
#         network=rbind(network,c(gene_name[index[i,1]],gene_name[index[j,1]],abs(Theta1[index[i,1],index[j,1]]),sign(Theta1[index[i,1],index[j,1]]),"TRUE"))
#         else
#         network=rbind(network,c(gene_name[index[i,1]],gene_name[index[j,1]],abs(Theta1[index[i,1],index[j,1]]),sign(Theta1[index[i,1],index[j,1]]),"FALSE"))
#     }
# }
# colnames(network)=c("source","target","value","sign","exists")









