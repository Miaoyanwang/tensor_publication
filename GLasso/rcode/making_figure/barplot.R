##### vertical barplots #####
### 03/12/21 ###

## Here is the example to make barplot for membership matrix U
## The code is originally used to plot a gene enrichment plot for the final project of BMI 826
## The original figure is an error bar figure using  geom_errorbar()

library(ggplot2)

## read data
load("corrected_expression_list_noY.RData")

brain_tissue_index = c(14,33,39:49)
brain_tissue_name = names(corrected_expression_list)[brain_tissue_index]

# check the clustering result quickly
U = read.csv("./results/r3rho1500/U_r3_rho1500.csv",header = F)

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

U_df = data.frame(tissue_name = tissue_name, U1 = U[,1], U2 = U[,2],U3 = U[,3],U4 = c(rep(0,12),1))
U_df$Tissues = rep(0,13)

U_df$Tissues[U_df$U1 > 0] = "Other Tissues"
U_df$Tissues[U_df$U1 < 0] = "Cortex"
U_df$Tissues[U_df$U2 > 0] = "Cerebellum"
U_df$Tissues[U_df$U2 < 0] = "Basal Ganglia"
U_df$Tissues[U_df$U3 > 0] = "Substantia nigra"
U_df$Tissues[U_df$U4 > 0] = "Global Only"

U_df$Tissues = as.factor(U_df$Tissues)
U_df = U_df[order(U_df$Tissues),]

selected_color = c("#6F9B3C","#F2AB1D",'#D26A1B',"#43ABE4","#214EB8","#C76AB2")
# check the color
barplot(1:length(selected_color), col=as.character(selected_color), names.arg = selected_color, las=2) 
U_df$color = rep(0,13)

for (i in 1:length(levels(U_df$Tissues))) {
  U_df$color[ U_df$Tissues %in% levels(U_df$Tissues)[i]  ] = selected_color[i]
}

# for group 1
gg1 = ggplot(U_df, aes(x = tissue_name, y = U1, fill = Tissues)) + 
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
  ggtitle(label = "Group 2") + 
  coord_flip()

gg1




