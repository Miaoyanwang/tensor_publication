## Code for Figure 5. Comparison of mean sqaured prediction error between our method and the other three previous methods
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
source("comparison.R")
library(tensorregress)
library(ggplot2)

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line #40 to run the code ##########

######## plot MSPE for three methods under different (signal, rank) setting #######

load("presaved/Figure5.RData")
final = final_mode
finalsd = final_sd_mode

library(ggplot2)
library(patchwork)
new_color = c("#069AA0","#CCC591","#BCA455","#D6CFC4")

signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))
pdf("Figure.pdf",width=12,height=6)

s=1;r=2;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,])/signal_range[s],sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/signal_range[s],Method=rep(c('STD (Our method)','HOLRR (Low-rank regression)','HOPLS (Partial least square)','TPG (Projection gradient)'),3),Category=c(rep(1,4),rep(2,4),rep(3,4)))
data[,3]=factor(data[,3],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p2=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE")+coord_cartesian(ylim = c(0, 0.1)) +  labs(title = "Low Signal, High Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
    scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))

s=2;r=1;
data=data.frame(PMSE=c(final[s,r,1,],final[s,r,2,],final[s,r,3,])/signal_range[s],sd=c(finalsd[s,r,1,],finalsd[s,r,2,],finalsd[s,r,3,])/signal_range[s],Method=rep(c('STD (Our method)','HOLRR (Low-rank regression)','HOPLS (Partial least square)','TPG (Projection gradient)'),3),Category=c(rep(1,4),rep(2,4),rep(3,4)))
data[,3]=factor(data[,3],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p3=ggplot(data=data, aes(x=as.factor(Category),y=PMSE, fill=Method))+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd), width=.2,position=position_dodge(.9))+labs(x="Number of modes with available features",y="MSPE")+coord_cartesian(ylim = c(0, 0.1)) +  labs(title = "High Signal, Low Rank",size = 5) +theme(plot.title = element_text(hjust = 0.5,size = 11))+
    scale_fill_manual(values=new_color)+theme(axis.text=element_text(size=12),axis.title=element_text(size=10))


library(ggplot2)
library(patchwork)
d=c(20,40,60,80,100)
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))

final = final_sample
finalsd = final_sd_sample

s=1;r=2;
data = data.frame(d = rep(d,4), Method = c(rep('STD (Our method)',length(d)),rep('HOLRR (Low-rank regression)',length(d)),rep('HOPLS (Partial least square)',length(d)),rep('TPG (Projection gradient)',length(d))),PMSE = c(final[s,r,,1],final[s,r,,2],final[s,r,,3],final[s,r,,4])/signal_range[s],sd = c(finalsd[s,r,,1],finalsd[s,r,,2],finalsd[s,r,,3],finalsd[s,r,,4])/signal_range[s])
data[,2]=factor(data[,2],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p5=ggplot(data, aes(x = d*400, y = PMSE)) + geom_line(aes(color = Method),size = 1.2)  + scale_colour_manual(values = new_color)+
geom_point(size = 2,aes(shape = Method))+ scale_shape_manual(values = c(16,5,17,4)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=10))+ xlab('Sample Size') + labs(title="Low Signal, High Rank") + ylab('MSPE')+ theme(plot.title = element_text(hjust = 0.5,size = 11))+coord_cartesian(ylim = c(0,0.12))
p5=p5+geom_errorbar(aes(ymin=PMSE-sd,ymax=PMSE+sd),width=0.5,position=position_dodge(0.05))


s=2;r=1;
data = data.frame(d = rep(d,4), Method = c(rep('STD (Our method)',length(d)),rep('HOLRR (Low-rank regression)',length(d)),rep('HOPLS (Partial least square)',length(d)),rep('TPG (Projection gradient)',length(d))),PMSE = c(final[s,r,,1],final[s,r,,2],final[s,r,,3],final[s,r,,4])/signal_range[s],sd = c(finalsd[s,r,,1],finalsd[s,r,,2],finalsd[s,r,,3],finalsd[s,r,,4])/signal_range[s])
data[,2]=factor(data[,2],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p6=ggplot(data, aes(x = d*400, y = PMSE)) + geom_line(aes(color = Method),size = 1.2) +scale_colour_manual(values = new_color)+
geom_point(size = 2,aes(shape = Method))+ scale_shape_manual(values = c(16,5,17,4)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=10))+ xlab('Sample Size') + labs(title="High Signal, Low Rank") + ylab('MSPE')+ theme(plot.title = element_text(hjust = 0.5,size = 11))+coord_cartesian(ylim = c(0, 0.12))
p6=p6+geom_errorbar(aes(ymin=PMSE-sd,ymax=PMSE+sd),width=0.5,position=position_dodge(0.05))

(p3|p2)/
(p6|p5)
dev.off()
