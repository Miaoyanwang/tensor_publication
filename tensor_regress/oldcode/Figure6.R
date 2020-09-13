## Code for Figure 6. Comparison of mean sqaured prediction error between our method and the other three previous methods
install.packages("tensorregress.tar.gz",repos = NULL, type="source")
source("simulation.R")
source("comparison.R")
library(tensorregress)

############# plot figures using saved data from prior simulation. If a new run of simulation is desired, please go to line #40 to run the code ##########

######## plot MSPE for three methods under different (signal, rank) setting #######
load("presaved/Figure6.RData")
library(ggplot2)
library(wesanderson)
library(patchwork)
d=c(20,40,60,80,100)
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))


signal_range=c(3,6)
new_color = c("#069AA0","#CCC591","#BCA455","#D6CFC4")
pdf("Figure6.pdf",width=21,height=3)
s=r=1;
data = data.frame(d = rep(d,4), Method = c(rep('STD (Our method)',length(d)),rep('HOLRR (Low-rank regression)',length(d)),rep('HOPLS (Partial least square)',length(d)),rep('TPG (Projection gradient)',length(d))),PMSE = c(final[s,r,,1],final[s,r,,2],final[s,r,,3],final[s,r,,4])/signal_range[s],sd = c(finalsd[s,r,,1],finalsd[s,r,,2],finalsd[s,r,,3],finalsd[s,r,,4])/signal_range[s])
data[,2]=factor(data[,2],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p1=ggplot(data, aes(x = d*400, y = PMSE)) + geom_line(aes(color = Method),size = 1.2) + scale_colour_manual(values = new_color) +
        geom_point(size = 2,aes(shape = Method)) + scale_shape_manual(values = c(16,5,17,4)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=10))+ xlab('Sample Size') + labs(title="Low Signal, Low Rank") + ylab('MSPE')+ theme(plot.title = element_text(hjust = 0.5,size = 11))+coord_cartesian(ylim = c(0, max(data[,3])+0.03))

p1=p1+geom_errorbar(aes(ymin=PMSE-sd, ymax=PMSE+sd),width=0.5,position=position_dodge(0.05))



r=2;
data = data.frame(d = rep(d,4), Method = c(rep('STD (Our method)',length(d)),rep('HOLRR (Low-rank regression)',length(d)),rep('HOPLS (Partial least square)',length(d)),rep('TPG (Projection gradient)',length(d))),PMSE = c(final[s,r,,1],final[s,r,,2],final[s,r,,3],final[s,r,,4])/signal_range[s],sd = c(finalsd[s,r,,1],finalsd[s,r,,2],finalsd[s,r,,3],finalsd[s,r,,4])/signal_range[s])
data[,2]=factor(data[,2],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p2=ggplot(data, aes(x = d*400, y = PMSE)) + geom_line(aes(color = Method),size = 1.2)  + scale_colour_manual(values = new_color)+
geom_point(size = 2,aes(shape = Method))+ scale_shape_manual(values = c(16,5,17,4)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=10))+ xlab('Sample Size') + labs(title="Low Signal, High Rank") + ylab('MSPE')+ theme(plot.title = element_text(hjust = 0.5,size = 11))+coord_cartesian(ylim = c(0,0.12))
p2=p2+geom_errorbar(aes(ymin=PMSE-sd,ymax=PMSE+sd),width=0.5,position=position_dodge(0.05))


s=2;r=1;
data = data.frame(d = rep(d,4), Method = c(rep('STD (Our method)',length(d)),rep('HOLRR (Low-rank regression)',length(d)),rep('HOPLS (Partial least square)',length(d)),rep('TPG (Projection gradient)',length(d))),PMSE = c(final[s,r,,1],final[s,r,,2],final[s,r,,3],final[s,r,,4])/signal_range[s],sd = c(finalsd[s,r,,1],finalsd[s,r,,2],finalsd[s,r,,3],finalsd[s,r,,4])/signal_range[s])
data[,2]=factor(data[,2],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p3=ggplot(data, aes(x = d*400, y = PMSE)) + geom_line(aes(color = Method),size = 1.2) +scale_colour_manual(values = new_color)+
geom_point(size = 2,aes(shape = Method))+ scale_shape_manual(values = c(16,5,17,4)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=10))+ xlab('Sample Size') + labs(title="High Signal, Low Rank") + ylab('MSPE')+ theme(plot.title = element_text(hjust = 0.5,size = 11))+coord_cartesian(ylim = c(0, 0.12))
p3=p3+geom_errorbar(aes(ymin=PMSE-sd,ymax=PMSE+sd),width=0.5,position=position_dodge(0.05))


s=r=2;
data = data.frame(d = rep(d,4), Method = c(rep('STD (Our method)',length(d)),rep('HOLRR (Low-rank regression)',length(d)),rep('HOPLS (Partial least square)',length(d)),rep('TPG (Projection gradient)',length(d))),PMSE = c(final[s,r,,1],final[s,r,,2],final[s,r,,3],final[s,r,,4])/signal_range[s],sd = c(finalsd[s,r,,1],finalsd[s,r,,2],finalsd[s,r,,3],finalsd[s,r,,4])/signal_range[s])
data[,2]=factor(data[,2],levels=c("STD (Our method)","HOLRR (Low-rank regression)","HOPLS (Partial least square)","TPG (Projection gradient)"))
p4=ggplot(data, aes(x = d*400, y = PMSE)) + geom_line(aes(color = Method),size = 1.2) +scale_colour_manual(values = new_color)+
geom_point(size = 2,aes(shape = Method))+ scale_shape_manual(values = c(16,5,17,4)) + theme(axis.text=element_text(size=12), axis.title=element_text(size=10))+ xlab('Sample Size') + labs(title="High Signal, High Rank") + ylab('MSPE')+ theme(plot.title = element_text(hjust = 0.5,size = 11))+coord_cartesian(ylim = c(0, max(data[,3])+0.05))
p4=p4+geom_errorbar(aes(ymin=PMSE-sd,ymax=PMSE+sd),width=0.5,position=position_dodge(0.05))

(p1|p2|p3|p4)
dev.off()

##################### end of plotting ##############################

### If a new run of simulation is desired, please run the code from here and save the results as .RData. Then run the above code to generate figures. ###
set.seed(0)
seed=0
## combinations of signal and rank
signal_range=c(3,6)
core_range=rbind(c(3,3,3),c(4,5,6))
d=c(20,40,60,80,100)

##
dup=30;
final=array(0,dim=c(2,2,5,4))
finalsd=array(0,dim=c(2,2,5,4))
err_res=err_holrr=err_pls=err_tpg=array(0,dim=c(2,2,5,dup))
dist="normal";


for(s in 1:2){ ## signal level
    for(r in 1:2){ ## rank
        for(i in 1:5){ ## dimension with informative mode
            whole_shape = c(d[i],20,20)
            
            for(n in 1:dup){ ## simulation replicates
        core_shape=core_range[r,]
        signal=signal_range[s]

            data=sim_data(whole_shape = whole_shape, core_shape=core_shape,p=c(0.4*d[i],0,0),dist=dist, dup=dup, signal=signal)
            res = tensor_regress(data$tsr[[n]],X_covar1 = data$X_covar1,core_shape=core_shape,Nsim=10, cons = 'non', dist = dist)### tensor_regress
    
        err_res[s,r,i,n]= mean((res$U - data$U)^2)
        
        ### holrr
        holrr = HOLRR(X = data$X_covar1, data$tsr[[n]], core_shape = core_shape) ## same rank as ours
        err_holrr[s,r,i,n] = mean((holrr$pre@data -data$U)^2)
        
        ## pls ## choose hyperparameter using the best one
        test=NULL
        for(R in 1:8){
        pls = HOPLS(X=data$X_covar1,data$tsr[[n]],R=R,Kn=c(4,6),tol=0.01) ## rank
        test=c(test,mean((pls$pre@data - data$U)^2))
        }
        err_pls[s,r,i,n]=min(test)
        
        ## Yu ## choose hyperparameter using the best one
        test=NULL
        for(R in 1:8){
            tpg=TPG(X=data$X_covar1,data$tsr[[n]],R=R,eta=.1,iter=15,tol=10-6)
            test=c(test,mean((tpg$U - data$U)^2))
        }
        err_tpg[s,r,i,n]=min(test)
    }
    final[s,r,i,1]=mean(err_res[s,r,i,])
    final[s,r,i,2]=mean(err_holrr[s,r,i,])
    final[s,r,i,3]=mean(err_pls[s,r,i,])
    final[s,r,i,4]=mean(err_tpg[s,r,i,])
    
    
    finalsd[s,r,i,1]=sd(err_res[s,r,i,])
    finalsd[s,r,i,2]=sd(err_holrr[s,r,i,])
    finalsd[s,r,i,3]=sd(err_pls[s,r,i,])
    finalsd[s,r,i,4]=sd(err_tpg[s,r,i,])
        }
    }
}

save(final,finalsd,file="Figure6.RData")
