#### Sushi Example ####

sushiA_score <- as.matrix(read.csv("~/Desktop/Paper2/Sushi/sushi3-2016/sushiA_score.csv")[,-1])
sushiA_order <- as.matrix(read.csv("~/Desktop/Paper2/Sushi/sushi3-2016/sushiA_order.csv")[,-1])
X <- sushiA_score
Pi <- sushiA_order
M <- 4
load("~/Desktop/Paper2/Sushi/sushi_priorpred_new.RData")
a <- results$a
b <- results$b
gamma1 <- results$gamma1
gamma2 <- results$gamma2
rm(results,sushiA_order,sushiA_score)


lambda <- 7
gamma_hyp1 <- 2
gamma_hyp2 <- 2

# hist(rgamma(10000,gamma_hyp1,gamma_hyp2),main="Prior Histogram on Gamma",xlab="gamma")
# hist(rpois(10000,lambda)+1,main="Prior Histogram on K",xlab="K")
plot_Kplus <- matrix(NA,nrow=0,ncol=3)
for(gamma in c(0.01,0.5,1,2,3,4)){
  print(gamma)
  pmfstatic2 <- nClusters(Kplus=1:20,N=nrow(X),type="static",gamma=gamma,maxK=50)
  dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = lambda))
  plot_Kplus <- rbind(plot_Kplus,matrix(c(rep(gamma,20),1:20,dens),ncol=3))
}
plot_Kplus <- as.data.frame(plot_Kplus)
names(plot_Kplus) <- c("gamma","K+","density")
plot_Kplus$gamma <- factor(plot_Kplus$gamma)
levels(plot_Kplus$gamma) <- c(expression("gamma: 0.01"),expression("gamma: 0.5"),
                              expression("gamma: 1"),expression("gamma: 2"),
                              expression("gamma: 3"),expression("gamma: 4"))
p1 <- ggplot(plot_Kplus,aes(`K+`,density))+geom_line()+geom_point()+
  facet_wrap(. ~gamma,labeller = label_parsed)+
  scale_x_continuous(limits = c(0,20),breaks=seq(0,20,by=2))+ylab("Prior Density")
p1
# ggsave("~/Desktop/PriorKplus_Sushi.pdf",p1)


# res <- btlb_mfm(Pi=Pi,X=X,M=M,lambda=lambda,a=a,b=b,
#                 gamma1=gamma1,gamma2=gamma2,gamma_hyp1=gamma_hyp1,gamma_hyp2=gamma_hyp2,
#                 Pi_full=NULL,mh_pjk = .20, mh_thetak = 5,mh_gamma=0.5,
#                 startK=5,max_iters=2,mh_iters=2,burn=0,thin=1,seed)


load("/Users/pearce790/Desktop/sushi_mfmm_TS_startK1_iters1000_seed1.RData")
res <- posterior
rm(posterior)


p2<-ggplot(data=data.frame(Kplus=res$Kplus),aes(Kplus,y=..density..))+
  geom_histogram()+
  #geom_histogram(breaks=c(8.5,9.5))+xlim(c(8,10))+
  xlab(expression("K+"))+ylab("Density")
p3<-ggplot(data=data.frame(gamma=res$gamma),aes(x=gamma,y=..density..))+
  #xlim(c(0,5))+
  geom_histogram(bins=20)+ylab("Density")+xlab(expression(gamma))
grid.arrange(p2,p3,nrow=1)
# ggsave("~/Desktop/res1_Sushi.pdf",grid.arrange(p2,p3,nrow=1),
#        width=11,height=4)


Zhat <- salso(x=res$Z[seq(1,5000,length=1000),seq(1,5000,length=500)],
              loss="binder",maxNClusters=max(res$Z))
sumZhat <- summary(Zhat)
plotZhat <- melt(sumZhat$psm[sumZhat$order,rev(sumZhat$order)])
p4 <- ggplot(plotZhat,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  scale_x_continuous(breaks=NULL,expand = c(0, 0))+
  scale_y_continuous(breaks=NULL,expand = c(0, 0))+
  labs(fill="Cluster\nSimilarity")+theme_bw()+
  xlab(NULL)+ylab(NULL)
# ggsave("~/Desktop/res2_Sushi.pdf",p4,width=8,height=7)





# 
# plot_p <- na.exclude(reshape2::melt(res$p[,1:2,res$Kplus %in% c(1,2)]))
# plot_p$Var2 <- factor(plot_p$Var2,levels=2:1,labels=1:2)
# p4<-ggplot(plot_p,aes(x=factor(Var1),y=value,fill=factor(Var2)))+geom_violin()+
#   theme(legend.position="none")+ylim(c(0,1))+
#   xlab("Item Quality Parameter")+ylab("Density")
# plot_theta <- na.exclude(reshape2::melt(res$theta[res$Kplus %in% c(1,2),1:2]))
# plot_theta$Var2 <- factor(plot_theta$Var2,levels=2:1,labels=1:2)
# p5<-ggplot(plot_theta,aes(x=factor(Var2),y=value,fill=factor(Var2)))+geom_violin()+
#   ylim(c(0,40))+xlab(expression(theta))+ylab("Density")+labs(fill="Cluster")
# 
# grid.arrange(p2,p1,p3,nrow=1,widths=c(0.25,0.25,0.5))
# grid.arrange(p4,p5,nrow=1,widths=c(.8,0.2))
# # ggsave("~/Desktop/AIBS_res1.pdf",grid.arrange(p2,p1,p3,nrow=1,widths=c(0.25,0.25,0.5)),
# #        width=11,height=4)
# # ggsave("~/Desktop/AIBS_res2.pdf",grid.arrange(p4,p5,nrow=1,widths=c(0.8,0.2)),
# #        width=11,height=4)
# 
# 
# View(plot_p %>% group_by(Var1,Var2) %>%summarize(mean(value)))
# plot_theta %>% group_by(Var2) %>% summarize(mean(value))
# 
# par(mfrow=c(2,3))
# plot(res$accept_p,ylim=c(0,1),type="l",ylab="Accept Prob for p")
# plot(res$accept_theta,ylim=c(0,1),type="l",ylab="Accept Prob for theta")
# plot(res$accept_gamma,ylim=c(0,1),type="l",ylab="Accept Prob for gamma")
# plot(res$K,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="K")
# plot(res$Kplus,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="Kplus")
# plot(res$gamma,type="l",ylim=c(0,max(res$gamma)+1),ylab="gamma")
# ggplot(reshape2::melt(res$pi),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
#   geom_line()+ylim(c(0,1))+theme(legend.position="bottom")+
#   ylab("Class Proportion Estimate")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: Class Proportions, pi")
# ggplot(reshape2::melt(res$Z),aes(x=Var1,y=jitter(value,.5),group=Var2,color=factor(Var2)))+
#   geom_line()+theme(legend.position="none")+ylim(c(0,max(res$Kplus)+1))+
#   ylab("Class Membership Indicators")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: Class Memberships, Z")
# ggplot(reshape2::melt(res$theta),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
#   geom_line()+theme(legend.position="bottom")+ylim(c(0,max(res$theta)+1))+
#   ylab("theta")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: theta")
# plot_p <- reshape2::melt(res$p)
# plot_p$jk <- as.factor(paste0(plot_p$Var1,"_",plot_p$Var2))
# ggplot(plot_p,aes(x=Var3,y=value,group=jk,color=factor(Var2)))+
#   facet_wrap(vars(Var1))+
#   geom_line()+theme(legend.position="bottom")+ylim(c(0,1))+
#   ylab("p")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: p")
# #save.image("~/BTLBinomial/AIBS_Panel2_res.RData")
# 
# 
# names(res)
# res$thin
# par(mfrow=c(2,2))
# plot(res$accept_p,ylim=c(0,1),type="l",ylab="Accept Prob for p")
# plot(res$accept_theta,ylim=c(0,1),type="l",ylab="Accept Prob for theta")
# plot(res$K,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="K")
# plot(res$Kplus,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="Kplus")
# ggplot(reshape2::melt(res$pi),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
#   geom_line()+ylim(c(0,1))+theme(legend.position="bottom")+
#   ylab("Class Proportion Estimate")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: Class Proportions, pi")
# ggplot(reshape2::melt(res$Z),aes(x=Var1,y=jitter(value,.5),group=Var2,color=factor(Var2)))+
#   geom_line()+theme(legend.position="none")+ylim(c(0,max(res$Kplus)+1))+
#   ylab("Class Membership Indicators")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: Class Memberships, Z")
# ggplot(reshape2::melt(res$theta),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
#   geom_line()+theme(legend.position="bottom")+ylim(c(0,max(res$theta)+1))+
#   ylab("theta")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: theta")
# plot_p <- reshape2::melt(res$p)
# plot_p$jk <- as.factor(paste0(plot_p$Var1,"_",plot_p$Var2))
# ggplot(plot_p,aes(x=Var3,y=value,group=jk,color=factor(Var2)))+
#   facet_wrap(vars(Var1))+
#   geom_line()+theme(legend.position="bottom")+ylim(c(0,1))+
#   ylab("p")+xlab("Iteration (after burn/thin)")+
#   labs(color="Class")+ggtitle("Trace Plot: p")







