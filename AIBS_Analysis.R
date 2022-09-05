#### AIBS Reproducibility ####

## Load Data and Choose Hyperparameters
load("AIBS_2021Panel2.RData")
fit <- fitdistrplus::fitdist(data=as.vector(na.exclude(c(X)))/40,distr="beta",method="mme")
a <- as.numeric(fit$estimate[1])
b <- as.numeric(fit$estimate[2])
gamma1 <- 10
gamma2 <- 0.5
gamma_hyp1 <- 3
gamma_hyp2 <- 2
lambda <- 1

## Explore Influence of Hyperparameters on K and K+
plot_Kplus <- matrix(NA,nrow=0,ncol=3)
for(gamma in c(0.01,1,2,3)){
  pmfstatic2 <- nClusters(Kplus=1:6,N=nrow(X),type="static",gamma=gamma,maxK=50)
  dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = lambda))
  plot_Kplus <- rbind(plot_Kplus,matrix(c(rep(gamma,6),1:6,dens),ncol=3))
}
plot_Kplus <- as.data.frame(plot_Kplus)
names(plot_Kplus) <- c("gamma","K+","density")
plot_Kplus$gamma <- factor(plot_Kplus$gamma)
levels(plot_Kplus$gamma) <- c(expression("gamma: 0.01"),expression("gamma: 1"),
                              expression("gamma: 2"),expression("gamma: 3"))
hist(rgamma(1000,gamma_hyp1,gamma_hyp2),main="Prior Histogram on Gamma",xlab="gamma")
hist(rpois(1000,lambda)+1,main="Prior Histogram on K",xlab="K")
p1 <- ggplot(plot_Kplus,aes(`K+`,density))+geom_line()+geom_point()+
  facet_wrap(. ~gamma,labeller = label_parsed)+
  scale_x_continuous(breaks=1:10)+ylab("Prior Density")
p1
#ggsave("~/Desktop/PriorKplus_AIBS.pdf",p1,)


res <- btlb_mfm(Pi=Pi,X=X,M=M,Pi_full=Pi_full,
                lambda=lambda,a=a,b=b,gamma1=gamma1,gamma2=gamma2,
                gamma_hyp1=gamma_hyp1,gamma_hyp2=gamma_hyp2,
                mh_pjk = .05, mh_thetak = 5, mh_gamma = 1,
                startK=17,max_iters=1000,mh_iters=10,burn=.50,thin=5,
                seed = 1)


p1<-ggplot(data=data.frame(Kplus=res$Kplus),aes(Kplus))+
  geom_histogram(breaks=c(1.75,2.25,2.75,3.25))+
  scale_x_continuous(breaks=c(2,3))+
  xlab(expression("K+"))+ylab("Density")+
  scale_y_continuous(breaks=seq(0,1000,length=11),labels=c(seq(0,1,length=11)))
p2<-ggplot(data=data.frame(gamma=res$gamma),aes(x=gamma,y=..density..))+
  geom_histogram(bins=25)+xlim(c(0,4))+ylab("Density")+
  xlab(expression(gamma))
plotZ <- reshape2::melt(res$Z[res$Kplus != 3,])
plotZ$value <- factor(plotZ$value,levels=2:1,labels=1:2)
p3<-ggplot(data=plotZ,aes(Var2,fill=factor(value)))+
  geom_bar()+xlab("Judge")+ylab("Density")+
  scale_x_continuous(breaks=1:17)+
  scale_y_continuous(breaks=seq(0,max(plotZ$Var1),length=11),labels=seq(0,1,length=11))+
  labs(fill="Cluster")

plot_p <- na.exclude(reshape2::melt(res$p[,1:2,res$Kplus %in% c(1,2)]))
plot_p$Var2 <- factor(plot_p$Var2,levels=2:1,labels=1:2)
p4<-ggplot(plot_p,aes(x=factor(Var1),y=value,fill=factor(Var2)))+geom_violin()+
  theme(legend.position="none")+ylim(c(0,1))+
  xlab("Item Quality Parameter")+ylab("Density")
plot_theta <- na.exclude(reshape2::melt(res$theta[res$Kplus %in% c(1,2),1:2]))
plot_theta$Var2 <- factor(plot_theta$Var2,levels=2:1,labels=1:2)
p5<-ggplot(plot_theta,aes(x=factor(Var2),y=value,fill=factor(Var2)))+geom_violin()+
  ylim(c(0,40))+xlab(expression(theta))+ylab("Density")+labs(fill="Cluster")

grid.arrange(p2,p1,p3,nrow=1,widths=c(0.25,0.25,0.5))
grid.arrange(p4,p5,nrow=1,widths=c(.8,0.2))
# ggsave("~/Desktop/AIBS_res1.pdf",grid.arrange(p2,p1,p3,nrow=1,widths=c(0.25,0.25,0.5)),
#        width=11,height=4)
# ggsave("~/Desktop/AIBS_res2.pdf",grid.arrange(p4,p5,nrow=1,widths=c(0.8,0.2)),
#        width=11,height=4)


View(plot_p %>% group_by(Var1,Var2) %>%summarize(mean(value)))
plot_theta %>% group_by(Var2) %>% summarize(mean(value))

par(mfrow=c(2,3))
plot(res$accept_p,ylim=c(0,1),type="l",ylab="Accept Prob for p")
plot(res$accept_theta,ylim=c(0,1),type="l",ylab="Accept Prob for theta")
plot(res$accept_gamma,ylim=c(0,1),type="l",ylab="Accept Prob for gamma")
plot(res$K,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="K")
plot(res$Kplus,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="Kplus")
plot(res$gamma,type="l",ylim=c(0,max(res$gamma)+1),ylab="gamma")
ggplot(reshape2::melt(res$pi),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
  geom_line()+ylim(c(0,1))+theme(legend.position="bottom")+
  ylab("Class Proportion Estimate")+xlab("Iteration (after burn/thin)")+
  labs(color="Class")+ggtitle("Trace Plot: Class Proportions, pi")
ggplot(reshape2::melt(res$Z),aes(x=Var1,y=jitter(value,.5),group=Var2,color=factor(Var2)))+
  geom_line()+theme(legend.position="none")+ylim(c(0,max(res$Kplus)+1))+
  ylab("Class Membership Indicators")+xlab("Iteration (after burn/thin)")+
  labs(color="Class")+ggtitle("Trace Plot: Class Memberships, Z")
ggplot(reshape2::melt(res$theta),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
  geom_line()+theme(legend.position="bottom")+ylim(c(0,max(res$theta)+1))+
  ylab("theta")+xlab("Iteration (after burn/thin)")+
  labs(color="Class")+ggtitle("Trace Plot: theta")
plot_p <- reshape2::melt(res$p)
plot_p$jk <- as.factor(paste0(plot_p$Var1,"_",plot_p$Var2))
ggplot(plot_p,aes(x=Var3,y=value,group=jk,color=factor(Var2)))+
  facet_wrap(vars(Var1))+
  geom_line()+theme(legend.position="bottom")+ylim(c(0,1))+
  ylab("p")+xlab("Iteration (after burn/thin)")+
  labs(color="Class")+ggtitle("Trace Plot: p")
#save.image("~/BTLBinomial/AIBS_Panel2_res.RData")

