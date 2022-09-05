#### AIBS Reproducibility ####
source("~/BTLBinomial/source.R")

#### Load Data and Choose Hyperparameters ####
load("AIBS_2021Panel2.RData")
fit <- fitdistrplus::fitdist(data=as.vector(na.exclude(c(X)))/40,distr="beta",method="mme")
a <- as.numeric(fit$estimate[1])
b <- as.numeric(fit$estimate[2])
gamma1 <- 10
gamma2 <- 0.5
gamma_hyp1 <- 3
gamma_hyp2 <- 2
lambda <- 1

#### Explore Influence of Hyperparameters on K and K+ ####
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

p1 <- ggplot(data=data.frame(gamma=rgamma(10000,gamma_hyp1,gamma_hyp2)),
             aes(x=gamma,y=..density..))+
  geom_histogram(bins=30)+ggtitle("Prior on Gamma")+xlab("Gamma")+ylab("Density")
p2 <- ggplot(data=data.frame(K=rpois(10000,lambda)+1),aes(K,y=..density..))+
  geom_histogram(bins=30)+ggtitle("Prior on K")+xlab("K")+ylab("Density")
p3 <- ggplot(plot_Kplus,aes(`K+`,density))+geom_line()+geom_point()+
  facet_wrap(. ~gamma,labeller = label_parsed)+
  scale_x_continuous(breaks=1:10)+ylab("Prior Density")+ggtitle("Prior on K+ | gamma")
grid.arrange(p1,p2,p3,layout_matrix=rbind(c(1,2),c(3,3)))
ggsave("Results_Plots/AIBS_Prior.pdf",
       grid.arrange(p1,p2,p3,layout_matrix=rbind(c(1,2),c(3,3))),
       width=11,height=6)

#### Run MFM Analyses ####
# res <- btlb_mfm(Pi=Pi,X=X,M=M,Pi_full=Pi_full,
#                 lambda=lambda,a=a,b=b,gamma1=gamma1,gamma2=gamma2,
#                 gamma_hyp1=gamma_hyp1,gamma_hyp2=gamma_hyp2,
#                 mh_pjk = .05, mh_thetak = 5, mh_gamma = 1,
#                 startK=17,max_iters=1000,mh_iters=10,burn=.50,thin=5,
#                 seed = 1)
load("~/BTLBinomial/AIBS_Panel2_res.RData")

#### Key MFM Results ####
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
p4<-ggplot(plot_p,aes(x=factor(Var1),y=value,fill=factor(Var2),color=factor(Var2)))+
  geom_violin()+
  theme(legend.position="none")+ylim(c(0,1))+
  xlab("Proposal Quality Parameter")+ylab("Density")
plot_theta <- na.exclude(reshape2::melt(res$theta[res$Kplus %in% c(1,2),1:2]))
plot_theta$Var2 <- factor(plot_theta$Var2,levels=2:1,labels=1:2)
p5<-ggplot(plot_theta,aes(x=factor(Var2),y=value,fill=factor(Var2),color=factor(Var2)))+
  geom_violin()+ylim(c(0,40))+xlab(expression(theta))+ylab("Density")+
  labs(fill="Cluster",color="Cluster")
grid.arrange(p2,p1,p3,p4,p5,layout_matrix=rbind(c(1,2,3,3,3),c(4,4,4,4,5)))
ggsave("Results_Plots/AIBS_res1.pdf",
       grid.arrange(p2,p1,p3,p4,p5,layout_matrix=rbind(c(1,2,3,3,3),c(4,4,4,4,5))),
       width=11,height=8)

TopProposals <- plot_p %>% group_by(Var1,Var2) %>% summarize(mean(value))
names(TopProposals) <- c("Proposal","Cluster","Posterior Mean p")
SummaryTheta <- plot_theta %>% group_by(Var2) %>% summarize(mean(value))
names(SummaryTheta) <- c("Cluster","Posterior Mean theta")
TopProposals[,c(2,1,3)] %>% group_by(Cluster) %>% slice_min(order_by = `Posterior Mean p`, n = 4)
SummaryTheta

#### Diagnostics and Other Trace Plots ####
plot_iters <- round(seq(res$burn*res$max_iters*res$mh_iters,res$max_iters*res$mh_iters,length=1000))
accept_data <- melt(data.frame(iters=plot_iters,
                               accept_p=res$accept_p[round(seq(1,length(res$accept_p),length=1000))],
                               accept_theta=res$accept_theta[round(seq(1,length(res$accept_theta),length=1000))],
                               accept_gamma=res$accept_theta[round(seq(1,length(res$accept_gamma),length=1000))]),
                    id.vars="iters")
accept_data$Parameter <- factor(accept_data$variable,labels=c("p","theta","gamma"))
p1 <- ggplot(data=accept_data,aes(x=iters,y=value,color=Parameter))+
  geom_line()+ylim(c(0,1))+ylab("Acceptance Probability")+xlab("Iteration")+
  theme(legend.position = "bottom")+ggtitle("Acceptance Probabilities")
p2 <- ggplot(data=data.frame(iters=plot_iters,K=res$K[-1]),aes(iters,K))+geom_line()+
  ylim(c(1,max(res$K,res$Kplus)))+ylab("K")+xlab("Iteration")+
  ggtitle("Trace Plot: K")
p3 <- ggplot(data=data.frame(iters=plot_iters,Kplus=res$Kplus[-1]),aes(iters,Kplus))+geom_line()+
  ylim(c(1,max(res$K,res$Kplus)))+ylab("K+")+xlab("Iteration")+
  ggtitle("Trace Plot: K+")
p4 <- ggplot(data=data.frame(iters=plot_iters,gamma=res$gamma[-1]),aes(iters,gamma))+geom_line()+
  ylim(c(0,4))+ylab("gamma")+xlab("Iteration")+
  ggtitle("Trace Plot: gamma")
p5 <- ggplot(reshape2::melt(res$pi),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
  geom_line()+ylim(c(0,1))+theme(legend.position="right")+
  ylab("Class Proportion Estimate")+xlab("Iteration")+
  scale_x_continuous(breaks=seq(0,1000,length=6),labels=paste0(5:10,"k"))+
  labs(color="Class")+ggtitle("Trace Plot: Class Proportions, pi")+
  theme(legend.position = "none")
p6 <- ggplot(reshape2::melt(res$Z),aes(x=Var1,y=jitter(value,.5),group=Var2,color=factor(Var2)))+
  geom_line()+theme(legend.position="none")+ylim(c(0,max(res$Kplus)+1))+
  ylab("Class")+xlab("Iteration")+
  scale_x_continuous(breaks=seq(0,1000,length=6),labels=paste0(5:10,"k"))+
  labs(color="Class")+ggtitle("Trace Plot: Class Memberships, Z")
p7 <- ggplot(reshape2::melt(res$theta),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
  geom_line()+theme(legend.position="right")+ylim(c(0,max(res$theta)+1))+
  ylab("theta")+xlab("Iteration")+
  scale_x_continuous(breaks=seq(0,1000,length=6),labels=paste0(5:10,"k"))+
  labs(color="Class")+ggtitle("Trace Plot: theta")
plot_p <- reshape2::melt(res$p)
plot_p$jk <- as.factor(paste0(plot_p$Var1,"_",plot_p$Var2))
p8 <- ggplot(plot_p,aes(x=Var3,y=value,group=jk,color=factor(Var2)))+
  facet_wrap(vars(Var1))+
  geom_line()+theme(legend.position="bottom")+ylim(c(0,1))+
  ylab("p")+xlab("Iteration")+
  scale_x_continuous(breaks=seq(0,1000,length=6),labels=paste0(5:10,"k"))+
  labs(color="Class")+ggtitle("Trace Plot: p")

ggsave("Results_Plots/AIBS_trace1.pdf",grid.arrange(p1,p2,p3,nrow=1),
       width=11,height=4)
ggsave("Results_Plots/AIBS_trace2.pdf",grid.arrange(p4,p5,p6,p7,layout_matrix=rbind(c(2,4),c(1,3))),
       width=11,height=8)
ggsave("Results_Plots/AIBS_trace3.pdf",p8,width=11,height=11)


