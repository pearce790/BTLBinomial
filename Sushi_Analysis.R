#### Sushi Reproducibility ####
source("~/BTLBinomial/source.R")

#### Load Data and Choose Hyperparameters ####
load("~/BTLBinomial/Sushi.RData")


#### Explore Influence of Hyperparameters on K and K+ ####
plot_Kplus <- matrix(NA,nrow=0,ncol=3)
for(gamma in c(0.01,0.5,1,2,3,4)){
  pmfstatic2 <- nClusters(Kplus=1:20,N=nrow(X),type="static",gamma=gamma,maxK=50)
  dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = lambda))
  plot_Kplus <- rbind(plot_Kplus,matrix(c(rep(gamma,20),1:20,dens),ncol=3))
}
plot_Kplus <- as.data.frame(plot_Kplus)
names(plot_Kplus) <- c("gamma","K+","density")
plot_Kplus$gamma <- factor(plot_Kplus$gamma)
levels(plot_Kplus$gamma) <- c(expression("gamma: 0.01"),expression("gamma: 0.5"),expression("gamma: 1"),
                              expression("gamma: 2"),expression("gamma: 3"),expression("gamma: 4"))

p1 <- ggplot(data=data.frame(gamma=rgamma(10000,gamma_hyp1,gamma_hyp2)),
             aes(x=gamma,y=..density..))+
  geom_histogram(bins=30)+ggtitle("Prior on Gamma")+xlab("Gamma")+ylab("Density")
p2 <- ggplot(data=data.frame(K=rpois(1000000,lambda)+1),aes(K,y=..density..))+
  geom_histogram(bins=100)+ggtitle("Prior on K")+xlab("K")+ylab("Density")
p3 <- ggplot(plot_Kplus,aes(`K+`,density))+geom_line()+geom_point()+
  facet_wrap(. ~gamma,labeller = label_parsed)+
  scale_x_continuous(breaks=1:20)+ylab("Prior Density")+ggtitle("Prior on K+ | gamma")
grid.arrange(p1,p2,p3,layout_matrix=rbind(c(1,2),c(3,3)))
ggsave("Results_Plots/Sushi_Prior.pdf",
       grid.arrange(p1,p2,p3,layout_matrix=rbind(c(1,2),c(3,3))),
       width=11,height=6)

#### Run MFM Analyses ####
# posterior <- btlb_mfm(Pi=Pi,X=X,M=M,lambda=lambda,a=a,b=b,
#                       gamma1=gamma1,gamma2=gamma2,gamma_hyp1=gamma_hyp1,gamma_hyp2=gamma_hyp2,
#                       Pi_full=NULL,mh_pjk = .2, mh_thetak = 5,mh_gamma=0.5,
#                       startK=1,max_iters=5000,mh_iters=10,burn=0.5,thin=5,seed=1)
load("~/BTLBinomial/Sushi_res.RData")


#### Key MFM Results ####
p1<-ggplot(data=data.frame(Kplus=res$Kplus),aes(Kplus,y=..density..))+
  geom_histogram(breaks=c(8.5,9.5))+xlim(c(7.5,10.5))
p2<-ggplot(data=data.frame(gamma=res$gamma),aes(x=gamma,y=..density..))+
  geom_histogram(bins=25)+xlim(c(0,5))+ylab("Density")+
  xlab(expression(gamma))
plot_pi <- melt(res$pi[,1:9])
names(plot_pi) <- c("It","Cluster","pi")
p3<-ggplot(plot_pi,aes(factor(Cluster),pi))+geom_violin()+
  ylab(expression(pi))+xlab("Cluster")+ylim(c(0,.3))
grid.arrange(p1,p2,p3,layout_matrix=rbind(c(1,2),c(3,3)))
ggsave("Results_Plots/Sushi_res1.pdf",
       grid.arrange(p1,p2,p3,layout_matrix=rbind(c(1,2),c(3,3))),
       width=11,height=8)

plot_p <- reshape2::melt(res$p[,1:9,])
plot_p$Var1 <- factor(plot_p$Var1,levels=1:10,
                      labels=c("Shrimp","Sea Eel","Tuna","Squid","Sea Urchin",
                               "Salmon Roe","Egg","Fatty Tuna","Tuna Roll","Cucumber Roll"))
plot_p$Var2 <- factor(plot_p$Var2)
plot_theta <- reshape2::melt(res$theta[,1:9])
plot_theta$Var2 <- factor(plot_theta$Var2)
TopTypes <- plot_p %>% group_by(Var1,Var2) %>% summarize(mean(value))
names(TopTypes) <- c("Type","Cluster","Posterior Mean p")
SummaryTheta <- plot_theta %>% group_by(Var2) %>% summarize(mean(value))
names(SummaryTheta) <- c("Cluster","Posterior Mean theta")
print(n=27,TopTypes[,c(2,1,3)] %>% group_by(Cluster) %>% slice_min(order_by = `Posterior Mean p`, n = 3))
SummaryTheta


Zhat <- salso(x=res$Z[round(seq(1,5000,length=1000)),seq(1,5000,length=500)],
              loss="binder",maxNClusters=max(res$Z))
sumZhat <- summary(Zhat)
plotZhat <- melt(sumZhat$psm[sumZhat$order,rev(sumZhat$order)])
p4 <- ggplot(plotZhat,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  scale_x_continuous(breaks=NULL,expand = c(0, 0))+
  scale_y_continuous(breaks=NULL,expand = c(0, 0))+
  labs(fill="Cluster\nSimilarity")+theme_bw()+
  xlab(NULL)+ylab(NULL)
ggsave("Results_Plots/Sushi_res2.pdf",p4,width=11,height=10)









