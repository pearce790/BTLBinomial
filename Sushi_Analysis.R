#### Sushi Reproducibility ####
source("~/BTLBinomial/source.R")

#### Load Data and Choose Hyperparameters ####
load("~/BTLBinomial/Sushi.RData")

#### EDA ####
plotX <- melt(X)
names(plotX) <- c("ID","Sushi Type","Rating")
plotX$`Sushi Type` <- factor(plotX$`Sushi Type`,levels=paste0("Item",1:10),
                             labels=c("Shrimp","Sea Eel","Tuna","Squid","Sea Urchin",
                                      "Salmon Roe","Egg","Fatty Tuna","Tuna Roll","Cucumber Roll"))
plotX$`Sushi Type` <- factor(plotX$`Sushi Type`,levels=c("Fatty Tuna","Tuna","Shrimp","Tuna Roll","Sea Eel","Salmon Roe",
                                                         "Squid","Egg","Sea Urchin","Cucumber Roll"))
p1 <- ggplot(plotX,aes(x=`Sushi Type`,y=Rating))+geom_jitter(size=.2)+xlab(NULL)+
  scale_y_continuous(breaks=0:4,labels=paste0(" ",0:4))

plotPi <- melt(Pi)
names(plotPi) <- c("ID","Rank","Sushi Type")
plotPi$`Sushi Type` <- factor(plotPi$`Sushi Type`,levels=paste0(1:10),
                             labels=c("Shrimp","Sea Eel","Tuna","Squid","Sea Urchin",
                                      "Salmon Roe","Egg","Fatty Tuna","Tuna Roll","Cucumber Roll"))
plotPi$`Sushi Type` <- factor(plotPi$`Sushi Type`,levels=c("Fatty Tuna","Tuna","Shrimp","Tuna Roll","Sea Eel","Salmon Roe",
                                                         "Squid","Egg","Sea Urchin","Cucumber Roll"))
plotPi$Rank <- factor(plotPi$Rank,levels=paste0("Place",10:1),labels=rev(c("First","Second","Third","Fourth","Fifth",
                                                                       "Sixth","Seventh","Eighth","Ninth","Tenth")))
p2 <- ggplot(plotPi,aes(x=`Sushi Type`,value=Rank,color=Rank,fill=Rank))+geom_bar(stat="count")+
  theme(legend.position = "bottom")+ylab("Count")+
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#555555","#999999"),
                    breaks=c("First","Second","Third","Fourth","Fifth","Sixth","Seventh","Eighth","Ninth","Tenth"))+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#555555","#999999"),
                     breaks=c("First","Second","Third","Fourth","Fifth","Sixth","Seventh","Eighth","Ninth","Tenth"))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow=1))+
  scale_y_continuous(breaks=seq(0,5000,length=6),labels=c("0",paste0(1:5,"k")))
grid.arrange(p1,p2,nrow=2,heights=c(.45,.55))
ggsave("Results_Plots/Sushi_EDA.pdf",grid.arrange(p1,p2,nrow=2,heights=c(.45,.55)),
       width=11,height=8)




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
SummaryTopTypes <- TopTypes[,c(2,1,3)] %>% group_by(Cluster) %>% slice_min(order_by = `Posterior Mean p`, n = 3)

sum_table <- data.frame(Cluster=1:9,Pi=round(apply(res$pi,2,mean)[1:9],3),Top3=NA,Theta=round(apply(res$theta,2,mean)[1:9],2))
for(Cluster in 1:9){
  sum_table[Cluster,"Top3"] <- paste0(unlist(c(SummaryTopTypes[SummaryTopTypes$Cluster == Cluster,"Type"])),collapse=", ")
}

print(xtable::xtable(sum_table,caption="Sushi Analysis Results"),include.rownames=FALSE)




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

#### Diagnostics and Other Trace Plots ####

gof <- get_gof(posterior=res,post_samples=200,reps=5000,X=X,Pi=Pi,Pi_full=NULL,seed=1)
p1 <- gof$p1 + ylim(c(0,M))+xlab("Sushi Type")+
  scale_x_discrete(labels=c("Shrimp","Sea Eel","Tuna","Squid",
                            "Sea Urchin","Salmon Roe","Egg",
                            "Fatty Tuna","Tuna Roll","Cucumber Roll"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- gof$p2 + xlab("Sushi Type")+
  scale_x_discrete(labels=c("Shrimp","Sea Eel","Tuna","Squid",
                            "Sea Urchin","Salmon Roe","Egg",
                            "Fatty Tuna","Tuna Roll","Cucumber Roll"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
ggsave("Results_Plots/Sushi_gof.pdf",grid.arrange(p1,p2,gof$p3,nrow=1),
       width=11,height=4)

plot_iters <- round(seq(res$burn*res$max_iters*res$mh_iters,res$max_iters*res$mh_iters,length=1000))
accept_data <- melt(data.frame(iters=plot_iters,
                               accept_p=res$accept_p[round(seq(1,length(res$accept_p),length=1000))],
                               accept_theta=res$accept_theta[round(seq(1,length(res$accept_theta),length=1000))],
                               accept_gamma=res$accept_theta[round(seq(1,length(res$accept_gamma),length=1000))]),
                    id.vars="iters")
accept_data$Parameter <- factor(accept_data$variable,labels=c("p","theta","gamma"))
# p1 <- ggplot(data=accept_data,aes(x=iters,y=value,color=Parameter))+
#   geom_line()+ylim(c(0,1))+ylab("Acceptance Probability")+xlab("Iteration")+
#   theme(legend.position = "bottom")+ggtitle("Acceptance Probabilities")
p2 <- ggplot(data=data.frame(iters=plot_iters,K=res$K[-1]),aes(iters,K))+geom_line()+
  ylab("K")+xlab("Iteration")+
  scale_x_continuous(breaks=seq(25000,50000,length=6),labels=paste0(seq(25,50,length=6),"k"))+
  scale_y_continuous(breaks=1:10,limits=c(1,max(res$K,res$Kplus)))+
  ggtitle("Trace Plot: K")
p3 <- ggplot(data=data.frame(iters=plot_iters,Kplus=res$Kplus[-1]),aes(iters,Kplus))+geom_line()+
  ylab("K+")+xlab("Iteration")+
  scale_x_continuous(breaks=seq(25000,50000,length=6),labels=paste0(seq(25,50,length=6),"k"))+
  scale_y_continuous(breaks=1:10,limits=c(1,max(res$K,res$Kplus)))+
  ggtitle("Trace Plot: K+")
p4 <- ggplot(data=data.frame(iters=plot_iters,gamma=res$gamma[-1]),aes(iters,gamma))+geom_line()+
  ylim(c(0,5))+ylab(expression(gamma))+xlab("Iteration")+
  scale_x_continuous(breaks=seq(25000,50000,length=6),labels=paste0(seq(25,50,length=6),"k"))+
  ggtitle(expression("Trace Plot: " ~gamma ))
p5 <- ggplot(reshape2::melt(res$pi),aes(x=Var1,y=value,group=Var2,color=factor(Var2)))+
  geom_line()+ylim(c(0,0.3))+theme(legend.position="right")+
  ylab(expression(pi))+xlab("Iteration")+
  scale_x_continuous(breaks=seq(0,5000,length=6),labels=paste0(seq(25,50,length=6),"k"))+
  labs(color="Class")+ggtitle(expression("Trace Plot: Class Proportions, "~pi))+
  theme(legend.position = "none")

ggsave("Results_Plots/Sushi_trace1.pdf",grid.arrange(p2,p3,p4,nrow=1),
       width=11,height=4)
ggsave("Results_Plots/Sushi_trace2.pdf",p5,
       width=11,height=4)




