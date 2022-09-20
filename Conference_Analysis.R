#### Sushi Reproducibility ####
source("~/BTLBinomial/source.R")

#### Choose Hyperparameters and Run Simulation ####
I <- 50
J <- 50

# results <- matrix(NA,nrow=0,ncol=4+2*50+1+50)
# for(iter in 1:100){for(M in c(4,9)){for(R in c(4,8,12,24)){for(theta in c(1,5,10,20,40)){
#   #initialize
#   set.seed(iter*M*R*theta)
#   print(c(M,R,theta,iter))
#   X <- matrix(NA,nrow=50,ncol=50)
#   Pi <- matrix(NA,nrow=50,ncol=4)
#   Pi_full <- matrix(NA,nrow=50,ncol=R)
#   
#   #make assignments at random
#   continue <- FALSE
#   while(continue==FALSE){
#     assignments <- matrix(0,nrow=50,ncol=50)
#     tryCatch({
#       for(i in 1:I){
#         assignments[i,sample(which(apply(assignments,2,sum)<R),R)] <- 1
#       }
#     },error=function(e){})
#     if(all(apply(assignments,2,sum)==R)){continue <- TRUE}
#   }
#   
#   #draw scores
#   p <- runif(50,min=0,max=1)
#   for(i in 1:50){
#     which_proposals <- which(assignments[i,]==1)
#     X[i,which_proposals] <- rbinom(n=R,size=M,prob=p[which_proposals])
#     Pi[i,1:4] <- sample(which_proposals,4,prob=exp(-theta*p[which_proposals]))
#     Pi_full[i,1:4] <- Pi[i,1:4]
#     if(R>4){Pi_full[i,5:R] <- setdiff(which_proposals,Pi[i,1:4])}
#   }
#   
#   tryCatch({
#     map_est <- btlb_map(Pi=Pi,X=X,M=M,Pi_full=Pi_full,K=1,gamma=1,a=1,b=1,gamma1=5,gamma2=0.25,tol=0.001,maxit=100,verbose=FALSE)
#     results <- rbind(results,c(M,R,theta,iter,p,map_est$p,map_est$theta,apply(X,2,function(x){mean(x,na.rm=TRUE)})))
#   },error=function(e){
#     print("Fail MAP Search for M,R,theta,iter:")
#     print(c(M,R,theta,iter))
#   })
# }}}
#   
# }
# results <- as.data.frame(results)
# names(results) <- c("M","R","theta","Iter",paste0("p",1:50),paste0("pmap",1:50),
#                     "thetamap",paste0("Xbar",1:50))

load("~/BTLBinomial/Conference_Analysis.RData")



#### Analyze Results ####

accuracy_p <- cbind(melt(results[,c(1,2,3,5:54)],id.vars = c(1:3))[,c(1,2,5)],
                    melt(results[,c(1,2,3,55:104)],id.vars = c(1:3))[,5])
names(accuracy_p) <- c("M","R","p","phat")
accuracy_theta <- results[,c(1,2,3,105)]
accuracy_theta$theta <- factor(accuracy_theta$theta)
accuracy_thetaplot <- data.frame(theta=c(1,5,10,20,40),thetamap=c(1,5,10,20,40))
accuracy_thetaplot$theta <- factor(accuracy_thetaplot$theta)
p1 <- ggplot(accuracy_p,aes(p,phat))+geom_point(size=.01,alpha=.1)+
  facet_grid(cols=vars(M),rows=vars(R),labeller = label_both)+
  geom_abline(slope=1,intercept=0,color="red",lty=2)+
  labs(x=expression(p[0]),y=expression(hat(p)))
p2 <- ggplot(accuracy_theta,aes(theta,thetamap))+geom_violin()+
  geom_point(data=accuracy_thetaplot,aes(theta,thetamap),color="red")+
  facet_grid(cols=vars(M),rows=vars(R),labeller = label_both)+
  labs(x=expression(theta[0]),y=expression(hat(theta)))
ggsave("Results_Plots/Conference_res1.pdf",grid.arrange(p1,p2,ncol=2),
       width=9,height=7)

diff_p <- cbind(results[,1:4],results[,55:104] - results[,5:54])
diff_p_long <- melt(diff_p,id.vars=c("M","R","theta","Iter"))
diff_theta <- cbind(results[,1:3],results[,105]-results[,3])
names(diff_p)[5:54] <- paste0("p",1:50)
names(diff_theta)[4] <- "diff_theta"


p3 <- ggplot(diff_p_long,aes(x=factor(theta),y=value))+geom_violin()+
  facet_grid(rows=vars(R),cols=vars(M),labeller = label_both)+
  labs(x=expression(theta[0]),y=expression(hat(p) - p[0]))+
  geom_abline(slope=0,intercept=0,color="red",lty=2)
p4 <- ggplot(diff_theta,aes(factor(theta),y=diff_theta))+geom_violin()+
  facet_grid(rows=vars(R),cols=vars(M),labeller=label_both)+
  geom_abline(slope=0,intercept=0,color="red",lty=2)+
  labs(x=expression(theta[0]),y=expression(hat(theta) - theta[0]))+
  theme(legend.position = "none")
ggsave("Results_Plots/Conference_res2.pdf",grid.arrange(p3,p4,ncol=2),
       width=9,height=4)
rm(diff_p,diff_p_long,diff_theta)


btlb_distance <- apply(results,1,function(res){
  rankrate::kendall(order(rank(res[55:104],ties.method="random")),order(res[5:54]))
})
standard_distance <- apply(results,1,function(res){
  rankrate::kendall(order(rank(res[106:155],ties.method="random")),order(res[5:54]))
})

diff_pi <- cbind(results[,1:4],btlb_distance,standard_distance)
diff_pi$btlb_distance <- diff_pi$btlb_distance/choose(50,2)
diff_pi$standard_distance <- diff_pi$standard_distance/choose(50,2)
p5 <- ggplot(diff_pi,aes(standard_distance,btlb_distance,color=factor(theta)))+
  geom_point(size=1,alpha=.3)+
  facet_grid(cols=vars(R),rows=vars(M),labeller=label_both)+
  geom_abline(slope=1,intercept=0,color="red",linetype=2)+
  labs(x=expression(paste("% Incorrect Pairs Between ",hat(pi)[0]^X," and ",pi[0])),
       y=expression(paste("% Incorrect Pairs Between ",hat(pi)[0]^BTLB," and ",pi[0])),
       color=expression(theta[0]))+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes=list(alpha=1,size=2)))
ggsave("Results_Plots/Conference_res3.pdf",p5,
       width=9,height=5)
