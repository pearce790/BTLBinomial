#### Packages ####
library(gtools)
library(matrixStats)
library(ggplot2)
library(fipp) #for calculating prior on Kplus
library(gridExtra)
library(dplyr)
library(salso)
library(reshape2)

#### BTL-B Density and Generation ####
dbtlb <- function(Pi,X,p,theta,M,log=FALSE,Pi_full=NULL){
  if(is.vector(Pi)){Pi <- matrix(Pi,nrow=1)}
  if(is.vector(X)){X <- matrix(X,nrow=1)}
  if(is.vector(Pi_full)){Pi_full <- matrix(Pi_full,nrow=1)}
  if(!is.matrix(X)){stop("X must be a matrix of ratings")}
  if(!is.matrix(Pi)){stop("Pi must be a matrix of rankings")}
  if(length(p)!=ncol(X)){stop("length(p) must equal ncol(X)")}
  
  worth <- exp(-theta*p)
  logd <- dbtl(Pi,worth,log=T,Pi_full=Pi_full)+sum(apply(X,1,function(x){dbinom(x,M,p,log=T)}),na.rm=T)
  
  if(log){
    return(logd)
  }else{
    return(exp(logd))
  }
}
dbtl <- function(Pi,worth,log=FALSE,Pi_full=NULL){
  if(is.vector(Pi)){Pi <- matrix(Pi,nrow=1)}
  if(is.vector(Pi_full)){Pi_full <- matrix(Pi_full,nrow=1)}
  if(!is.matrix(Pi)){stop("Pi must be a matrix of (partial) rankings")}
  if(!is.null(Pi_full)){if(nrow(Pi)!=nrow(Pi_full)){stop("nrow(Pi) must equal nrow(Pi_full).")}}
  if(length(worth)<max(c(0,Pi,Pi_full),na.rm=T)){stop("Ensure one worth parameter for every object in Pi and/or Pi_full")}
  
  if(sum(apply(Pi,1,function(pi){!all(is.na(pi))}))==0){ #if all rankings empty, return likelihood of 1 (convenience)
    if(log){return(0)}else{return(1)}
  }
  
  if(!is.null(Pi_full)){
    logd <- sum(unlist(lapply(1:nrow(Pi),function(i){
      pi <- na.exclude(Pi[i,])
      if(length(pi)==0){return(0)
      }else{
        num <- worth[pi]
        denom_terms <- na.exclude(worth[Pi_full[i,]])
        return(sum(log(num))-sum(unlist(lapply(1:length(num),function(j){log(sum(denom_terms[j:length(denom_terms)]))}))))
      }
    })))
  }else{
    logd <- sum(apply(Pi,1,function(pi){
      if(length(na.exclude(pi))==0){return(0)
      }else{
        pi <- na.exclude(pi)
        num <- worth[pi]
        denom_terms <- c(worth[pi],worth[setdiff(1:length(worth),pi)])
        return(sum(log(num)-log(rev(cumsum(rev(denom_terms)))[1:length(pi)])))
      }
    }))
  }
  
  if(log){return(logd)
  }else{return(exp(logd))}
}
rbtlb <- function(I,p,theta,M,R=length(p)){
  J <- length(p)
  
  X <- as.matrix(t(replicate(I,rbinom(n=J,size=M,prob=p))))
  Pi <- as.matrix(t(replicate(I,{sample(1:J,R,prob=exp(-theta*p))})))
  
  return(list(X=X,Pi=Pi))
}

#### MAP Estimation ####
set.seed(1)
J <- 10
M <- 4
dat1 <- rbtlb(I=20,p=runif(J),theta=10,M=M)
dat2 <- rbtlb(I=10,p=runif(J),theta=20,M=M)
X <- rbind(dat1$X,dat2$X)
Pi <- rbind(dat1$Pi,dat2$Pi)
rm(dat1,dat2,J)
Pi_full <- NULL
K <- 3
gamma <- 1
a <- 1
b <- 1
gamma1 <- 10
gamma2 <- 0.5

map_btlb <- function(Pi,X,M,Pi_full=NULL,K,gamma,
                     a,b,gamma1,gamma2,tol=1,maxit=50,verbose=TRUE,seed=NULL){
  
  
  # function to calculate objective function
  get_obj <- function(ptheta,pi){
    lik_terms <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      pk <- ptheta[1:J,k]
      thetak <- ptheta[J+1,k]
      worthk <- exp(-thetak*pk)
      lik_terms[,k] <- log(pi[k])+unlist(lapply(1:I,function(i){
        dbtl(Pi=Pi[i,],worth=worthk,Pi_full=Pi_full[i,],log=T)+sum(dbinom(x=X[i,],size=M,prob=pk,log=T),na.rm=T)
      }))
    }
    lik_term <- sum(apply(lik_terms,1,logSumExp))
    prior_term <- log(ddirichlet(pi,rep(gamma,K)))+sum(dbeta(ptheta[1:J,],a,b,log=TRUE))+
      sum(dgamma(ptheta[J+1,],gamma1,gamma2,log=TRUE))
    currobj <- lik_term+prior_term
    return(currobj)
  }
  
  # model constants
  I <- nrow(X)
  J <- ncol(X)
  
  # initializers
  if(!is.null(seed)){set.seed(seed)}
  pi <- rep(1/K,K)
  ptheta <- matrix(c(rbeta(J*K,a,b),rgamma(K,gamma1,gamma2)),byrow=TRUE,nrow=J+1,ncol=K)
  Z <- matrix(NA,nrow=I,ncol=K)
  
  #get current loglikelihood
  curr_obj <- get_obj(ptheta,pi)
  
  currdiff <- Inf
  iter <- 0
  while(currdiff>tol & iter<=maxit){
    iter <- iter + 1
    
    # E-Step: Update Z
    t1 <- Sys.time()
    for(i in 1:I){
      lprobs <- unlist(lapply(1:K,function(k){
        log(pi[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=ptheta[1:J,k],
                         theta=ptheta[J+1,k],M=M,log=T)
      }))
      Z[i,] <- exp(lprobs-logSumExp(lprobs))
    }

    
    # M-Step:
    ## Update (p,theta)
    t2 <- Sys.time()
    for(k in 1:K){
      #print(k)
      zhat_k <- Z[,k]
      if(length(which_vals)==0){which_vals <- 1:I}
      obj <- function(par){
        worthk <- exp(-par[1:J]*par[J+1])
        lik_term <- sum(unlist(lapply(which_vals,function(i){
          zhat_k[i]*(dbtl(Pi=Pi[i,],worth=worthk,log=T)+sum(dbinom(X[i,],M,par[1:J],log=T)))
        })))
        prior_term <- sum(dbeta(par[1:J],a,b,log=TRUE))+dgamma(par[J+1],gamma1,gamma2,log=TRUE)
        return(-lik_term-prior_term)
      }
      ptheta[,k] <- optim(ptheta[,k],obj,method="L-BFGS-B",lower=rep(0.00001,J+1),upper=c(rep(1-.00001,J),qgamma(.9999,gamma1,gamma2)))$par
    }
    
    ## Update pi
    t3 <- Sys.time()
    pi <- (gamma-1+apply(Z,2,sum))/(I - K + gamma*K)
    if(round(sum(pi),5)!=1){print("Error! Something with pi is wrong")}

   
    # Iteration Updates
    t4 <- Sys.time()
    new_obj <- get_obj(ptheta,pi)
    currdiff <- new_obj-curr_obj
    curr_obj <- new_obj
    t5 <- Sys.time()
    
    if(verbose){
      print(paste0("End of Iteration ",iter))
      print("Time to update Z, ptheta,pi,objective: ")
      print(round(c(t2-t1,t3-t2,t4-t3,t5-t4),2))
      print(paste0("Difference: ",round(currdiff,3)))
    }
  }
  
  return(list(p=ptheta[1:J,],theta=ptheta[J+1,],pihat=pi,Z=Z,obj=curr_obj))
  
}

map_btlb(Pi,X,M,K=2,gamma=gamma,a=a,b=b,gamma1=gamma1,gamma2=gamma2,seed=2)


#### MAP Sandbox ####
#set.seed(1)
J <- 10
M <- 4
dat1 <- rbtlb(I=20,p=runif(J),theta=10,M=M)
dat2 <- rbtlb(I=10,p=runif(J),theta=20,M=M)
X <- rbind(dat1$X,dat2$X)
Pi <- rbind(dat1$Pi,dat2$Pi)
rm(dat1,dat2,J)
Pi_full <- Pi
gamma <- 1
a <- 1
b <- 1
gamma1 <- 10
gamma2 <- 0.5
K <- 2
tol <- 1
maxit <- 50

map_btlb(2,Pi,X,M,Pi_full=NULL,gamma,a,b,gamma1,gamma2,tol=.01,maxit=50)$obj


#### MFMM Estimation ####
btlb_mfm <- function(Pi,X,M,Pi_full=NULL,lambda,a,b,gamma1,gamma2,gamma_hyp1,gamma_hyp2,
                     startK = 1, mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                     max_iters = 100, mh_iters = 5, burn = 0.5,thin = 2,seed = NULL){
  
  print("Initializing Chain")
  if(!is.null(seed)){set.seed(seed)}
  
  ## Determine Constants
  I <- nrow(X)
  J <- ncol(X)
  R <- ncol(Pi)
  
  ## Set Up Data Storage
  which_keep <- seq(from=round(burn*max_iters*mh_iters),to=max_iters * mh_iters,by=thin)
  K_all <- c()
  Kplus_all <- c()
  gamma_all <- c()
  pi_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  p_all <- array(NA,dim=c(J,I,length(which_keep)))
  theta_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  Z_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  accept_p <- c()
  accept_theta <- c()
  accept_gamma <- c()
  
  ## Initialize at Random
  K <- startK
  gamma <- rgamma(1,gamma_hyp1,gamma_hyp2)
  pi <- c(rdirichlet(1,rep(gamma,K)))
  ptheta <- matrix(NA,nrow=J+1,ncol=K)
  for(k in 1:K){ptheta[,k] <- c(rbeta(J,a,b),rgamma(1,gamma1,gamma2))}
  Z <- rep(NA,I)
  
  print("Starting Chain")
  iter <- 1
  progress_iters <- round(seq(0,max(which_keep),length=11))[-1]
  while(iter <= max(which_keep)){
    
    ## Step 1: Update Labels
    probs <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      pk <- ptheta[1:J,k]
      thetak <- ptheta[J+1,k]
      worthk <- exp(-thetak*pk)
      probs[,k] <- exp(log(pi[k])+unlist(lapply(1:I,function(i){
        dbtl(Pi=Pi[i,],worth=worthk,Pi_full=Pi_full[i,],log=T)+sum(dbinom(x=X[i,],size=M,prob=pk,log=T),na.rm=T)
      })))
    }
    Z <- apply(probs,1,function(prob){sample.int(K,1,prob=prob)})
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)})) #update Nk, Kplus, and reorder
    Kplus <- sum(Nk>0)
    ptheta <- as.matrix(ptheta[,which(Nk>0)])
    Ztmp <- Z
    for(k in 1:K){Ztmp[which(Z==k)] <- which(c(which(Nk>0),setdiff(1:K,which(Nk>0)))==k)}
    Z <- Ztmp
    rm(Ztmp)
    
    
    ## Step 2: Update non-empty component parameters
    for(mh_iter in 0:(mh_iters-1)){
      for(k in 1:Kplus){
        whichk <- which(Z == k)
        constant1 <- a+apply(matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        constant2 <- b+apply(M-matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        Pi_mat <- matrix(Pi[whichk,],nrow=length(whichk))
        if(is.null(Pi_full)){Pi_full_mat <- NULL}else{Pi_full_mat <- matrix(Pi_full[whichk,],nrow=length(whichk))}
        
        for(j in 1:J){ #update each p_jk
          prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)
          while(prop_pjk<=0 | prop_pjk>=1){prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)}
          prop_p <- ptheta[1:J,k]
          prop_p[j] <- prop_pjk
          curr_worth <- exp(-ptheta[J+1,k]*ptheta[1:J,k])
          prop_worth <- exp(-ptheta[J+1,k]*prop_p)
          
          logprob_prop <- dbtl(Pi=Pi_mat,worth=prop_worth,log=T,Pi_full=Pi_full_mat)+
            (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
          logprob_curr <- dbtl(Pi=Pi_mat,worth=curr_worth,log=T,Pi_full=Pi_full_mat)+
            (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          
          u <- runif(1)
          if(log(u) <  logprob_prop-logprob_curr){
            accept_p <- c(accept_p,1)
            ptheta[j,k] <- prop_pjk
          }else{accept_p <- c(accept_p,0)}
        }
        #update theta_k
        prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)
        while(prop_thetak<=0){prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)}
        prop_worth <- exp(-prop_thetak*ptheta[1:J,k])
        curr_worth <- exp(-ptheta[J+1,k]*ptheta[1:J,k])
        
        logprob_prop <- dbtl(Pi_mat,worth=prop_worth,log=T,Pi_full=Pi_full_mat)+
          (gamma1-1)*log(prop_thetak)-gamma2*prop_thetak
        logprob_curr <- dbtl(Pi_mat,worth=curr_worth,log=T,Pi_full=Pi_full_mat)+
          (gamma1-1)*log(ptheta[J+1,k])-gamma2*ptheta[J+1,k]
        
        u <- runif(1)
        if(log(u) <  logprob_prop-logprob_curr){
          accept_theta <- c(accept_theta,1)
          ptheta[J+1,k] <- prop_thetak
        }else{accept_theta <- c(accept_theta,0)}
      }
      if(iter+mh_iter %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        p_all[,1:Kplus,it] <- ptheta[1:J,] #update only the first Kplus parameters,
        theta_all[it,1:Kplus] <- ptheta[J+1,] #(will update the rest later!)
      }
    }
    
    
    ## Step 3: Update K and gamma
    
    #update K
    probs <- unlist(lapply(Kplus:(Kplus+100),function(k){
      dpois(k-1,lambda,log=T)+lfactorial(k)-lfactorial(k-Kplus)+lgamma(gamma*k)-lgamma(I+gamma*k)
    }))
    K <- sample(Kplus:(Kplus+100),1,prob = exp(probs-logSumExp(probs)))
    
    #update gamma
    prop_gamma <- rnorm(1,gamma,mh_gamma)
    while(prop_gamma<=0){prop_gamma <- rnorm(1,gamma,mh_gamma)}
    logprob_prop <- dgamma(prop_gamma,gamma_hyp1,gamma_hyp2,log=TRUE)+lgamma(prop_gamma*K)-lgamma(I+prop_gamma*K)+
      sum(lgamma(Nk[1:Kplus]+prop_gamma)-lgamma(prop_gamma))
    logprob_curr <- dgamma(gamma,3,2,log=TRUE)+lgamma(gamma*K)-lgamma(I+gamma*K)+
      sum(lgamma(Nk[1:Kplus]+gamma)-lgamma(gamma))
    u <- runif(1)
    if(log(u) <  logprob_prop-logprob_curr){
      accept_gamma <- c(accept_gamma,1)
      gamma <- prop_gamma
    }else{accept_gamma <- c(accept_gamma,0)}
    
    ## Step 4: Update empty components
    if(K>Kplus){
      for(k in (Kplus+1):K){ptheta <- cbind(ptheta,c(rbeta(J,a,b),rgamma(1,gamma1,gamma2)))}
    }
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)}))
    pi <- c(rdirichlet(1,gamma+Nk))
    
    if(length(pi)!=K | length(Nk)!=K | ncol(ptheta)!=K ){stop("something wrong!")}
    if(Kplus > K){stop("something wrong with Kplus!")}
    
    ## Save Values and Update Counter
    for(mh_iter in 0:(mh_iters-1)){
      if((iter+mh_iter) %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        K_all <- c(K_all,K)
        Kplus_all <- c(Kplus_all,Kplus)
        gamma_all <- c(gamma_all,gamma)
        pi_all[it,1:K] <- pi
        Z_all[it,] <- Z
        if(K>Kplus){
          p_all[,(Kplus+1):K,it] <- ptheta[1:J,(Kplus+1):K] #update only the first Kplus parameters,
          theta_all[it,(Kplus+1):K] <- ptheta[J+1,(Kplus+1):K] #(will update the rest later!)
        }
      }}
    
    iter <- iter + mh_iters
    if(length(progress_iters)>0 & iter >= progress_iters[1]){
      print(paste0((11-length(progress_iters))*10,"% Complete: Iteration ",iter-1," out of ",max(which_keep)))
      print(paste0("Current K = ",K,"; Kplus = ",Kplus))
      progress_iters <- progress_iters[-1]
    }
  }
  print(paste0("Done! Saving ",length(which_keep)," estimate iterations after burning/thinning"))
  return(list(Z=Z_all,p=p_all[,1:max(K_all),],theta=theta_all[,1:max(K_all)],
              pi=pi_all[,1:max(K_all)],K=K_all,Kplus=Kplus_all,gamma=gamma_all,
              accept_p=(cumsum(accept_p)/1:length(accept_p))[seq(round(length(accept_p)*burn),length(accept_p),by=thin)],
              accept_theta=(cumsum(accept_theta)/1:length(accept_theta))[seq(round(length(accept_theta)*burn),length(accept_theta),by=thin)],
              accept_gamma=(cumsum(accept_gamma)/1:length(accept_gamma))[seq(round(length(accept_gamma)*burn),length(accept_gamma),by=thin)],
              max_iters=max_iters,mh_iters=mh_iters,burn=burn,thin=thin,seed=seed))
}

#### Prior Distribution on K+ Example ####

pmfstatic2 <- nClusters(Kplus=1:20,N=10,type="static",gamma=5,maxK=50)
dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = 7))
plot(dens)
sum((1:length(dens))*dens)
#### MFMM Sandbox ####

set.seed(1)
J <- 10
M <- 4
dat1 <- rbtlb(I=20,p=runif(J),theta=10,M=M)
dat2 <- rbtlb(I=10,p=runif(J),theta=20,M=M)
X <- rbind(dat1$X,dat2$X)
Pi <- rbind(dat1$Pi,dat2$Pi)
rm(dat1,dat2)
I <- nrow(X)
Pi_full <- NULL
a <- 1
b <- 1
gamma1 <- 10
gamma2 <- 0.5
lambda <- 3

res <- btlb_mfm(Pi=Pi,X=X,M=M,lambda=lambda,a=a,b=b,
                gamma1=gamma1,gamma2=gamma2,gamma_hyp1=3,gamma_hyp2=2,
                Pi_full=NULL,mh_pjk = .05, mh_thetak = 5,mh_gamma=0.5,
                startK=15,max_iters=100,mh_iters=10,burn=0,thin=2)
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
ggplot(plot_p,aes(x=Var3,y=value,group=jk,color=factor(jk)))+
  geom_line()+theme(legend.position="none")+ylim(c(0,1))+
  ylab("p")+xlab("Iteration (after burn/thin)")+
  ggtitle("Trace Plot: p")




#### AIBS Example ####

load("/Users/pearce790/Desktop/Elena Errant Files/Experiments/Panel2_2021.RData")
coi <- coi2
X <- as.matrix(X2)
Pi <- as.matrix(Pi2)
M <- 40
J <- ncol(X)
Pi_full <- matrix(NA,nrow=nrow(X),ncol=J)
for(i in 1:nrow(X)){
  tmp <- na.exclude(c(Pi[i,],setdiff(1:J,c(Pi[i,],which(coi[i,]==1)))))
  Pi_full[i,1:length(tmp)] <- tmp
}
fit <- fitdistrplus::fitdist(data=as.vector(na.exclude(c(X)))/40,distr="beta",method="mme")
a <- as.numeric(fit$estimate[1])
b <- as.numeric(fit$estimate[2])
gamma1 <- 10
gamma2 <- 0.5
rm(coi,coi2,fit,Pi2,X2,i,J,tmp,Proposals2,Reviewers2)

gamma_hyp1 <- 3
gamma_hyp2 <- 2
lambda <- 1


#hist(rgamma(1000,gamma_hyp1,gamma_hyp2),main="Prior Histogram on Gamma",xlab="gamma")
#hist(rpois(1000,lambda)+1,main="Prior Histogram on K",xlab="K")
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







