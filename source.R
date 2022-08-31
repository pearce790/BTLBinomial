#### Packages ####
library(gtools)
library(matrixStats)
library(ggplot2)

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
map_btlb <- function(K,Pi,X,M,a,b,gamma1,gamma2,tol=1,maxit=50){
  I <- nrow(X)
  J <- ncol(X)
  
  ptheta <- matrix(NA,nrow=J+1,ncol=K)
  ptheta[1:J,] <- rbeta(J*K,a,b)
  ptheta[J+1,] <- rgamma(K,gamma1,gamma2)
  pi <- rdirichlet(1,rep(gamma,K))
  
  lik_term <- sum(unlist(lapply(1:I,function(i){
    logSumExp(unlist(lapply(1:K,function(k){
      log(pi[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=ptheta[1:J,k],theta=ptheta[J+1,k],M=M,log=TRUE)
    })))
  })))
  prior_term <- log(ddirichlet(pi,rep(gamma,K)))+sum(dbeta(ptheta[1:J,],a,b,log=TRUE))+
    sum(dgamma(ptheta[J+1,],gamma1,gamma2,log=TRUE))
  currobj <- lik_term+prior_term
  
  Zhat <- matrix(NA,nrow=I,ncol=K)
  
  diff <- Inf
  nsteps <- 1
  while(diff > tol & nsteps <= maxit){
    print(paste0("EM Iteration: ",nsteps))
    
    # E Step
    for(i in 1:I){
      lprobs <- unlist(lapply(1:K,function(k){
        log(pi[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=ptheta[1:J,k],
                         theta=ptheta[J+1,k],M=M,log=T)
      }))
      Zhat[i,] <- exp(lprobs-logSumExp(lprobs))
    }
    
    # M Step
    pi <- (gamma-1+apply(Zhat,2,sum))/(I - K + gamma*K)
    if(round(sum(pi),5)!=1){print("Error! Something with pi is wrong")}
    
    for(k in 1:K){
      #print(k)
      zhat_k <- Zhat[,k]
      which_vals <- which(zhat_k>.001)
      if(length(which_vals)==0){which_vals <- 1:I}
      obj <- function(par){
        lik_term <- sum(unlist(lapply(which_vals,function(i){
          zhat_k[i]*dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=par[1:J],theta=par[J+1],M=M,log=TRUE)
        })))
        prior_term <- sum(dbeta(par[1:J],a,b,log=TRUE))+dgamma(par[J+1],gamma1,gamma2,log=TRUE)
        return(-lik_term-prior_term)
      }
      ptheta[,k] <- optim(ptheta[,k],obj,method="L-BFGS-B",lower=rep(0.00001,J+1),upper=c(rep(1-.00001,J),qgamma(.9999,gamma1,gamma2)),
                          control=list(maxit=10))$par
    }
    
    lik_term <- sum(unlist(lapply(1:I,function(i){
      logSumExp(unlist(lapply(1:K,function(k){
        log(pi[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=ptheta[1:J,k],theta=ptheta[J+1,k],M=M,log=TRUE)
      })))
    })))
    prior_term <- log(ddirichlet(pi,rep(gamma,K)))+sum(dbeta(ptheta[1:J,],a,b,log=TRUE))+
      sum(dgamma(ptheta[J+1,],gamma1,gamma2,log=TRUE))
    newobj <- lik_term+prior_term
    diff <- abs(newobj - currobj)
    
    print(diff)
    print(c(newobj,currobj))
    currobj <- newobj
    nsteps <- nsteps+1
    #print(ptheta)
  }
  
  list(phat=ptheta[1:J,],thetahat=ptheta[J+1,],pihat=pi,zhat=Zhat,obj=currobj)
  
  
}

#### MFMM Estimation ####
btlb_mfm <- function(Pi,X,M,gamma,lambda,a,b,gamma1,gamma2,Pi_full=NULL,
                     startK = 1, mh_pjk = 0.01, mh_thetak = 1,
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
  pi_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  p_all <- array(NA,dim=c(J,I,length(which_keep)))
  theta_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  Z_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  accept_p <- c()
  accept_theta <- c()
  
  ## Initialize at Random
  K <- startK
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
    
    
    ## Step 3: Update K (and gamma; skipped here)
    probs <- unlist(lapply(Kplus:(Kplus+100),function(k){
      dpois(k-1,lambda,log=T)+lfactorial(k)-lfactorial(k-Kplus)+lgamma(gamma*k)-lgamma(I+gamma*k)
    }))
    K <- sample(Kplus:(Kplus+100),1,prob = exp(probs-logSumExp(probs)))
    
    
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
        pi_all[it,1:K] <- pi
        Z_all[it,] <- Z
        if(K>Kplus){
          p_all[,(Kplus+1):K,it] <- ptheta[1:J,(Kplus+1):K] #update only the first Kplus parameters,
          theta_all[it,(Kplus+1):K] <- ptheta[J+1,(Kplus+1):K] #(will update the rest later!)
        }
      }}
    
    iter <- iter + mh_iters
    if(iter >= progress_iters[1]){
      print(paste0((11-length(progress_iters))*10,"% Complete: Iteration ",iter-1," out of ",max(which_keep)))
      print(paste0("Current K = ",K,"; Kplus = ",Kplus))
      progress_iters <- progress_iters[-1]
    }
  }
  print(paste0("Done! Saving ",length(which_keep)," estimate iterations after burning/thinning"))
  return(list(Z=Z_all,p=p_all[,1:max(K_all),],theta=theta_all[,1:max(K_all)],
              pi=pi_all[,1:max(K_all)],K=K_all,Kplus=Kplus_all,
              accept_p=(cumsum(accept_p)/1:length(accept_p))[which_keep],
              accept_theta=(cumsum(accept_theta)/1:length(accept_p))[which_keep],
              max_iters=max_iters,mh_iters=mh_iters,burn=burn,thin=thin,seed=seed))
}

#### Sandbox ####

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

gamma <- 1
lambda <- 3

res <- btlb_mfm(Pi=Pi,X=X,M=M,gamma=gamma,lambda=lambda,a=a,b=b,
                gamma1=gamma1,gamma2=gamma2,
                Pi_full=NULL,mh_pjk = .05, mh_thetak = 5,
                startK=15,max_iters=200,burn=0,thin=10)
par(mfrow=c(2,2))
plot(res$accept_p,ylim=c(0,1),type="l",ylab="Accept Prob for p")
plot(res$accept_theta,ylim=c(0,1),type="l",ylab="Accept Prob for theta")
plot(res$K,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="K")
plot(res$Kplus,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="Kplus")
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

gamma <- 1
lambda <- 1

res <- btlb_mfm(Pi=Pi,X=X,M=M,gamma=gamma,lambda=lambda,a=a,b=b,
                gamma1=gamma1,gamma2=gamma2,
                Pi_full=Pi_full,mh_pjk = .05, mh_thetak = 5,
                startK=17,max_iters=50,mh_iters=10,burn=.10,thin=10)
par(mfrow=c(2,2))
plot(res$accept_p,ylim=c(0,1),type="l",ylab="Accept Prob for p")
plot(res$accept_theta,ylim=c(0,1),type="l",ylab="Accept Prob for theta")
plot(res$K,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="K")
plot(res$Kplus,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="Kplus")
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

gamma <- .1
lambda <- 5

res <- btlb_mfm(Pi=Pi,X=X,M=M,gamma=gamma,lambda=lambda,a=a,b=b,
                gamma1=gamma1,gamma2=gamma2,
                Pi_full=NULL,mh_pjk = .20, mh_thetak = 5,
                startK=50,max_iters=20,mh_iters=5,burn=0,thin=1)
names(res)
res$thin
par(mfrow=c(2,2))
plot(res$accept_p,ylim=c(0,1),type="l",ylab="Accept Prob for p")
plot(res$accept_theta,ylim=c(0,1),type="l",ylab="Accept Prob for theta")
plot(res$K,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="K")
plot(res$Kplus,type="l",ylim=c(1,max(res$K,res$Kplus)),ylab="Kplus")
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




