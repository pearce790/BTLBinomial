#### Packages ####
library(gtools) # for dirichlet distribution functions
library(matrixStats) #for logSumExp function (helps with rounding errors)
library(fipp) #for calculating prior on Kplus
library(salso) #for cluster assignment plotting
library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)

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
map_btlb <- function(Pi,X,M,Pi_full=NULL,K,gamma,a,b,gamma1,gamma2,
                     tol=1,maxit=50,verbose=TRUE,seed=NULL){
  
  
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
                         theta=ptheta[J+1,k],M=M,log=T,Pi_full=Pi_full[i,])
      }))
      Z[i,] <- exp(lprobs-logSumExp(lprobs))
    }

    
    # M-Step:
    ## Update (p,theta)
    t2 <- Sys.time()
    for(k in 1:K){
      #print(k)
      zhat_k <- Z[,k]
      # which_vals <- which(zhat_k>.001)
      # if(length(which_vals)==0){which_vals <- 1:I}
      obj <- function(par){
        worthk <- exp(-par[1:J]*par[J+1])
        lik_term <- sum(unlist(lapply(1:I,function(i){
          zhat_k[i]*(dbtl(Pi=Pi[i,],worth=worthk,log=T,Pi_full=Pi_full[i,])+sum(X[i,]*log(par[1:J])+(M-X[i,])*log(1-par[1:J])))
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

# #### Prior Distribution on K+ Example ####
# 
# pmfstatic2 <- nClusters(Kplus=1:20,N=10,type="static",gamma=5,maxK=50)
# dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = 7))
# plot(dens)
# sum((1:length(dens))*dens)
# #### MFMM Sandbox ####
# 
# set.seed(1)
# J <- 10
# M <- 4
# dat1 <- rbtlb(I=20,p=runif(J),theta=10,M=M)
# dat2 <- rbtlb(I=10,p=runif(J),theta=20,M=M)
# X <- rbind(dat1$X,dat2$X)
# Pi <- rbind(dat1$Pi,dat2$Pi)
# rm(dat1,dat2)
# I <- nrow(X)
# Pi_full <- NULL
# a <- 1
# b <- 1
# gamma1 <- 10
# gamma2 <- 0.5
# lambda <- 3
# 
# res <- btlb_mfm(Pi=Pi,X=X,M=M,lambda=lambda,a=a,b=b,
#                 gamma1=gamma1,gamma2=gamma2,gamma_hyp1=3,gamma_hyp2=2,
#                 Pi_full=NULL,mh_pjk = .05, mh_thetak = 5,mh_gamma=0.5,
#                 startK=15,max_iters=100,mh_iters=10,burn=0,thin=2)
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
# ggplot(plot_p,aes(x=Var3,y=value,group=jk,color=factor(jk)))+
#   geom_line()+theme(legend.position="none")+ylim(c(0,1))+
#   ylab("p")+xlab("Iteration (after burn/thin)")+
#   ggtitle("Trace Plot: p")
# 
# 
# 
# 
# #### MAP Sandbox ####
# set.seed(1)
# J <- 10
# M <- 4
# dat1 <- rbtlb(I=20,p=runif(J),theta=10,M=M)
# dat2 <- rbtlb(I=10,p=runif(J),theta=20,M=M)
# X <- rbind(dat1$X,dat2$X)
# Pi <- rbind(dat1$Pi,dat2$Pi)
# rm(dat1,dat2,J)
# Pi_full <- NULL
# K <- 2
# gamma <- 1
# a <- 1
# b <- 1
# gamma1 <- 10
# gamma2 <- 0.5
# map_btlb(Pi=Pi,X=X,M=M,Pi_full=Pi_full,K=K,gamma=gamma,a=a,b=b,gamma1=gamma1,gamma2=gamma2,seed=1)
# 
# 
