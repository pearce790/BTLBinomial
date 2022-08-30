#### Packages ####

library(gtools)
library(matrixStats)

#### Functions ####
dbtlb <- function(Pi,X,p,theta,M,log=FALSE,Pi_full=NULL){
  if(!is.matrix(X)){stop("X must be a matrix of ratings")}
  if(length(p)!=ncol(X)){stop("length(p) must equal ncol(X)")}
  worth <- exp(-theta*p)
  logd <- dbtl(Pi,worth,log=T,Pi_full=Pi_full)+sum(apply(X,1,function(x){dbinom(x,M,p,log=T)}),na.rm=T)
  if(log){return(logd)
  }else{return(exp(logd))}
}
dbtl <- function(Pi,worth,log=FALSE,Pi_full=NULL){
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

