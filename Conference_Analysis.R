#### Sushi Reproducibility ####
#source("~/BTLBinomial/source.R")
source("source.R")

#### Choose Hyperparameters and Run Simulation ####
I <- 50
J <- 50
results <- matrix(NA,nrow=0,ncol=4+2*J+1+J)
iter <- as.numeric(commandArgs(trailingOnly=T)[1])

for(M in c(4,9)){
  for(R in c(4,8,12,24)){
    for(theta in c(1,5,10,20,40)){
      #initialize
      print("initializing")
      set.seed(iter*M*R*theta)
      print(c(M,R,theta,iter))
      X <- matrix(NA,nrow=I,ncol=J)
      Pi <- matrix(NA,nrow=I,ncol=4)
      Pi_full <- matrix(NA,nrow=I,ncol=R)
      
      #make assignments at random
      print("assignments")
      continue <- FALSE
      while(continue==FALSE){
        assignments <- matrix(0,nrow=I,ncol=J)
        tryCatch({
          for(i in 1:I){
            assignments[i,sample(which(apply(assignments,2,sum)<R),R)] <- 1
          }
        },error=function(e){})
        if(all(apply(assignments,2,sum)==R)){continue <- TRUE}
      }
      
      #draw data
      print("draw data")
      p <- runif(J,min=0,max=1)
      for(i in 1:I){
        which_proposals <- which(assignments[i,]==1)
        X[i,which_proposals] <- rbinom(n=R,size=M,prob=p[which_proposals])
        Pi[i,1:4] <- sample(which_proposals,4,prob=exp(-theta*p[which_proposals]))
        Pi_full[i,1:4] <- Pi[i,1:4]
        if(R>4){Pi_full[i,5:R] <- setdiff(which_proposals,Pi[i,1:4])}
      }
      
      print("estimate")
      tryCatch({
        map_est <- btlb_map(Pi=Pi,X=X,M=M,Pi_full=Pi_full,K=1,gamma=1,a=1,b=1,gamma1=5,gamma2=0.25,tol=0.001,maxit=100,verbose=FALSE)
        results <- rbind(results,c(M,R,theta,iter,p,map_est$p,map_est$theta,apply(X,2,function(x){mean(x,na.rm=TRUE)})))
      },error=function(e){
        print("Fail MAP Search for M,R,theta,iter:")
        print(c(M,R,theta,iter))
      })
    }
  }
}

results <- as.data.frame(results)
names(results) <- c("M","R","theta","Iter",paste0("p",1:J),paste0("pmap",1:J),
                    "thetamap",paste0("Xbar",1:J))

save(results,file=paste0("results/conf_simulation",iter,".RData"))
