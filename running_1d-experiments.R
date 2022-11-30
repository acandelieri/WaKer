rm(list=ls()); graphics.off(); cat("\014")

source("core.R")
source("test_problems.R")

dgt = 2
m = 2 # 1D probability simplex

costM = 4 * (matrix(1,m,m) - diag(m))
p = 2

# set.seed(42)


cat("> Selecting the 1D test problem:\n")
cat("[1] Test problem #02\n")
cat("[2] Test problem #13\n")
cat("[3] Test problem #15\n")
cat("[4] Test problem 'modified Xiong' function\n\n")
reply=""
while( !(reply %in% c("1","2","3","4")) ) {
  reply = readline("> Your choice:")
}

f = switch (reply,
  "1" = test.problem02,
  "2" = test.problem13,
  "3" = test.problem15,
  "4" = test.modifiedXiong
)

runs = as.integer(readline("> Number of independent runs: "))
N = as.integer(readline("> Size of the dataset, N: "))
S = as.integer(readline("> Size of the testset, S: "))


As = round(seq(0,1,length.out=S),dgt); As = cbind(As,1-As)
# *******************************************************************
# Patch for the computation of WST with transport library
ixs = which(As[,1]==1)
if( length(ixs)>0 ) {
  As[ixs,1] = 1-10^-dgt; As[ixs,2] = 10^-dgt 
}
ixs = which(As[,2]==1)
if( length(ixs)>0 ) {
  As[ixs,2] = 1-10^-dgt; As[ixs,1] = 10^-dgt 
}
# *******************************************************************
ys = apply(As,1,f)

RMSE.kse = matrix(NA,runs,N-m)
RMSE.kwse = matrix(NA,runs,N-m)


for( run in 1:runs ) {
  
  set.seed(run)
  
  cat("\014> Run",run,"out of",runs,"...\n")
  
  # sample N distributions
  A = round(runif(2*N,0,1),dgt); A = cbind(A,1-A)
  A = unique(A)
  stopifnot(nrow(A)>=N)
  A = A[1:N,]
  
  # *******************************************************************
  # Patch for the computation of WST with transport library
  ixs = which(A[,1]==1)
  if( length(ixs)>0 ) {
    A[ixs,1] = 1-10^-dgt; A[ixs,2] = 10^-dgt 
  }
  ixs = which(A[,2]==1)
  if( length(ixs)>0 ) {
    A[ixs,2] = 1-10^-dgt; A[ixs,1] = 10^-dgt 
  }
  # *******************************************************************
  
  y = apply(A,1,f)

  
  
  n = m+1
  stop.kse = stop.kwse = F
  
  while( n<=N && !(stop.kse && stop.kwse) ) {
    
    cat("\n***** [",n,"/",N,"] *****\n\n",sep="")
    
    gp = NULL
    if( !stop.kse )
      try( expr=(gp = train.kse.GP( X=A[1:n,], y=y[1:n], isotropic=F, ym=0, nrep=10 )),
           silent=T ) 
    if( is.null(gp) ) {
      stop.kse = T
    } else {
      # predict with the gp
      gp.ys = numeric(length(ys))
      for( i in 1:length(ys) ) {
        pred = predict.kse.GP( x=As[i,], gp.hpars=gp, X=A[1:n,], y=y[1:n] )  
        gp.ys[i] = pred$mu
      }
      RMSE.kse[run,n-m] = sqrt(mean((ys-gp.ys)^2))
    }
    
    
    
    cat("\n")
    
    wgp = NULL
    if( !stop.kwse )
      try( expr=(wgp = train.kwse.GP( X=A[1:n,], y=y[1:n], isotropic=F, ym=0, nrep=10, p=p, costM=costM )),
           silent=T )
    if( is.null(wgp) ) {
      stop.kwse = T
    } else {
      # predict with the wasserstein-gp
      wgp.ys = numeric(length(ys))
      for( i in 1:length(ys) ) {
        pred = predict.kwse.GP( x=As[i,], gp.hpars=gp, X=A[1:n,], y=y[1:n], p=p, costM=costM )
        wgp.ys[i] = pred$mu
      }
      RMSE.kwse[run,n-m] = sqrt(mean((ys-wgp.ys)^2))
    }
    
    cat("\n")
    
    n = n+1
  }
}

RMSEs = list(RMSE.kse=RMSE.kse,RMSE.kwse=RMSE.kwse)

save.image( file=paste0("Results_1D_",reply,"_",runs,"_",N,"_",S,".RData") )
