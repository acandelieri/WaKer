library(transport)

# Squred Exponential kernel
kse <- function( x, y, sf2=1, l ) {
  stopifnot( length(x)==length(y) && (length(l)==1 || length(l)==length(x)) && length(sf2)==1 )
  kse = sf2 * exp( -0.5 * sum( ((x-y)/l)^2 ) )
}

# Exponential kernel
kexp <- function( x, y, sf2=1, l ) {
  stopifnot( length(x)==length(y) && (length(l)==1 || length(l)==length(x)) && length(sf2)==1 )
  kexp = sf2 * exp( - sum( abs((x-y)/l) ) )
}

# Wasserstein Squared Exponential kernel
kwse <- function( x, y, sf2=1, r, p, costM ) {
  stopifnot( length(x)==length(y) && (length(r)==1 || length(r)==length(x)) && length(sf2)==1 )
  
  if( length(r)==1 ) {
    # isotropic
    Wp = wasserstein( a=x, b=y, p=p, costm=costM )
    kwse = sf2 * exp( -0.5 * (Wp/r)^2 )
  } else {
    # anisotropic0
    res = transport( a=x, b=y, p=p, costm=costM )
    TP = matrix(0,length(x),length(y))
    for( ii in 1:nrow(res) ) 
      TP[res$from[ii],res$to[ii]] = res$mass[ii]
    M = TP * costM * (1/r)^2 %*% t(rep(1,length(x)))

    Wp = sum(M)^(1/p)
    kwse = sf2 * exp( -0.5 * Wp^2 )
  }
}

# Function computing the kernel matrix
K.matrix <- function( X, kernel, sf2=1, ... ) {
  
  stopifnot( identical(kernel,kse) || identical(kernel,kwse) )
  args = list(...)
  
  if( identical(kernel,kwse) ) {
    r = args$r
    p = args$p
    costM = args$costM
  } else {
    l = args$l
  }

  K = matrix(NA,nrow(X),nrow(X))
  for( i in 1:nrow(K) ) {
    for( j in i:ncol(K) ) {
      if( identical(kernel,kwse) ) {
        Kij = kwse( x=X[i,], y=X[j,], sf2=sf2, r=r, p=p, costM=costM )
      } else {
        Kij = kse( x=X[i,], y=X[j,], sf2=sf2, l=l )
      }
      K[i,j] = Kij
      if( i != j )
        K[j,i] = Kij
    }
  }
  K.matrix = K
}

# Marginal Likelihood, given a Kernel matrix
MLE <- function( K, y, ym ) {
  
  # inversion by Cholesky decomposition
  # L = chol(K); L=t(L) # to have lower triangular as in Williams and Rasmussen
  # M = backsolve(L,y-ym,upper.tri=F)
  # alpha = backsolve(t(L),M,upper.tri=T)
  # term1 = -0.5 * t(y-ym) %*% alpha
  # term2 = -0.5 * log(aum(diag(L)))
  # term3 = -0.5 * length(y) * log(2*pi)
  
  # inversion by solve
  K.inv = NULL
  try( expr=(K.inv = solve(K)), silent=T )
  if( is.null(K.inv) ) {
    mle = NA
  } else {
    term1 = -0.5 * t(y-ym) %*% K.inv %*% (y-ym)
    term2 = -0.5 * log(det(K))
    term3 = -0.5 * length(y) * log(2*pi)
  
    mle = term1 + term2 + term3
  }
  
  MLE = mle
}

# Marginal Likelihood as an objective function to maximize
MLE.obj <- function( theta, X, y, ym, k.fun, ... ) {
  
  stopifnot( identical(k.fun,kse) || identical(k.fun,kwse) )
  
  sf2 = theta[1]
  l = theta[-1]
  if( identical(k.fun,kse) ) {
    K = K.matrix( X=X, kernel=kse, sf2=sf2, l=l )
  } else {
    args = list(...)
    K = K.matrix( X=X, kernel=kwse, sf2=sf2, r=l, p=args$p, costM=args$costM )
  }
  
  mle = MLE( K=K, y=y, ym=ym )
  
  if( is.na(mle) ) {
    # cat("***** Nugget estimation required!\n")
    
    # res = optimize( f=MLE.with.nugget, interval=c(.Machine$double.eps,10^-5), maximum=T,
    #                 K=K, y=y, ym=ym )
    res = optimize( f=MLE.with.nugget, interval=c(.Machine$double.eps,var(y)), maximum=T,
                    K=K, y=y, ym=ym )
    
    NUGGET.GLOBAL <<- res$maximum
    mle = res$objective
    # cat(" --> ", mle,"\n")
  }
  
  MLE.obj = mle
}

# Marginal likelihood optimized only depending on nugget
MLE.with.nugget <- function( nugget, K, y, ym ) {

  K = K + nugget * diag(nrow(K))

  mle = MLE( K=K, y=y, ym=ym )
  if( is.na(mle) )
    stop("* ATTENZIONE! MLE.with.nugget è NA!\n")

  MLE.with.nugget = mle
}





# training a GP with SE kernel
train.kse.GP <- function( X, y, isotropic=F, ym=0, nrep=10 ) {
  
  cat("> Training a GP with SE kernel...\n")
  
  l.lo = 10^-1; l.up = 1
  
  sols = vals = nuggets = NULL
  
  for( iter in 1:nrep ) {
  
    NUGGET.GLOBAL <<- 0 # it is changed in case that nugget estimation is needed
    cat(".")
    
    sf2_0 = runif(1,10^-5,var(y))
    
    if(isotropic) {
      # theta0 = c(sf2_0,runif(1,10^-ncol(X),2*diff(apply(X,2,range))))
      theta0 = c(sf2_0,runif(1,l.lo,l.up))
    } else {
      # theta0 = c(sf2_0,runif(ncol(X),10^-ncol(X),2*diff(apply(X,2,range))))
      theta0 = c(sf2_0,runif(ncol(X),l.lo,l.up))
    }
    
    # cat("- initial theta:",theta0,"\n")
    res = optim( par=theta0, fn=MLE.obj, gr=NULL, method="L-BFGS-B",
                 # lower=c(10^-5,rep(10^-ncol(X),ncol(X))), upper=c(var(y),2*diff(apply(X,2,range))),
                 lower=c(10^-5,rep(l.lo,ncol(X))), upper=c(var(y),rep(l.up,ncol(X))),
                 control=list(fnscale=-1,trace=0),
                 X=X, y=y, k.fun=kse, ym=ym )
    # cat("* optimal theta:",res$par,"\n")
    # cat("* with value:",res$value,"\n")
    # cat("* and with nugget:",NUGGET.GLOBAL,"\n")
    sols = rbind( sols, res$par )
    vals = c( vals, res$value )
    
    nugget = NUGGET.GLOBAL
    nuggets = c( nuggets, nugget )
  
  }
  cat("Done!\n")
  
  sols = round(sols,8)
  vals = round(vals,8)
  
  ixs = which(vals==max(vals))
  if( length(ixs)>1 )
    ixs = ixs[ which.min(nuggets[ixs]) ]

  sol = sols[ixs,]
  
  cat("\n[Evaluated Configurations]\n")
  print( data.frame( theta=sols, MLE=vals, nugget=nuggets ) )
  
  # ## Just to check!
  # cat("\n$$$ JUST FOR CHECK:\n")
  # cat("ix =",ixs,"\n")
  # cat("MLE =", MLE.obj( theta=sol, X=X, y=y, ym=ym, k.fun=kse ),"\n" )
  # cat("NUGGET.GLOBAL =", NUGGET.GLOBAL,"\n" )
  # print( K.matrix(X=X,kernel=kse,sf2=sol[1],l=sol[-1] ) + nuggets[ixs] * diag(length(y)) )
  
  train.kse.GP = list( sf2=sol[1], l=sol[-1], nugget=nuggets[ixs] )
  
}

# training a GP with WSE kernel
train.kwse.GP <- function( X, y, isotropic=F, p, costM, ym=0, nrep=10 ) {

  cat("> Training a GP with (isotropic) Wasserstein SE kernel...\n")

  l.lo = 10^-1; l.up = 1

  sols = vals = nuggets = NULL
  
  
  for( iter in 1:nrep ) {
    
    NUGGET.GLOBAL <<- 0 # it is changed in case that nugget estimation is needed
    cat(".")
  
    sf2_0 = runif(1,10^-5,var(y))
    if(isotropic) {
      # theta0 = c(sf2_0,runif(1,10^-ncol(X),2*diff(apply(X,2,range))))
      theta0 = c(sf2_0,runif(1,l.lo,l.up))
    } else {
      # theta0 = c(sf2_0,runif(ncol(X),10^-ncol(X),2*diff(apply(X,2,range))))
      theta0 = c(sf2_0,runif(ncol(X),l.lo,l.up))
    }
    # cat("- initial theta:",theta0,"\n")
    res = optim( par=theta0, fn=MLE.obj, gr=NULL, method="L-BFGS-B",
                 # lower=c(10^-5,rep(10^-ncol(X),ncol(X))), upper=c(var(y),2*diff(apply(X,2,range))),
                 lower=c(10^-5,rep(l.lo,ncol(X))), upper=c(var(y),rep(l.up,ncol(X))),
                 control=list(fnscale=-1,trace=0),
                 X=X, y=y, k.fun=kwse, ym=ym, p=p, costM=costM )
    # cat("* optimal theta:",res$par,"\n")
    sols = rbind( sols, res$par )
    vals = c( vals, res$value )
    nugget = NUGGET.GLOBAL

    nuggets = c( nuggets, nugget )

  }
  cat("Done!\n")

  sols = round(sols,8)
  vals = round(vals,8)

  ixs = which(vals==max(vals))
  if( length(ixs)>1 )
    ixs = ixs[ which.min(nuggets[ixs]) ]

  sol = sols[ixs,]

  cat("\n[Evaluated Configurations]\n")
  print( data.frame( theta=sols, MLE=vals, nugget=nuggets ) )

  train.kwse.GP = list( sf2=sol[1], l=sol[-1], nugget=nuggets[ixs] )

}

# predicting through a GP with a SE kernel
predict.kse.GP <- function( x, gp.hpars, X, y ) {
  
  stopifnot( is.vector(x) && is.matrix(X) )
  
  K = K.matrix( X=X, kernel=kse, sf2=gp.hpars$sf2, l=gp.hpars$l )
  K = K + gp.hpars$nugget * diag(nrow(K))
  
  K.inv = solve(K)
  
  kv = NULL
  for( i in 1:nrow(X) )
    kv = c(kv,kse(x=x,y=X[i,],sf2=gp.hpars$sf2,l=gp.hpars$l))
  
  mu = t(kv) %*% K.inv %*% y
  sg = gp.hpars$sf2 - t(kv) %*% K.inv %*% kv
  
  
  return( list( mu=mu, sg=sg ) )
}

# predicting through a GP with a WSE kernel
predict.kwse.GP <- function( x, gp.hpars, X, y, p, costM ) {

  stopifnot( is.vector(x) && is.matrix(X) )

  K = K.matrix( X=X, kernel=kwse, sf2=gp.hpars$sf2, r=gp.hpars$l, p=p, costM=costM )
  K = K + gp.hpars$nugget * diag(nrow(K))

  K.inv = solve(K)

  kv = NULL
  for( i in 1:nrow(X) )
    kv = c(kv,kwse(x=x,y=X[i,],sf2=gp.hpars$sf2,r=gp.hpars$l,p=p,costM=costM))

  mu = t(kv) %*% K.inv %*% y
  sg = gp.hpars$sf2 - t(kv) %*% K.inv %*% kv


  return( list( mu=mu, sg=sg ) )
}