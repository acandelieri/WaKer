rm(list=ls()); graphics.off(); cat("\014")

library( plot3D )

source("test_problems.R")
source("core.R")


set.seed(1)


f = test.problem02
# f = test.problem13
# f = test.problem15


p = 2 # MUST BE p=2 !!!!!
costM.bin = 1 - diag(2)
costM.euc = costM.bin * 4
costM = costM.euc

n0 = 7

random.sampling = T


if( !random.sampling ) {
  A = seq(0,1,length.out=n0); A = cbind(A,1-A)
} else {
  A = round(runif(10*n0),4); A = cbind(A,1-A)
}
A = unique(A)
A = A[sample(1:nrow(A),n0,F),]

y = apply( A, 1, f )

ym = 0
# ym = mean(y)


# ---------------------------------------------------------
# PATCH!
# ---------------------------------------------------------
A[which(A==1,arr.ind=T)] = 1-10^-8
A[which(A==0,arr.ind=T)] = 10^-8
# ---------------------------------------------------------



l.min=10^-1; l.max=max(2*diff(apply(A,2,range))); ngrid=50


# isotropic kernel (1D MLE charts: MLE wrt length scale)
cat("> isotropic kernel (1D MLE charts: MLE wrt length scale)\n")
sf2 = 1
ls = seq(l.min,l.max,length.out=ngrid)

MLE.kse.values = MLE.kwse.values = NULL
ct.kse = ct.kwse = NULL
for( l in ls ) {
  
  K.se = K.matrix( X=A, kernel=kse, sf2=sf2, l=rep(l,ncol(A)) )
  mle = NA
  elapse = Sys.time()
  try( expr=(mle=MLE(K=K.se,y=y,ym=ym)), silent=T )
  ct.kse = c(ct.kse,difftime(Sys.time(),elapse,units="secs"))
  MLE.kse.values = c(MLE.kse.values,mle)
  
  K.wse = K.matrix( X=A, kernel=kwse, sf2=sf2, r=rep(l,ncol(A)), p=p, costM=costM )
  mle = NA
  elapse = Sys.time()
  try( expr=(mle=MLE(K=K.wse,y=y,ym=ym)), silent=T )
  ct.kwse = c(ct.kwse,as.numeric(difftime(Sys.time(),elapse,units="secs")))
  MLE.kwse.values = c(MLE.kwse.values,mle)
}
cat("> Computational time MLE for SE Kernel:",sum(ct.kse),"\n")
cat("> Computational time MLE for WSE Kernel:",sum(ct.kwse),"\n")

par(mfrow=c(2,1))
MIN.Y = min( MLE.kse.values, MLE.kwse.values, na.rm=T )
MAX.Y = max( MLE.kse.values, MLE.kwse.values, na.rm=T )
# plot( ls, MLE.kse.values, type="l", col="blue", lwd=3, ylim=c(MIN.Y,MAX.Y) )
# plot( ls, MLE.kwse.values, type="l", col="red", lwd=3, ylim=c(MIN.Y,MAX.Y) )
plot( ls, MLE.kse.values, type="l", col="blue", lwd=3 )
plot( ls, MLE.kwse.values, type="l", col="red", lwd=3 )
par(mfrow=c(1,1))




# isotropic kernel (2D MLE charts: MLE wrt length scale and sf2)
cat("> isotropic kernel (2D MLE charts: MLE wrt length scale and sf2)\n")
ls = seq(l.min,l.max,length.out=ngrid)
sf2s = seq(10^-1,var(y),length.out=ngrid)
MLE.kse.values = MLE.kwse.values = NULL
ct.kse = ct.kwse = NULL
for( l in ls ) {
  for( sf2 in sf2s ) {
    
    K.se = K.matrix( X=A, kernel=kse, sf2=sf2, l=rep(l,ncol(A)) )
    mle = NA
    elapse = Sys.time()
    try( expr=(mle=MLE(K=K.se,y=y,ym=ym)), silent=T )
    ct.kse = c(ct.kse,difftime(Sys.time(),elapse,units="secs"))
    MLE.kse.values = c(MLE.kse.values,mle)
    
    K.wse = K.matrix( X=A, kernel=kwse, sf2=sf2, r=rep(l,ncol(A)), p=p, costM=costM )
    mle = NA
    elapse = Sys.time()
    try( expr=(mle=MLE(K=K.wse,y=y,ym=ym)), silent=T )
    ct.kwse = c(ct.kwse,as.numeric(difftime(Sys.time(),elapse,units="secs")))
    MLE.kwse.values = c(MLE.kwse.values,mle)
  }
}
cat("> Computational time MLE for SE Kernel:",sum(ct.kse),"\n")
cat("> Computational time MLE for WSE Kernel:",sum(ct.kwse),"\n")


par(mfrow=c(2,2))
persp3D( x=ls, y=sf2s, z=matrix(MLE.kse.values,length(ls)), col=rev(heat.colors(100)) )
persp3D( x=ls, y=sf2s, z=matrix(MLE.kwse.values,length(ls)), col=rev(heat.colors(100)) )

image2D( x=ls, y=sf2s, z=matrix(MLE.kse.values,length(ls)), NAcol="grey", col=rev(heat.colors(100)) )
contour2D( x=ls, y=sf2s, z=matrix(MLE.kse.values,length(ls)), col="black", nlevels=10, add=T )
image2D( x=ls, y=sf2s, z=matrix(MLE.kwse.values,length(ls)), NAcol="grey", col=rev(heat.colors(100)) )
contour2D( x=ls, y=sf2s, z=matrix(MLE.kwse.values,length(ls)), col="black", nlevels=10, add=T )
par(mfrow=c(1,1))





# anisotropic kernel (2D MLE charts: MLE wrt 2D length scale - sf2=1)
cat("> anisotropic kernel (2D MLE charts: MLE wrt 2D length scale - sf2=1)\n")
ls = seq(l.min,l.max,length.out=ngrid)
sf2 = 1
MLE.kse.values = MLE.kwse.values = NULL
ct.kse = ct.kwse = NULL
nugget = 0
if( nugget != 0 ) {
  invisible(readline(("\n> STAI UANDO UN NUGGET !=0... SE VA BENE, PREMI [INVIO]")))
}
for( l1 in ls ) {
  for( l2 in ls ) {
    
    K.se = K.matrix( X=A, kernel=kse, sf2=sf2, l=c(l1,l2) )
    K.se = K.se + nugget * diag(nrow(K.se))
    mle = NA
    elapse = Sys.time()
    try( expr=(mle=MLE(K=K.se,y=y,ym=ym)), silent=T )
    ct.kse = c(ct.kse,difftime(Sys.time(),elapse,units="secs"))
    MLE.kse.values = c(MLE.kse.values,mle)
    
    K.wse = K.matrix( X=A, kernel=kwse, sf2=sf2, r=c(l1,l2), p=p, costM=costM )
    K.wse = K.wse + nugget * diag(nrow(K.wse))
    mle = NA
    elapse = Sys.time()
    try( expr=(mle=MLE(K=K.wse,y=y,ym=ym)), silent=T )
    ct.kwse = c(ct.kwse,as.numeric(difftime(Sys.time(),elapse,units="secs")))
    MLE.kwse.values = c(MLE.kwse.values,mle)
  }
}
cat("> Computational time MLE for SE Kernel:",sum(ct.kse),"\n")
cat("> Computational time MLE for WSE Kernel:",sum(ct.kwse),"\n")

mar0 = par("mar")
par(mfrow=c(2,2))
par(mar=c(1.1,3.1,1.1,1.1) )
persp3D( x=ls, y=ls, z=matrix(MLE.kse.values,length(ls)), col=rev(heat.colors(100)),
         xlab="", ylab="", zlab="MLE", colkey=F, cex.lab=1.5, cex.axis=1.5 )
persp3D( x=ls, y=ls, z=matrix(MLE.kwse.values,length(ls)), col=rev(heat.colors(100)),
         xlab="", ylab="", zlab="MLE", colkey=F, cex.lab=1.5, cex.axis=1.5 )

par(mar=c(4.6,4.6,2.1,1.1) )
image2D( x=ls, y=ls, z=matrix(MLE.kse.values,length(ls)), NAcol="grey", col=rev(heat.colors(100)),
         xlab=expression(l[1]), ylab=expression(l[2]), colkey=F,
         cex.lab=2, cex.axis=2, cex.main=2, main="SE kernel" )
contour2D( x=ls, y=ls, z=matrix(MLE.kse.values,length(ls)), col="black", nlevels=10, add=T )
image2D( x=ls, y=ls, z=matrix(MLE.kwse.values,length(ls)), NAcol="grey", col=rev(heat.colors(100)),
         xlab=expression(rho[1]), ylab=expression(rho[2]), colkey=F, 
         cex.lab=2, cex.axis=2, cex.main=2, main="WSE kernel" )
contour2D( x=ls, y=ls, z=matrix(MLE.kwse.values,length(ls)), col="black", nlevels=10, add=T )
par(mar=mar0)
par(mfrow=c(1,1))
