rm(list=ls()); graphics.off(); cat("\014")

# filename template is "Results_<x>D_<test problem>_<independent runs>_<dataset size>_<test set size for RMSE computation>"

cat("> Select one of the experiments reported in the paper:\n")
cat(" [1] test problem 02\n")
cat(" [2] test problem 13\n")
cat(" [3] test problem 15\n")
cat(" [4] modified Xiong function\n")
cat(" [5] Bird function\n")

reply = ""
while( !(reply %in% c("1","2","3","4","5")) ) {
  reply = readline("> Choice:")
}

switch (reply,
  "1" = {load("~/Progetti R/WaKer/Results_1D_1_500_20_50.RData"); testProblemName = "test problem 02"},
  "2" = {load("~/Progetti R/WaKer/Results_1D_2_500_20_50.RData"); testProblemName = "test problem 13"},
  "3" = {load("~/Progetti R/WaKer/Results_1D_3_500_20_50.RData"); testProblemName = "test problem 15"},
  "4" = {load("~/Progetti R/WaKer/Results_1D_4_500_20_50.RData"); testProblemName = "modified Xiong function"},
  "5" = {load("~/Progetti R/WaKer/Results_2D_Bird_500_20_50.RData"); testProblemName = "Bird function"},
)


stopifnot( dim(RMSE.kse)==dim(RMSE.kwse) )

na.kse.perc = na.kwse.perc = NULL

for( i in 1:ncol(RMSE.kse) ) {
  na.kse.perc = c(na.kse.perc,round(100*length(which(is.na(RMSE.kse[,i])))/nrow(RMSE.kse),2))
  na.kwse.perc = c(na.kwse.perc,round(100*length(which(is.na(RMSE.kwse[,i])))/nrow(RMSE.kwse),2))
}

curr.mar = par("mar")
par(mar=c(4.1,4.6,2.1,1.1))
plot( (m+1):N, na.kwse.perc, type="o", lwd=3, col="red", ylim=c(0,max(na.kwse.perc,na.kse.perc)),
      xlab="dataset size", ylab="runs with failed GP learning [%]",
      cex.axis=2, cex.lab=2, cex.main=2, main=testProblemName)
lines( (m+1):N, na.kse.perc, type="o", lwd=3, col="blue" )
legend( "right", legend=c("WSE-GP via MLE","WSE-GP via Th. 4"), col=c("red","blue"), lwd=3, cex=2 )





kse.avg = kse.sd = numeric(ncol(RMSE.kse))+NA
for( i in 1:ncol(RMSE.kse) ) {
  kse.avg[i] = mean(RMSE.kse[,i],na.rm=T)
  kse.sd[i] = sd(RMSE.kse[,i],na.rm=T)
}

kwse.avg = kwse.sd = numeric(ncol(RMSE.kwse))+NA
for( i in 1:ncol(RMSE.kwse) ) {
  kwse.avg[i] = mean(RMSE.kwse[,i],na.rm=T)
  kwse.sd[i] = sd(RMSE.kwse[,i],na.rm=T)
}

kse.up = kse.avg+kse.sd; kse.lo = kse.avg-kse.sd
kwse.up = kwse.avg+kwse.sd; kwse.lo = kwse.avg-kwse.sd

kse.up[which( kse.up<0 )] = 0
kse.lo[which( kse.lo<0 )] = 0
kwse.up[which( kwse.up<0 )] = 0
kwse.lo[which( kwse.lo<0 )] = 0


reply = invisible(readline("> Press [ENTER] to visualize RMSE results."))

plot( kse.avg, type="p", col="blue", pch=15, main=testProblemName,
      cex.main=2, cex.axis=2, cex.lab=2, xlab = "dataset size", ylab="RMSE",
      ylim=c(min(kse.lo,kwse.lo,na.rm=T),max(kse.up,kwse.up,na.rm=T))  )
points( kwse.avg,pch=15, col="red" )
for( i in 1:length(kse.avg) )
  lines( rep(i,2), c(kse.up[i],kse.lo[i]), col="blue" )
for( i in 1:length(kwse.avg) )
  lines( rep(i,2), c(kwse.up[i],kwse.lo[i]), col="red" )

if( testProblemName=="Bird function" ) {
  leg.pos = "top"
} else {
  leg.pos = "topright"
}
  
legend( leg.pos, legend=c("WSE-GP via MLE", "WSE-GP via Th. 4"),
        col=c("red","blue"), lwd=4, cex=1.7 )


par(mar=curr.mar)