rm(list=ls()); graphics.off(); cat("\014")

# filename template is "Results_<x>D_<test problem>_<independent runs>_<dataset size>_<test set size for RMSE computation>"

notcomputable.perc <- function(RMSE) {
  na.perc = numeric(ncol(RMSE))
  for( i in 1:length(na.perc) ) {
    na.perc[i] = round(100*length(which(is.na(RMSE[,i])))/nrow(RMSE),2)
  }
  notcomputable.perc = na.perc
}

# ***********************************************************************************
# Large size Experiments
# ***********************************************************************************
load("~/Progetti R/WaKer/Results_1D_1_500_20_50.RData")

na.kse = na.kwse = list()
testNames = "test problem 02"
na.kwse[[1]] = notcomputable.perc(RMSE.kwse)

load("~/Progetti R/WaKer/Results_1D_2_500_20_50.RData")
testNames = c(testNames, "test problem 13")
na.kwse[[2]] = notcomputable.perc(RMSE.kwse)

load("~/Progetti R/WaKer/Results_1D_3_500_20_50.RData")
testNames = c(testNames, "test problem 15")
na.kwse[[3]] = notcomputable.perc(RMSE.kwse)

load("~/Progetti R/WaKer/Results_1D_4_500_20_50.RData")
testNames = c(testNames, "modified Xiong function")
na.kwse[[4]] = notcomputable.perc(RMSE.kwse)

# ***********************************************************************************


curr.mar = par("mar")
par(mar=c(4.1,4.6,2.1,1.1))
clrs = rainbow(4)
plot( (m+1):N, na.kwse[[1]], type="o", lwd=3, col="red", ylim=c(0,100),
      xlab="dataset size", ylab="runs with failed GP learning [%]",
      cex.axis=2, cex.lab=2, cex.main=2 )
for( i in 1:length(na.kwse) ) {
  lines( (m+1):N, na.kwse[[i]], type="o", lwd=3, col=clrs[i] )
  points( (m+1):N, na.kwse[[i]], col=clrs[i], pch=19 )
}
legend( "bottomright", legend=testNames, col=clrs, lwd=4, cex=2 )
par(mar=curr.mar)