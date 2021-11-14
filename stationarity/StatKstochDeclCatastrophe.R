#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

################################################################################
##  K DECLINING at -0.001 STOCHASTIC (5% SD initial size, Gaussian-resampled) ##
################################################################################

rm(list = ls())
library(ggplot2)
library(plotly)
library(expss)
library(ggridges)
library(ggpubr)

return.time <- function(ats) {
  lngTS <- length(ats)
  mnTS <- mean(ats, na.rm = TRUE)
  inout <- RTi <- numeric()
  for(i in 2:lngTS) {
    inout[i - 1] <- ifelse(mnTS > min(ats[i - 1], ats[i]) & mnTS < max(ats[i - 1], ats[i]), 0, 1)
  }
  ints <- inout.runs <- rep(rle(inout)$lengths, rle(inout)$lengths) * inout
  lng.runs <- length(ints)
  for(i in 2:lng.runs) {
    ints[i] <- ifelse(ints[i - 1] > 1, ints[i - 1] - 1, ints[i])
  }
  for(i in 2:lngTS) {
    RTi[i - 1] <- ints[i - 1] + ((ats[(i - 1) + ints[i - 1]] - mnTS) / (ats[(i - 1) + ints[i - 1]] - ats[i + ints[i - 1]]))
  }
  RTi <- RTi[!is.na(RTi)]
  lngRTi <- length(RTi[!is.na(RTi)])
  mRT <- mean(RTi, na.rm = TRUE)
  vRT <- var(RTi, na.rm = TRUE)
  return(list(mRT = mRT, vRT = vRT))
}

## load matrices
load("DPdeclKsm.RData")
load("PAdeclKsm.RData")
load("ZTdeclKsm.RData")
load("PHdeclKsm.RData")
load("VUdeclKsm.RData")
load("PGdeclKsm.RData")
load("SSdeclKsm.RData")
load("PTdeclKsm.RData")
load("SOdeclKsm.RData")
load("MNdeclKsm.RData")
load("ORdeclKsm.RData")
load("NRdeclKsm.RData")
load("GNdeclKsm.RData")
load("DNdeclKsm.RData")
load("ALdeclKsm.RData")
load("TCdeclKsm.RData")
load("THdeclKsm.RData")
load("SHdeclKsm.RData")
load("DMdeclKsm.RData")
load("TAdeclKsm.RData")
load("MRdeclKsm.RData")


# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.KstochDeclSm.iter <- dim(DP.Nmat2)[1]
DP.KstochDeclSm.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.KstochDeclSm.MRT <- DP.KstochDeclSm.VRT <- DP.KstochDeclSm.stat <- rep(NA,DP.KstochDeclSm.iter)
for(i in 1:DP.KstochDeclSm.iter) {
  DP.KstochDeclSm.MRT[i] <- DP.KstochDeclSm.RT[[i]]$mRT
  DP.KstochDeclSm.VRT[i] <- DP.KstochDeclSm.RT[[i]]$vRT
  DP.KstochDeclSm.stat[i] <- DP.KstochDeclSm.MRT[i] / DP.KstochDeclSm.VRT[i]
}
hist(DP.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.KstochDeclSm.iter <- dim(PA.Nmat2)[1]
PA.KstochDeclSm.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.KstochDeclSm.MRT <- PA.KstochDeclSm.VRT <- PA.KstochDeclSm.stat <- rep(NA,PA.KstochDeclSm.iter)
for(i in 1:PA.KstochDeclSm.iter) {
  PA.KstochDeclSm.MRT[i] <- PA.KstochDeclSm.RT[[i]]$mRT
  PA.KstochDeclSm.VRT[i] <- PA.KstochDeclSm.RT[[i]]$vRT
  PA.KstochDeclSm.stat[i] <- PA.KstochDeclSm.MRT[i] / PA.KstochDeclSm.VRT[i]
}
hist(PA.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.KstochDeclSm.iter <- dim(ZT.Nmat2)[1]
ZT.KstochDeclSm.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.KstochDeclSm.MRT <- ZT.KstochDeclSm.VRT <- ZT.KstochDeclSm.stat <- rep(NA,ZT.KstochDeclSm.iter)
for(i in 1:ZT.KstochDeclSm.iter) {
  ZT.KstochDeclSm.MRT[i] <- ZT.KstochDeclSm.RT[[i]]$mRT
  ZT.KstochDeclSm.VRT[i] <- ZT.KstochDeclSm.RT[[i]]$vRT
  ZT.KstochDeclSm.stat[i] <- ZT.KstochDeclSm.MRT[i] / ZT.KstochDeclSm.VRT[i]
}
hist(ZT.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.KstochDeclSm.iter <- dim(PH.Nmat2)[1]
PH.KstochDeclSm.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.KstochDeclSm.MRT <- PH.KstochDeclSm.VRT <- PH.KstochDeclSm.stat <- rep(NA,PH.KstochDeclSm.iter)
for(i in 1:PH.KstochDeclSm.iter) {
  PH.KstochDeclSm.MRT[i] <- PH.KstochDeclSm.RT[[i]]$mRT
  PH.KstochDeclSm.VRT[i] <- PH.KstochDeclSm.RT[[i]]$vRT
  PH.KstochDeclSm.stat[i] <- PH.KstochDeclSm.MRT[i] / PH.KstochDeclSm.VRT[i]
}
hist(PH.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.KstochDeclSm.iter <- dim(VU.Nmat2)[1]
VU.KstochDeclSm.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.KstochDeclSm.MRT <- VU.KstochDeclSm.VRT <- VU.KstochDeclSm.stat <- rep(NA,VU.KstochDeclSm.iter)
for(i in 1:VU.KstochDeclSm.iter) {
  VU.KstochDeclSm.MRT[i] <- VU.KstochDeclSm.RT[[i]]$mRT
  VU.KstochDeclSm.VRT[i] <- VU.KstochDeclSm.RT[[i]]$vRT
  VU.KstochDeclSm.stat[i] <- VU.KstochDeclSm.MRT[i] / VU.KstochDeclSm.VRT[i]
}
hist(VU.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.KstochDeclSm.iter <- dim(PG.Nmat2)[1]
PG.KstochDeclSm.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.KstochDeclSm.MRT <- PG.KstochDeclSm.VRT <- PG.KstochDeclSm.stat <- rep(NA,PG.KstochDeclSm.iter)
for(i in 1:PG.KstochDeclSm.iter) {
  PG.KstochDeclSm.MRT[i] <- PG.KstochDeclSm.RT[[i]]$mRT
  PG.KstochDeclSm.VRT[i] <- PG.KstochDeclSm.RT[[i]]$vRT
  PG.KstochDeclSm.stat[i] <- PG.KstochDeclSm.MRT[i] / PG.KstochDeclSm.VRT[i]
}
hist(PG.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.KstochDeclSm.iter <- dim(SS.Nmat2)[1]
SS.KstochDeclSm.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.KstochDeclSm.MRT <- SS.KstochDeclSm.VRT <- SS.KstochDeclSm.stat <- rep(NA,SS.KstochDeclSm.iter)
for(i in 1:SS.KstochDeclSm.iter) {
  SS.KstochDeclSm.MRT[i] <- SS.KstochDeclSm.RT[[i]]$mRT
  SS.KstochDeclSm.VRT[i] <- SS.KstochDeclSm.RT[[i]]$vRT
  SS.KstochDeclSm.stat[i] <- SS.KstochDeclSm.MRT[i] / SS.KstochDeclSm.VRT[i]
}
hist(SS.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.KstochDeclSm.iter <- dim(PT.Nmat2)[1]
PT.KstochDeclSm.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.KstochDeclSm.MRT <- PT.KstochDeclSm.VRT <- PT.KstochDeclSm.stat <- rep(NA,PT.KstochDeclSm.iter)
for(i in 1:PT.KstochDeclSm.iter) {
  PT.KstochDeclSm.MRT[i] <- PT.KstochDeclSm.RT[[i]]$mRT
  PT.KstochDeclSm.VRT[i] <- PT.KstochDeclSm.RT[[i]]$vRT
  PT.KstochDeclSm.stat[i] <- PT.KstochDeclSm.MRT[i] / PT.KstochDeclSm.VRT[i]
}
hist(PT.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.KstochDeclSm.iter <- dim(SO.Nmat2)[1]
SO.KstochDeclSm.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.KstochDeclSm.MRT <- SO.KstochDeclSm.VRT <- SO.KstochDeclSm.stat <- rep(NA,SO.KstochDeclSm.iter)
for(i in 1:SO.KstochDeclSm.iter) {
  SO.KstochDeclSm.MRT[i] <- SO.KstochDeclSm.RT[[i]]$mRT
  SO.KstochDeclSm.VRT[i] <- SO.KstochDeclSm.RT[[i]]$vRT
  SO.KstochDeclSm.stat[i] <- SO.KstochDeclSm.MRT[i] / SO.KstochDeclSm.VRT[i]
}
hist(SO.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.KstochDeclSm.iter <- dim(MN.Nmat2)[1]
MN.KstochDeclSm.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.KstochDeclSm.MRT <- MN.KstochDeclSm.VRT <- MN.KstochDeclSm.stat <- rep(NA,MN.KstochDeclSm.iter)
for(i in 1:MN.KstochDeclSm.iter) {
  MN.KstochDeclSm.MRT[i] <- MN.KstochDeclSm.RT[[i]]$mRT
  MN.KstochDeclSm.VRT[i] <- MN.KstochDeclSm.RT[[i]]$vRT
  MN.KstochDeclSm.stat[i] <- MN.KstochDeclSm.MRT[i] / MN.KstochDeclSm.VRT[i]
}
hist(MN.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.KstochDeclSm.iter <- dim(OR.Nmat2)[1]
OR.KstochDeclSm.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.KstochDeclSm.MRT <- OR.KstochDeclSm.VRT <- OR.KstochDeclSm.stat <- rep(NA,OR.KstochDeclSm.iter)
for(i in 1:OR.KstochDeclSm.iter) {
  OR.KstochDeclSm.MRT[i] <- OR.KstochDeclSm.RT[[i]]$mRT
  OR.KstochDeclSm.VRT[i] <- OR.KstochDeclSm.RT[[i]]$vRT
  OR.KstochDeclSm.stat[i] <- OR.KstochDeclSm.MRT[i] / OR.KstochDeclSm.VRT[i]
}
hist(OR.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.KstochDeclSm.iter <- dim(NR.Nmat2)[1]
NR.KstochDeclSm.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.KstochDeclSm.MRT <- NR.KstochDeclSm.VRT <- NR.KstochDeclSm.stat <- rep(NA,NR.KstochDeclSm.iter)
for(i in 1:NR.KstochDeclSm.iter) {
  NR.KstochDeclSm.MRT[i] <- NR.KstochDeclSm.RT[[i]]$mRT
  NR.KstochDeclSm.VRT[i] <- NR.KstochDeclSm.RT[[i]]$vRT
  NR.KstochDeclSm.stat[i] <- NR.KstochDeclSm.MRT[i] / NR.KstochDeclSm.VRT[i]
}
hist(NR.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.KstochDeclSm.iter <- dim(GN.Nmat2)[1]
GN.KstochDeclSm.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.KstochDeclSm.MRT <- GN.KstochDeclSm.VRT <- GN.KstochDeclSm.stat <- rep(NA,GN.KstochDeclSm.iter)
for(i in 1:GN.KstochDeclSm.iter) {
  GN.KstochDeclSm.MRT[i] <- GN.KstochDeclSm.RT[[i]]$mRT
  GN.KstochDeclSm.VRT[i] <- GN.KstochDeclSm.RT[[i]]$vRT
  GN.KstochDeclSm.stat[i] <- GN.KstochDeclSm.MRT[i] / GN.KstochDeclSm.VRT[i]
}
hist(GN.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.KstochDeclSm.iter <- dim(DN.Nmat2)[1]
DN.KstochDeclSm.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.KstochDeclSm.MRT <- DN.KstochDeclSm.VRT <- DN.KstochDeclSm.stat <- rep(NA,DN.KstochDeclSm.iter)
for(i in 1:DN.KstochDeclSm.iter) {
  DN.KstochDeclSm.MRT[i] <- DN.KstochDeclSm.RT[[i]]$mRT
  DN.KstochDeclSm.VRT[i] <- DN.KstochDeclSm.RT[[i]]$vRT
  DN.KstochDeclSm.stat[i] <- DN.KstochDeclSm.MRT[i] / DN.KstochDeclSm.VRT[i]
}
hist(DN.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.KstochDeclSm.iter <- dim(AL.Nmat2)[1]
AL.KstochDeclSm.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.KstochDeclSm.MRT <- AL.KstochDeclSm.VRT <- AL.KstochDeclSm.stat <- rep(NA,AL.KstochDeclSm.iter)
for(i in 1:AL.KstochDeclSm.iter) {
  AL.KstochDeclSm.MRT[i] <- AL.KstochDeclSm.RT[[i]]$mRT
  AL.KstochDeclSm.VRT[i] <- AL.KstochDeclSm.RT[[i]]$vRT
  AL.KstochDeclSm.stat[i] <- AL.KstochDeclSm.MRT[i] / AL.KstochDeclSm.VRT[i]
}
hist(AL.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.KstochDeclSm.iter <- dim(TC.Nmat2)[1]
TC.KstochDeclSm.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.KstochDeclSm.MRT <- TC.KstochDeclSm.VRT <- TC.KstochDeclSm.stat <- rep(NA,TC.KstochDeclSm.iter)
for(i in 1:TC.KstochDeclSm.iter) {
  TC.KstochDeclSm.MRT[i] <- TC.KstochDeclSm.RT[[i]]$mRT
  TC.KstochDeclSm.VRT[i] <- TC.KstochDeclSm.RT[[i]]$vRT
  TC.KstochDeclSm.stat[i] <- TC.KstochDeclSm.MRT[i] / TC.KstochDeclSm.VRT[i]
}
hist(TC.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.KstochDeclSm.iter <- dim(TH.Nmat2)[1]
TH.KstochDeclSm.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.KstochDeclSm.MRT <- TH.KstochDeclSm.VRT <- TH.KstochDeclSm.stat <- rep(NA,TH.KstochDeclSm.iter)
for(i in 1:TH.KstochDeclSm.iter) {
  TH.KstochDeclSm.MRT[i] <- TH.KstochDeclSm.RT[[i]]$mRT
  TH.KstochDeclSm.VRT[i] <- TH.KstochDeclSm.RT[[i]]$vRT
  TH.KstochDeclSm.stat[i] <- TH.KstochDeclSm.MRT[i] / TH.KstochDeclSm.VRT[i]
}
hist(TH.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.KstochDeclSm.iter <- dim(SH.Nmat2)[1]
SH.KstochDeclSm.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.KstochDeclSm.MRT <- SH.KstochDeclSm.VRT <- SH.KstochDeclSm.stat <- rep(NA,SH.KstochDeclSm.iter)
for(i in 1:SH.KstochDeclSm.iter) {
  SH.KstochDeclSm.MRT[i] <- SH.KstochDeclSm.RT[[i]]$mRT
  SH.KstochDeclSm.VRT[i] <- SH.KstochDeclSm.RT[[i]]$vRT
  SH.KstochDeclSm.stat[i] <- SH.KstochDeclSm.MRT[i] / SH.KstochDeclSm.VRT[i]
}
hist(SH.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.KstochDeclSm.iter <- dim(DM.Nmat2)[1]
DM.KstochDeclSm.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.KstochDeclSm.MRT <- DM.KstochDeclSm.VRT <- DM.KstochDeclSm.stat <- rep(NA,DM.KstochDeclSm.iter)
for(i in 1:DM.KstochDeclSm.iter) {
  DM.KstochDeclSm.MRT[i] <- DM.KstochDeclSm.RT[[i]]$mRT
  DM.KstochDeclSm.VRT[i] <- DM.KstochDeclSm.RT[[i]]$vRT
  DM.KstochDeclSm.stat[i] <- DM.KstochDeclSm.MRT[i] / DM.KstochDeclSm.VRT[i]
}
hist(DM.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.KstochDeclSm.iter <- dim(TA.Nmat2)[1]
TA.KstochDeclSm.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.KstochDeclSm.MRT <- TA.KstochDeclSm.VRT <- TA.KstochDeclSm.stat <- rep(NA,TA.KstochDeclSm.iter)
for(i in 1:TA.KstochDeclSm.iter) {
  TA.KstochDeclSm.MRT[i] <- TA.KstochDeclSm.RT[[i]]$mRT
  TA.KstochDeclSm.VRT[i] <- TA.KstochDeclSm.RT[[i]]$vRT
  TA.KstochDeclSm.stat[i] <- TA.KstochDeclSm.MRT[i] / TA.KstochDeclSm.VRT[i]
}
hist(TA.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.KstochDeclSm.iter <- dim(MR.Nmat2)[1]
MR.KstochDeclSm.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.KstochDeclSm.MRT <- MR.KstochDeclSm.VRT <- MR.KstochDeclSm.stat <- rep(NA,MR.KstochDeclSm.iter)
for(i in 1:MR.KstochDeclSm.iter) {
  MR.KstochDeclSm.MRT[i] <- MR.KstochDeclSm.RT[[i]]$mRT
  MR.KstochDeclSm.VRT[i] <- MR.KstochDeclSm.RT[[i]]$vRT
  MR.KstochDeclSm.stat[i] <- MR.KstochDeclSm.MRT[i] / MR.KstochDeclSm.VRT[i]
}
hist(MR.KstochDeclSm.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.KstochDeclSm.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.KstochDeclSm.stat.med <- median(DP.KstochDeclSm.stat, na.rm=T)
PA.KstochDeclSm.stat.med <- median(PA.KstochDeclSm.stat, na.rm=T)
ZT.KstochDeclSm.stat.med <- median(ZT.KstochDeclSm.stat, na.rm=T)
PH.KstochDeclSm.stat.med <- median(PH.KstochDeclSm.stat, na.rm=T)
VU.KstochDeclSm.stat.med <- median(VU.KstochDeclSm.stat, na.rm=T)
PG.KstochDeclSm.stat.med <- median(PG.KstochDeclSm.stat, na.rm=T)
SS.KstochDeclSm.stat.med <- median(SS.KstochDeclSm.stat, na.rm=T)
PT.KstochDeclSm.stat.med <- median(PT.KstochDeclSm.stat, na.rm=T)
SO.KstochDeclSm.stat.med <- median(SO.KstochDeclSm.stat, na.rm=T)
MN.KstochDeclSm.stat.med <- median(MN.KstochDeclSm.stat, na.rm=T)
OR.KstochDeclSm.stat.med <- median(OR.KstochDeclSm.stat, na.rm=T)
NR.KstochDeclSm.stat.med <- median(NR.KstochDeclSm.stat, na.rm=T)
GN.KstochDeclSm.stat.med <- median(GN.KstochDeclSm.stat, na.rm=T)
DN.KstochDeclSm.stat.med <- median(DN.KstochDeclSm.stat, na.rm=T)
AL.KstochDeclSm.stat.med <- median(AL.KstochDeclSm.stat, na.rm=T)
TC.KstochDeclSm.stat.med <- median(TC.KstochDeclSm.stat, na.rm=T)
TH.KstochDeclSm.stat.med <- median(TH.KstochDeclSm.stat, na.rm=T)
SH.KstochDeclSm.stat.med <- median(SH.KstochDeclSm.stat, na.rm=T)
DM.KstochDeclSm.stat.med <- median(DM.KstochDeclSm.stat, na.rm=T)
TA.KstochDeclSm.stat.med <- median(TA.KstochDeclSm.stat, na.rm=T)
MR.KstochDeclSm.stat.med <- median(MR.KstochDeclSm.stat, na.rm=T)

DP.KstochDeclSm.stat.up <- quantile(DP.KstochDeclSm.stat, probs=0.975, na.rm=T)
PA.KstochDeclSm.stat.up <- quantile(PA.KstochDeclSm.stat, probs=0.975, na.rm=T)
ZT.KstochDeclSm.stat.up <- quantile(ZT.KstochDeclSm.stat, probs=0.975, na.rm=T)
PH.KstochDeclSm.stat.up <- quantile(PH.KstochDeclSm.stat, probs=0.975, na.rm=T)
VU.KstochDeclSm.stat.up <- quantile(VU.KstochDeclSm.stat, probs=0.975, na.rm=T)
PG.KstochDeclSm.stat.up <- quantile(PG.KstochDeclSm.stat, probs=0.975, na.rm=T)
SS.KstochDeclSm.stat.up <- quantile(SS.KstochDeclSm.stat, probs=0.975, na.rm=T)
PT.KstochDeclSm.stat.up <- quantile(PT.KstochDeclSm.stat, probs=0.975, na.rm=T)
SO.KstochDeclSm.stat.up <- quantile(SO.KstochDeclSm.stat, probs=0.975, na.rm=T)
MN.KstochDeclSm.stat.up <- quantile(MN.KstochDeclSm.stat, probs=0.975, na.rm=T)
OR.KstochDeclSm.stat.up <- quantile(OR.KstochDeclSm.stat, probs=0.975, na.rm=T)
NR.KstochDeclSm.stat.up <- quantile(NR.KstochDeclSm.stat, probs=0.975, na.rm=T)
GN.KstochDeclSm.stat.up <- quantile(GN.KstochDeclSm.stat, probs=0.975, na.rm=T)
DN.KstochDeclSm.stat.up <- quantile(DN.KstochDeclSm.stat, probs=0.975, na.rm=T)
AL.KstochDeclSm.stat.up <- quantile(AL.KstochDeclSm.stat, probs=0.975, na.rm=T)
TC.KstochDeclSm.stat.up <- quantile(TC.KstochDeclSm.stat, probs=0.975, na.rm=T)
TH.KstochDeclSm.stat.up <- quantile(TH.KstochDeclSm.stat, probs=0.975, na.rm=T)
SH.KstochDeclSm.stat.up <- quantile(SH.KstochDeclSm.stat, probs=0.975, na.rm=T)
DM.KstochDeclSm.stat.up <- quantile(DM.KstochDeclSm.stat, probs=0.975, na.rm=T)
TA.KstochDeclSm.stat.up <- quantile(TA.KstochDeclSm.stat, probs=0.975, na.rm=T)
MR.KstochDeclSm.stat.up <- quantile(MR.KstochDeclSm.stat, probs=0.975, na.rm=T)

DP.KstochDeclSm.stat.lo <- quantile(DP.KstochDeclSm.stat, probs=0.025, na.rm=T)
PA.KstochDeclSm.stat.lo <- quantile(PA.KstochDeclSm.stat, probs=0.025, na.rm=T)
ZT.KstochDeclSm.stat.lo <- quantile(ZT.KstochDeclSm.stat, probs=0.025, na.rm=T)
PH.KstochDeclSm.stat.lo <- quantile(PH.KstochDeclSm.stat, probs=0.025, na.rm=T)
VU.KstochDeclSm.stat.lo <- quantile(VU.KstochDeclSm.stat, probs=0.025, na.rm=T)
PG.KstochDeclSm.stat.lo <- quantile(PG.KstochDeclSm.stat, probs=0.025, na.rm=T)
SS.KstochDeclSm.stat.lo <- quantile(SS.KstochDeclSm.stat, probs=0.025, na.rm=T)
PT.KstochDeclSm.stat.lo <- quantile(PT.KstochDeclSm.stat, probs=0.025, na.rm=T)
SO.KstochDeclSm.stat.lo <- quantile(SO.KstochDeclSm.stat, probs=0.025, na.rm=T)
MN.KstochDeclSm.stat.lo <- quantile(MN.KstochDeclSm.stat, probs=0.025, na.rm=T)
OR.KstochDeclSm.stat.lo <- quantile(OR.KstochDeclSm.stat, probs=0.025, na.rm=T)
NR.KstochDeclSm.stat.lo <- quantile(NR.KstochDeclSm.stat, probs=0.025, na.rm=T)
GN.KstochDeclSm.stat.lo <- quantile(GN.KstochDeclSm.stat, probs=0.025, na.rm=T)
DN.KstochDeclSm.stat.lo <- quantile(DN.KstochDeclSm.stat, probs=0.025, na.rm=T)
AL.KstochDeclSm.stat.lo <- quantile(AL.KstochDeclSm.stat, probs=0.025, na.rm=T)
TC.KstochDeclSm.stat.lo <- quantile(TC.KstochDeclSm.stat, probs=0.025, na.rm=T)
TH.KstochDeclSm.stat.lo <- quantile(TH.KstochDeclSm.stat, probs=0.025, na.rm=T)
SH.KstochDeclSm.stat.lo <- quantile(SH.KstochDeclSm.stat, probs=0.025, na.rm=T)
DM.KstochDeclSm.stat.lo <- quantile(DM.KstochDeclSm.stat, probs=0.025, na.rm=T)
TA.KstochDeclSm.stat.lo <- quantile(TA.KstochDeclSm.stat, probs=0.025, na.rm=T)
MR.KstochDeclSm.stat.lo <- quantile(MR.KstochDeclSm.stat, probs=0.025, na.rm=T)

KstochDeclSm.stat.med <- c(DP.KstochDeclSm.stat.med, PA.KstochDeclSm.stat.med, ZT.KstochDeclSm.stat.med, PH.KstochDeclSm.stat.med, VU.KstochDeclSm.stat.med, PG.KstochDeclSm.stat.med,
                           SS.KstochDeclSm.stat.med, PT.KstochDeclSm.stat.med, SO.KstochDeclSm.stat.med, MN.KstochDeclSm.stat.med, OR.KstochDeclSm.stat.med, NR.KstochDeclSm.stat.med,
                           GN.KstochDeclSm.stat.med, DN.KstochDeclSm.stat.med, AL.KstochDeclSm.stat.med, TC.KstochDeclSm.stat.med, TH.KstochDeclSm.stat.med, SH.KstochDeclSm.stat.med,
                           DM.KstochDeclSm.stat.med, TA.KstochDeclSm.stat.med, MR.KstochDeclSm.stat.med)
KstochDeclSm.stat.up <- c(DP.KstochDeclSm.stat.up, PA.KstochDeclSm.stat.up, ZT.KstochDeclSm.stat.up, PH.KstochDeclSm.stat.up, VU.KstochDeclSm.stat.up, PG.KstochDeclSm.stat.up,
                          SS.KstochDeclSm.stat.up, PT.KstochDeclSm.stat.up, SO.KstochDeclSm.stat.up, MN.KstochDeclSm.stat.up, OR.KstochDeclSm.stat.up, NR.KstochDeclSm.stat.up,
                          GN.KstochDeclSm.stat.up, DN.KstochDeclSm.stat.up, AL.KstochDeclSm.stat.up, TC.KstochDeclSm.stat.up, TH.KstochDeclSm.stat.up, SH.KstochDeclSm.stat.up,
                          DM.KstochDeclSm.stat.up, TA.KstochDeclSm.stat.up, MR.KstochDeclSm.stat.up)
KstochDeclSm.stat.lo <- c(DP.KstochDeclSm.stat.lo, PA.KstochDeclSm.stat.lo, ZT.KstochDeclSm.stat.lo, PH.KstochDeclSm.stat.lo, VU.KstochDeclSm.stat.lo, PG.KstochDeclSm.stat.lo,
                          SS.KstochDeclSm.stat.lo, PT.KstochDeclSm.stat.lo, SO.KstochDeclSm.stat.lo, MN.KstochDeclSm.stat.lo, OR.KstochDeclSm.stat.lo, NR.KstochDeclSm.stat.lo,
                          GN.KstochDeclSm.stat.lo, DN.KstochDeclSm.stat.lo, AL.KstochDeclSm.stat.lo, TC.KstochDeclSm.stat.lo, TH.KstochDeclSm.stat.lo, SH.KstochDeclSm.stat.lo,
                          DM.KstochDeclSm.stat.lo, TA.KstochDeclSm.stat.lo, MR.KstochDeclSm.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), KstochDeclSm.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
KstochDeclSm.dat <- data.frame(spp.mass.vec, KstochDeclSm.stat.med, KstochDeclSm.stat.up, KstochDeclSm.stat.lo)
colnames(KstochDeclSm.dat) <- c("M", "statM", "statUP", "statLO")
rownames(KstochDeclSm.dat) <- labs.vec

p <- ggplot(KstochDeclSm.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPKstochDeclSmstat <- data.frame(rep("DP",length(DP.KstochDeclSm.stat)), DP.KstochDeclSm.stat)
PAKstochDeclSmstat <- data.frame(rep("PA",length(PA.KstochDeclSm.stat)), PA.KstochDeclSm.stat)
ZTKstochDeclSmstat <- data.frame(rep("ZT",length(ZT.KstochDeclSm.stat)), ZT.KstochDeclSm.stat)
PHKstochDeclSmstat <- data.frame(rep("PH",length(PH.KstochDeclSm.stat)), PH.KstochDeclSm.stat)
VUKstochDeclSmstat <- data.frame(rep("VU",length(VU.KstochDeclSm.stat)), VU.KstochDeclSm.stat)
PGKstochDeclSmstat <- data.frame(rep("PG",length(PG.KstochDeclSm.stat)), PG.KstochDeclSm.stat)
SSKstochDeclSmstat <- data.frame(rep("SS",length(SS.KstochDeclSm.stat)), SS.KstochDeclSm.stat)
PTKstochDeclSmstat <- data.frame(rep("PT",length(PT.KstochDeclSm.stat)), PT.KstochDeclSm.stat)
SOKstochDeclSmstat <- data.frame(rep("SO",length(SO.KstochDeclSm.stat)), SO.KstochDeclSm.stat)
MNKstochDeclSmstat <- data.frame(rep("MN",length(MN.KstochDeclSm.stat)), MN.KstochDeclSm.stat)
ORKstochDeclSmstat <- data.frame(rep("OR",length(OR.KstochDeclSm.stat)), OR.KstochDeclSm.stat)
NRKstochDeclSmstat <- data.frame(rep("NR",length(NR.KstochDeclSm.stat)), NR.KstochDeclSm.stat)
GNKstochDeclSmstat <- data.frame(rep("GN",length(GN.KstochDeclSm.stat)), GN.KstochDeclSm.stat)
DNKstochDeclSmstat <- data.frame(rep("DN",length(DN.KstochDeclSm.stat)), DN.KstochDeclSm.stat)
ALKstochDeclSmstat <- data.frame(rep("AL",length(AL.KstochDeclSm.stat)), AL.KstochDeclSm.stat)
TCKstochDeclSmstat <- data.frame(rep("TC",length(TC.KstochDeclSm.stat)), TC.KstochDeclSm.stat)
THKstochDeclSmstat <- data.frame(rep("TH",length(TH.KstochDeclSm.stat)), TH.KstochDeclSm.stat)
SHKstochDeclSmstat <- data.frame(rep("SH",length(SH.KstochDeclSm.stat)), SH.KstochDeclSm.stat)
DMKstochDeclSmstat <- data.frame(rep("DM",length(DM.KstochDeclSm.stat)), DM.KstochDeclSm.stat)
TAKstochDeclSmstat <- data.frame(rep("TA",length(TA.KstochDeclSm.stat)), TA.KstochDeclSm.stat)
MRKstochDeclSmstat <- data.frame(rep("MR",length(MR.KstochDeclSm.stat)), MR.KstochDeclSm.stat)

colnames(DPKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(PAKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(ZTKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(PHKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(VUKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(PGKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(SSKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(PTKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(SOKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(MNKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(ORKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(NRKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(GNKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(DNKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(ALKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(TCKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(THKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(SHKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(DMKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(TAKstochDeclSmstat) <- c("SP","mRTvRT")
colnames(MRKstochDeclSmstat) <- c("SP","mRTvRT")

KstochDeclSmstat.dat <- rbind(DPKstochDeclSmstat, PAKstochDeclSmstat, ZTKstochDeclSmstat, PHKstochDeclSmstat, VUKstochDeclSmstat, PGKstochDeclSmstat, SSKstochDeclSmstat, PTKstochDeclSmstat,
                              SOKstochDeclSmstat, MNKstochDeclSmstat, ORKstochDeclSmstat, NRKstochDeclSmstat, GNKstochDeclSmstat, DNKstochDeclSmstat, ALKstochDeclSmstat, TCKstochDeclSmstat,
                              THKstochDeclSmstat, SHKstochDeclSmstat, DMKstochDeclSmstat, TAKstochDeclSmstat, MRKstochDeclSmstat)
KstochDeclSmstat.dat$SP <- as.factor(KstochDeclSmstat.dat$SP)
nspp <- length(table(KstochDeclSmstat.dat$SP))

mu <- KstochDeclSmstat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(KstochDeclSmstat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(KstochDeclSmstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(KstochDeclSmstat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(KstochDeclSmstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("KstochDeclSmstatdistrib.RData")
