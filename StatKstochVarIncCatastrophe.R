#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

##########################################################
##  K STOCHASTIC (5% SD initial size, Gaussian-resampled)
##  increasing to 10% over course of projection
##########################################################

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
load("DPKstochVarInc.RData")
load("PAKstochVarInc.RData")
load("ZTKstochVarInc.RData")
load("PHKstochVarInc.RData")
load("VUKstochVarInc.RData")
load("PGKstochVarInc.RData")
load("SSKstochVarInc.RData")
load("PTKstochVarInc.RData")
load("SOKstochVarInc.RData")
load("MNKstochVarInc.RData")
load("ORKstochVarInc.RData")
load("NRKstochVarInc.RData")
load("GNKstochVarInc.RData")
load("DNKstochVarInc.RData")
load("ALKstochVarInc.RData")
load("TCKstochVarInc.RData")
load("THKstochVarInc.RData")
load("SHKstochVarInc.RData")
load("DMKstochVarInc.RData")
load("TAKstochVarInc.RData")
load("MRKstochVarInc.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.KstochVarInc.iter <- dim(DP.Nmat2)[1]
DP.KstochVarInc.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.KstochVarInc.MRT <- DP.KstochVarInc.VRT <- DP.KstochVarInc.stat <- rep(NA,DP.KstochVarInc.iter)
for(i in 1:DP.KstochVarInc.iter) {
  DP.KstochVarInc.MRT[i] <- DP.KstochVarInc.RT[[i]]$mRT
  DP.KstochVarInc.VRT[i] <- DP.KstochVarInc.RT[[i]]$vRT
  DP.KstochVarInc.stat[i] <- DP.KstochVarInc.MRT[i] / DP.KstochVarInc.VRT[i]
}
hist(DP.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.KstochVarInc.iter <- dim(PA.Nmat2)[1]
PA.KstochVarInc.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.KstochVarInc.MRT <- PA.KstochVarInc.VRT <- PA.KstochVarInc.stat <- rep(NA,PA.KstochVarInc.iter)
for(i in 1:PA.KstochVarInc.iter) {
  PA.KstochVarInc.MRT[i] <- PA.KstochVarInc.RT[[i]]$mRT
  PA.KstochVarInc.VRT[i] <- PA.KstochVarInc.RT[[i]]$vRT
  PA.KstochVarInc.stat[i] <- PA.KstochVarInc.MRT[i] / PA.KstochVarInc.VRT[i]
}
hist(PA.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.KstochVarInc.iter <- dim(ZT.Nmat2)[1]
ZT.KstochVarInc.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.KstochVarInc.MRT <- ZT.KstochVarInc.VRT <- ZT.KstochVarInc.stat <- rep(NA,ZT.KstochVarInc.iter)
for(i in 1:ZT.KstochVarInc.iter) {
  ZT.KstochVarInc.MRT[i] <- ZT.KstochVarInc.RT[[i]]$mRT
  ZT.KstochVarInc.VRT[i] <- ZT.KstochVarInc.RT[[i]]$vRT
  ZT.KstochVarInc.stat[i] <- ZT.KstochVarInc.MRT[i] / ZT.KstochVarInc.VRT[i]
}
hist(ZT.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.KstochVarInc.iter <- dim(PH.Nmat2)[1]
PH.KstochVarInc.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.KstochVarInc.MRT <- PH.KstochVarInc.VRT <- PH.KstochVarInc.stat <- rep(NA,PH.KstochVarInc.iter)
for(i in 1:PH.KstochVarInc.iter) {
  PH.KstochVarInc.MRT[i] <- PH.KstochVarInc.RT[[i]]$mRT
  PH.KstochVarInc.VRT[i] <- PH.KstochVarInc.RT[[i]]$vRT
  PH.KstochVarInc.stat[i] <- PH.KstochVarInc.MRT[i] / PH.KstochVarInc.VRT[i]
}
hist(PH.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.KstochVarInc.iter <- dim(VU.Nmat2)[1]
VU.KstochVarInc.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.KstochVarInc.MRT <- VU.KstochVarInc.VRT <- VU.KstochVarInc.stat <- rep(NA,VU.KstochVarInc.iter)
for(i in 1:VU.KstochVarInc.iter) {
  VU.KstochVarInc.MRT[i] <- VU.KstochVarInc.RT[[i]]$mRT
  VU.KstochVarInc.VRT[i] <- VU.KstochVarInc.RT[[i]]$vRT
  VU.KstochVarInc.stat[i] <- VU.KstochVarInc.MRT[i] / VU.KstochVarInc.VRT[i]
}
hist(VU.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.KstochVarInc.iter <- dim(PG.Nmat2)[1]
PG.KstochVarInc.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.KstochVarInc.MRT <- PG.KstochVarInc.VRT <- PG.KstochVarInc.stat <- rep(NA,PG.KstochVarInc.iter)
for(i in 1:PG.KstochVarInc.iter) {
  PG.KstochVarInc.MRT[i] <- PG.KstochVarInc.RT[[i]]$mRT
  PG.KstochVarInc.VRT[i] <- PG.KstochVarInc.RT[[i]]$vRT
  PG.KstochVarInc.stat[i] <- PG.KstochVarInc.MRT[i] / PG.KstochVarInc.VRT[i]
}
hist(PG.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.KstochVarInc.iter <- dim(SS.Nmat2)[1]
SS.KstochVarInc.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.KstochVarInc.MRT <- SS.KstochVarInc.VRT <- SS.KstochVarInc.stat <- rep(NA,SS.KstochVarInc.iter)
for(i in 1:SS.KstochVarInc.iter) {
  SS.KstochVarInc.MRT[i] <- SS.KstochVarInc.RT[[i]]$mRT
  SS.KstochVarInc.VRT[i] <- SS.KstochVarInc.RT[[i]]$vRT
  SS.KstochVarInc.stat[i] <- SS.KstochVarInc.MRT[i] / SS.KstochVarInc.VRT[i]
}
hist(SS.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.KstochVarInc.iter <- dim(PT.Nmat2)[1]
PT.KstochVarInc.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.KstochVarInc.MRT <- PT.KstochVarInc.VRT <- PT.KstochVarInc.stat <- rep(NA,PT.KstochVarInc.iter)
for(i in 1:PT.KstochVarInc.iter) {
  PT.KstochVarInc.MRT[i] <- PT.KstochVarInc.RT[[i]]$mRT
  PT.KstochVarInc.VRT[i] <- PT.KstochVarInc.RT[[i]]$vRT
  PT.KstochVarInc.stat[i] <- PT.KstochVarInc.MRT[i] / PT.KstochVarInc.VRT[i]
}
hist(PT.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.KstochVarInc.iter <- dim(SO.Nmat2)[1]
SO.KstochVarInc.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.KstochVarInc.MRT <- SO.KstochVarInc.VRT <- SO.KstochVarInc.stat <- rep(NA,SO.KstochVarInc.iter)
for(i in 1:SO.KstochVarInc.iter) {
  SO.KstochVarInc.MRT[i] <- SO.KstochVarInc.RT[[i]]$mRT
  SO.KstochVarInc.VRT[i] <- SO.KstochVarInc.RT[[i]]$vRT
  SO.KstochVarInc.stat[i] <- SO.KstochVarInc.MRT[i] / SO.KstochVarInc.VRT[i]
}
hist(SO.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.KstochVarInc.iter <- dim(MN.Nmat2)[1]
MN.KstochVarInc.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.KstochVarInc.MRT <- MN.KstochVarInc.VRT <- MN.KstochVarInc.stat <- rep(NA,MN.KstochVarInc.iter)
for(i in 1:MN.KstochVarInc.iter) {
  MN.KstochVarInc.MRT[i] <- MN.KstochVarInc.RT[[i]]$mRT
  MN.KstochVarInc.VRT[i] <- MN.KstochVarInc.RT[[i]]$vRT
  MN.KstochVarInc.stat[i] <- MN.KstochVarInc.MRT[i] / MN.KstochVarInc.VRT[i]
}
hist(MN.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.KstochVarInc.iter <- dim(OR.Nmat2)[1]
OR.KstochVarInc.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.KstochVarInc.MRT <- OR.KstochVarInc.VRT <- OR.KstochVarInc.stat <- rep(NA,OR.KstochVarInc.iter)
for(i in 1:OR.KstochVarInc.iter) {
  OR.KstochVarInc.MRT[i] <- OR.KstochVarInc.RT[[i]]$mRT
  OR.KstochVarInc.VRT[i] <- OR.KstochVarInc.RT[[i]]$vRT
  OR.KstochVarInc.stat[i] <- OR.KstochVarInc.MRT[i] / OR.KstochVarInc.VRT[i]
}
hist(OR.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.KstochVarInc.iter <- dim(NR.Nmat2)[1]
NR.KstochVarInc.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.KstochVarInc.MRT <- NR.KstochVarInc.VRT <- NR.KstochVarInc.stat <- rep(NA,NR.KstochVarInc.iter)
for(i in 1:NR.KstochVarInc.iter) {
  NR.KstochVarInc.MRT[i] <- NR.KstochVarInc.RT[[i]]$mRT
  NR.KstochVarInc.VRT[i] <- NR.KstochVarInc.RT[[i]]$vRT
  NR.KstochVarInc.stat[i] <- NR.KstochVarInc.MRT[i] / NR.KstochVarInc.VRT[i]
}
hist(NR.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.KstochVarInc.iter <- dim(GN.Nmat2)[1]
GN.KstochVarInc.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.KstochVarInc.MRT <- GN.KstochVarInc.VRT <- GN.KstochVarInc.stat <- rep(NA,GN.KstochVarInc.iter)
for(i in 1:GN.KstochVarInc.iter) {
  GN.KstochVarInc.MRT[i] <- GN.KstochVarInc.RT[[i]]$mRT
  GN.KstochVarInc.VRT[i] <- GN.KstochVarInc.RT[[i]]$vRT
  GN.KstochVarInc.stat[i] <- GN.KstochVarInc.MRT[i] / GN.KstochVarInc.VRT[i]
}
hist(GN.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.KstochVarInc.iter <- dim(DN.Nmat2)[1]
DN.KstochVarInc.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.KstochVarInc.MRT <- DN.KstochVarInc.VRT <- DN.KstochVarInc.stat <- rep(NA,DN.KstochVarInc.iter)
for(i in 1:DN.KstochVarInc.iter) {
  DN.KstochVarInc.MRT[i] <- DN.KstochVarInc.RT[[i]]$mRT
  DN.KstochVarInc.VRT[i] <- DN.KstochVarInc.RT[[i]]$vRT
  DN.KstochVarInc.stat[i] <- DN.KstochVarInc.MRT[i] / DN.KstochVarInc.VRT[i]
}
hist(DN.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.KstochVarInc.iter <- dim(AL.Nmat2)[1]
AL.KstochVarInc.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.KstochVarInc.MRT <- AL.KstochVarInc.VRT <- AL.KstochVarInc.stat <- rep(NA,AL.KstochVarInc.iter)
for(i in 1:AL.KstochVarInc.iter) {
  AL.KstochVarInc.MRT[i] <- AL.KstochVarInc.RT[[i]]$mRT
  AL.KstochVarInc.VRT[i] <- AL.KstochVarInc.RT[[i]]$vRT
  AL.KstochVarInc.stat[i] <- AL.KstochVarInc.MRT[i] / AL.KstochVarInc.VRT[i]
}
hist(AL.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.KstochVarInc.iter <- dim(TC.Nmat2)[1]
TC.KstochVarInc.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.KstochVarInc.MRT <- TC.KstochVarInc.VRT <- TC.KstochVarInc.stat <- rep(NA,TC.KstochVarInc.iter)
for(i in 1:TC.KstochVarInc.iter) {
  TC.KstochVarInc.MRT[i] <- TC.KstochVarInc.RT[[i]]$mRT
  TC.KstochVarInc.VRT[i] <- TC.KstochVarInc.RT[[i]]$vRT
  TC.KstochVarInc.stat[i] <- TC.KstochVarInc.MRT[i] / TC.KstochVarInc.VRT[i]
}
hist(TC.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.KstochVarInc.iter <- dim(TH.Nmat2)[1]
TH.KstochVarInc.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.KstochVarInc.MRT <- TH.KstochVarInc.VRT <- TH.KstochVarInc.stat <- rep(NA,TH.KstochVarInc.iter)
for(i in 1:TH.KstochVarInc.iter) {
  TH.KstochVarInc.MRT[i] <- TH.KstochVarInc.RT[[i]]$mRT
  TH.KstochVarInc.VRT[i] <- TH.KstochVarInc.RT[[i]]$vRT
  TH.KstochVarInc.stat[i] <- TH.KstochVarInc.MRT[i] / TH.KstochVarInc.VRT[i]
}
hist(TH.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.KstochVarInc.iter <- dim(SH.Nmat2)[1]
SH.KstochVarInc.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.KstochVarInc.MRT <- SH.KstochVarInc.VRT <- SH.KstochVarInc.stat <- rep(NA,SH.KstochVarInc.iter)
for(i in 1:SH.KstochVarInc.iter) {
  SH.KstochVarInc.MRT[i] <- SH.KstochVarInc.RT[[i]]$mRT
  SH.KstochVarInc.VRT[i] <- SH.KstochVarInc.RT[[i]]$vRT
  SH.KstochVarInc.stat[i] <- SH.KstochVarInc.MRT[i] / SH.KstochVarInc.VRT[i]
}
hist(SH.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.KstochVarInc.iter <- dim(DM.Nmat2)[1]
DM.KstochVarInc.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.KstochVarInc.MRT <- DM.KstochVarInc.VRT <- DM.KstochVarInc.stat <- rep(NA,DM.KstochVarInc.iter)
for(i in 1:DM.KstochVarInc.iter) {
  DM.KstochVarInc.MRT[i] <- DM.KstochVarInc.RT[[i]]$mRT
  DM.KstochVarInc.VRT[i] <- DM.KstochVarInc.RT[[i]]$vRT
  DM.KstochVarInc.stat[i] <- DM.KstochVarInc.MRT[i] / DM.KstochVarInc.VRT[i]
}
hist(DM.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.KstochVarInc.iter <- dim(TA.Nmat2)[1]
TA.KstochVarInc.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.KstochVarInc.MRT <- TA.KstochVarInc.VRT <- TA.KstochVarInc.stat <- rep(NA,TA.KstochVarInc.iter)
for(i in 1:TA.KstochVarInc.iter) {
  TA.KstochVarInc.MRT[i] <- TA.KstochVarInc.RT[[i]]$mRT
  TA.KstochVarInc.VRT[i] <- TA.KstochVarInc.RT[[i]]$vRT
  TA.KstochVarInc.stat[i] <- TA.KstochVarInc.MRT[i] / TA.KstochVarInc.VRT[i]
}
hist(TA.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.KstochVarInc.iter <- dim(MR.Nmat2)[1]
MR.KstochVarInc.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.KstochVarInc.MRT <- MR.KstochVarInc.VRT <- MR.KstochVarInc.stat <- rep(NA,MR.KstochVarInc.iter)
for(i in 1:MR.KstochVarInc.iter) {
  MR.KstochVarInc.MRT[i] <- MR.KstochVarInc.RT[[i]]$mRT
  MR.KstochVarInc.VRT[i] <- MR.KstochVarInc.RT[[i]]$vRT
  MR.KstochVarInc.stat[i] <- MR.KstochVarInc.MRT[i] / MR.KstochVarInc.VRT[i]
}
hist(MR.KstochVarInc.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.KstochVarInc.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.KstochVarInc.stat.med <- median(DP.KstochVarInc.stat, na.rm=T)
PA.KstochVarInc.stat.med <- median(PA.KstochVarInc.stat, na.rm=T)
ZT.KstochVarInc.stat.med <- median(ZT.KstochVarInc.stat, na.rm=T)
PH.KstochVarInc.stat.med <- median(PH.KstochVarInc.stat, na.rm=T)
VU.KstochVarInc.stat.med <- median(VU.KstochVarInc.stat, na.rm=T)
PG.KstochVarInc.stat.med <- median(PG.KstochVarInc.stat, na.rm=T)
SS.KstochVarInc.stat.med <- median(SS.KstochVarInc.stat, na.rm=T)
PT.KstochVarInc.stat.med <- median(PT.KstochVarInc.stat, na.rm=T)
SO.KstochVarInc.stat.med <- median(SO.KstochVarInc.stat, na.rm=T)
MN.KstochVarInc.stat.med <- median(MN.KstochVarInc.stat, na.rm=T)
OR.KstochVarInc.stat.med <- median(OR.KstochVarInc.stat, na.rm=T)
NR.KstochVarInc.stat.med <- median(NR.KstochVarInc.stat, na.rm=T)
GN.KstochVarInc.stat.med <- median(GN.KstochVarInc.stat, na.rm=T)
DN.KstochVarInc.stat.med <- median(DN.KstochVarInc.stat, na.rm=T)
AL.KstochVarInc.stat.med <- median(AL.KstochVarInc.stat, na.rm=T)
TC.KstochVarInc.stat.med <- median(TC.KstochVarInc.stat, na.rm=T)
TH.KstochVarInc.stat.med <- median(TH.KstochVarInc.stat, na.rm=T)
SH.KstochVarInc.stat.med <- median(SH.KstochVarInc.stat, na.rm=T)
DM.KstochVarInc.stat.med <- median(DM.KstochVarInc.stat, na.rm=T)
TA.KstochVarInc.stat.med <- median(TA.KstochVarInc.stat, na.rm=T)
MR.KstochVarInc.stat.med <- median(MR.KstochVarInc.stat, na.rm=T)

DP.KstochVarInc.stat.up <- quantile(DP.KstochVarInc.stat, probs=0.975, na.rm=T)
PA.KstochVarInc.stat.up <- quantile(PA.KstochVarInc.stat, probs=0.975, na.rm=T)
ZT.KstochVarInc.stat.up <- quantile(ZT.KstochVarInc.stat, probs=0.975, na.rm=T)
PH.KstochVarInc.stat.up <- quantile(PH.KstochVarInc.stat, probs=0.975, na.rm=T)
VU.KstochVarInc.stat.up <- quantile(VU.KstochVarInc.stat, probs=0.975, na.rm=T)
PG.KstochVarInc.stat.up <- quantile(PG.KstochVarInc.stat, probs=0.975, na.rm=T)
SS.KstochVarInc.stat.up <- quantile(SS.KstochVarInc.stat, probs=0.975, na.rm=T)
PT.KstochVarInc.stat.up <- quantile(PT.KstochVarInc.stat, probs=0.975, na.rm=T)
SO.KstochVarInc.stat.up <- quantile(SO.KstochVarInc.stat, probs=0.975, na.rm=T)
MN.KstochVarInc.stat.up <- quantile(MN.KstochVarInc.stat, probs=0.975, na.rm=T)
OR.KstochVarInc.stat.up <- quantile(OR.KstochVarInc.stat, probs=0.975, na.rm=T)
NR.KstochVarInc.stat.up <- quantile(NR.KstochVarInc.stat, probs=0.975, na.rm=T)
GN.KstochVarInc.stat.up <- quantile(GN.KstochVarInc.stat, probs=0.975, na.rm=T)
DN.KstochVarInc.stat.up <- quantile(DN.KstochVarInc.stat, probs=0.975, na.rm=T)
AL.KstochVarInc.stat.up <- quantile(AL.KstochVarInc.stat, probs=0.975, na.rm=T)
TC.KstochVarInc.stat.up <- quantile(TC.KstochVarInc.stat, probs=0.975, na.rm=T)
TH.KstochVarInc.stat.up <- quantile(TH.KstochVarInc.stat, probs=0.975, na.rm=T)
SH.KstochVarInc.stat.up <- quantile(SH.KstochVarInc.stat, probs=0.975, na.rm=T)
DM.KstochVarInc.stat.up <- quantile(DM.KstochVarInc.stat, probs=0.975, na.rm=T)
TA.KstochVarInc.stat.up <- quantile(TA.KstochVarInc.stat, probs=0.975, na.rm=T)
MR.KstochVarInc.stat.up <- quantile(MR.KstochVarInc.stat, probs=0.975, na.rm=T)

DP.KstochVarInc.stat.lo <- quantile(DP.KstochVarInc.stat, probs=0.025, na.rm=T)
PA.KstochVarInc.stat.lo <- quantile(PA.KstochVarInc.stat, probs=0.025, na.rm=T)
ZT.KstochVarInc.stat.lo <- quantile(ZT.KstochVarInc.stat, probs=0.025, na.rm=T)
PH.KstochVarInc.stat.lo <- quantile(PH.KstochVarInc.stat, probs=0.025, na.rm=T)
VU.KstochVarInc.stat.lo <- quantile(VU.KstochVarInc.stat, probs=0.025, na.rm=T)
PG.KstochVarInc.stat.lo <- quantile(PG.KstochVarInc.stat, probs=0.025, na.rm=T)
SS.KstochVarInc.stat.lo <- quantile(SS.KstochVarInc.stat, probs=0.025, na.rm=T)
PT.KstochVarInc.stat.lo <- quantile(PT.KstochVarInc.stat, probs=0.025, na.rm=T)
SO.KstochVarInc.stat.lo <- quantile(SO.KstochVarInc.stat, probs=0.025, na.rm=T)
MN.KstochVarInc.stat.lo <- quantile(MN.KstochVarInc.stat, probs=0.025, na.rm=T)
OR.KstochVarInc.stat.lo <- quantile(OR.KstochVarInc.stat, probs=0.025, na.rm=T)
NR.KstochVarInc.stat.lo <- quantile(NR.KstochVarInc.stat, probs=0.025, na.rm=T)
GN.KstochVarInc.stat.lo <- quantile(GN.KstochVarInc.stat, probs=0.025, na.rm=T)
DN.KstochVarInc.stat.lo <- quantile(DN.KstochVarInc.stat, probs=0.025, na.rm=T)
AL.KstochVarInc.stat.lo <- quantile(AL.KstochVarInc.stat, probs=0.025, na.rm=T)
TC.KstochVarInc.stat.lo <- quantile(TC.KstochVarInc.stat, probs=0.025, na.rm=T)
TH.KstochVarInc.stat.lo <- quantile(TH.KstochVarInc.stat, probs=0.025, na.rm=T)
SH.KstochVarInc.stat.lo <- quantile(SH.KstochVarInc.stat, probs=0.025, na.rm=T)
DM.KstochVarInc.stat.lo <- quantile(DM.KstochVarInc.stat, probs=0.025, na.rm=T)
TA.KstochVarInc.stat.lo <- quantile(TA.KstochVarInc.stat, probs=0.025, na.rm=T)
MR.KstochVarInc.stat.lo <- quantile(MR.KstochVarInc.stat, probs=0.025, na.rm=T)

KstochVarInc.stat.med <- c(DP.KstochVarInc.stat.med, PA.KstochVarInc.stat.med, ZT.KstochVarInc.stat.med, PH.KstochVarInc.stat.med, VU.KstochVarInc.stat.med, PG.KstochVarInc.stat.med,
                           SS.KstochVarInc.stat.med, PT.KstochVarInc.stat.med, SO.KstochVarInc.stat.med, MN.KstochVarInc.stat.med, OR.KstochVarInc.stat.med, NR.KstochVarInc.stat.med,
                           GN.KstochVarInc.stat.med, DN.KstochVarInc.stat.med, AL.KstochVarInc.stat.med, TC.KstochVarInc.stat.med, TH.KstochVarInc.stat.med, SH.KstochVarInc.stat.med,
                           DM.KstochVarInc.stat.med, TA.KstochVarInc.stat.med, MR.KstochVarInc.stat.med)
KstochVarInc.stat.up <- c(DP.KstochVarInc.stat.up, PA.KstochVarInc.stat.up, ZT.KstochVarInc.stat.up, PH.KstochVarInc.stat.up, VU.KstochVarInc.stat.up, PG.KstochVarInc.stat.up,
                          SS.KstochVarInc.stat.up, PT.KstochVarInc.stat.up, SO.KstochVarInc.stat.up, MN.KstochVarInc.stat.up, OR.KstochVarInc.stat.up, NR.KstochVarInc.stat.up,
                          GN.KstochVarInc.stat.up, DN.KstochVarInc.stat.up, AL.KstochVarInc.stat.up, TC.KstochVarInc.stat.up, TH.KstochVarInc.stat.up, SH.KstochVarInc.stat.up,
                          DM.KstochVarInc.stat.up, TA.KstochVarInc.stat.up, MR.KstochVarInc.stat.up)
KstochVarInc.stat.lo <- c(DP.KstochVarInc.stat.lo, PA.KstochVarInc.stat.lo, ZT.KstochVarInc.stat.lo, PH.KstochVarInc.stat.lo, VU.KstochVarInc.stat.lo, PG.KstochVarInc.stat.lo,
                          SS.KstochVarInc.stat.lo, PT.KstochVarInc.stat.lo, SO.KstochVarInc.stat.lo, MN.KstochVarInc.stat.lo, OR.KstochVarInc.stat.lo, NR.KstochVarInc.stat.lo,
                          GN.KstochVarInc.stat.lo, DN.KstochVarInc.stat.lo, AL.KstochVarInc.stat.lo, TC.KstochVarInc.stat.lo, TH.KstochVarInc.stat.lo, SH.KstochVarInc.stat.lo,
                          DM.KstochVarInc.stat.lo, TA.KstochVarInc.stat.lo, MR.KstochVarInc.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), KstochVarInc.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
KstochVarInc.dat <- data.frame(spp.mass.vec, KstochVarInc.stat.med, KstochVarInc.stat.up, KstochVarInc.stat.lo)
colnames(KstochVarInc.dat) <- c("M", "statM", "statUP", "statLO")
rownames(KstochVarInc.dat) <- labs.vec

p <- ggplot(KstochVarInc.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPKstochVarIncstat <- data.frame(rep("DP",length(DP.KstochVarInc.stat)), DP.KstochVarInc.stat)
PAKstochVarIncstat <- data.frame(rep("PA",length(PA.KstochVarInc.stat)), PA.KstochVarInc.stat)
ZTKstochVarIncstat <- data.frame(rep("ZT",length(ZT.KstochVarInc.stat)), ZT.KstochVarInc.stat)
PHKstochVarIncstat <- data.frame(rep("PH",length(PH.KstochVarInc.stat)), PH.KstochVarInc.stat)
VUKstochVarIncstat <- data.frame(rep("VU",length(VU.KstochVarInc.stat)), VU.KstochVarInc.stat)
PGKstochVarIncstat <- data.frame(rep("PG",length(PG.KstochVarInc.stat)), PG.KstochVarInc.stat)
SSKstochVarIncstat <- data.frame(rep("SS",length(SS.KstochVarInc.stat)), SS.KstochVarInc.stat)
PTKstochVarIncstat <- data.frame(rep("PT",length(PT.KstochVarInc.stat)), PT.KstochVarInc.stat)
SOKstochVarIncstat <- data.frame(rep("SO",length(SO.KstochVarInc.stat)), SO.KstochVarInc.stat)
MNKstochVarIncstat <- data.frame(rep("MN",length(MN.KstochVarInc.stat)), MN.KstochVarInc.stat)
ORKstochVarIncstat <- data.frame(rep("OR",length(OR.KstochVarInc.stat)), OR.KstochVarInc.stat)
NRKstochVarIncstat <- data.frame(rep("NR",length(NR.KstochVarInc.stat)), NR.KstochVarInc.stat)
GNKstochVarIncstat <- data.frame(rep("GN",length(GN.KstochVarInc.stat)), GN.KstochVarInc.stat)
DNKstochVarIncstat <- data.frame(rep("DN",length(DN.KstochVarInc.stat)), DN.KstochVarInc.stat)
ALKstochVarIncstat <- data.frame(rep("AL",length(AL.KstochVarInc.stat)), AL.KstochVarInc.stat)
TCKstochVarIncstat <- data.frame(rep("TC",length(TC.KstochVarInc.stat)), TC.KstochVarInc.stat)
THKstochVarIncstat <- data.frame(rep("TH",length(TH.KstochVarInc.stat)), TH.KstochVarInc.stat)
SHKstochVarIncstat <- data.frame(rep("SH",length(SH.KstochVarInc.stat)), SH.KstochVarInc.stat)
DMKstochVarIncstat <- data.frame(rep("DM",length(DM.KstochVarInc.stat)), DM.KstochVarInc.stat)
TAKstochVarIncstat <- data.frame(rep("TA",length(TA.KstochVarInc.stat)), TA.KstochVarInc.stat)
MRKstochVarIncstat <- data.frame(rep("MR",length(MR.KstochVarInc.stat)), MR.KstochVarInc.stat)

colnames(DPKstochVarIncstat) <- c("SP","mRTvRT")
colnames(PAKstochVarIncstat) <- c("SP","mRTvRT")
colnames(ZTKstochVarIncstat) <- c("SP","mRTvRT")
colnames(PHKstochVarIncstat) <- c("SP","mRTvRT")
colnames(VUKstochVarIncstat) <- c("SP","mRTvRT")
colnames(PGKstochVarIncstat) <- c("SP","mRTvRT")
colnames(SSKstochVarIncstat) <- c("SP","mRTvRT")
colnames(PTKstochVarIncstat) <- c("SP","mRTvRT")
colnames(SOKstochVarIncstat) <- c("SP","mRTvRT")
colnames(MNKstochVarIncstat) <- c("SP","mRTvRT")
colnames(ORKstochVarIncstat) <- c("SP","mRTvRT")
colnames(NRKstochVarIncstat) <- c("SP","mRTvRT")
colnames(GNKstochVarIncstat) <- c("SP","mRTvRT")
colnames(DNKstochVarIncstat) <- c("SP","mRTvRT")
colnames(ALKstochVarIncstat) <- c("SP","mRTvRT")
colnames(TCKstochVarIncstat) <- c("SP","mRTvRT")
colnames(THKstochVarIncstat) <- c("SP","mRTvRT")
colnames(SHKstochVarIncstat) <- c("SP","mRTvRT")
colnames(DMKstochVarIncstat) <- c("SP","mRTvRT")
colnames(TAKstochVarIncstat) <- c("SP","mRTvRT")
colnames(MRKstochVarIncstat) <- c("SP","mRTvRT")

KstochVarIncstat.dat <- rbind(DPKstochVarIncstat, PAKstochVarIncstat, ZTKstochVarIncstat, PHKstochVarIncstat, VUKstochVarIncstat, PGKstochVarIncstat, SSKstochVarIncstat, PTKstochVarIncstat,
                              SOKstochVarIncstat, MNKstochVarIncstat, ORKstochVarIncstat, NRKstochVarIncstat, GNKstochVarIncstat, DNKstochVarIncstat, ALKstochVarIncstat, TCKstochVarIncstat,
                              THKstochVarIncstat, SHKstochVarIncstat, DMKstochVarIncstat, TAKstochVarIncstat, MRKstochVarIncstat)
KstochVarIncstat.dat$SP <- as.factor(KstochVarIncstat.dat$SP)
nspp <- length(table(KstochVarIncstat.dat$SP))

mu <- KstochVarIncstat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(KstochVarIncstat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(KstochVarIncstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(KstochVarIncstat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(KstochVarIncstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("KstochVarIncstatdistrib.RData")

