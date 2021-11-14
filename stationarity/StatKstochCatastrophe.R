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
load("DPKstoch.RData")
load("PAKstoch.RData")
load("ZTKstoch.RData")
load("PHKstoch.RData")
load("VUKstoch.RData")
load("PGKstoch.RData")
load("SSKstoch.RData")
load("PTKstoch.RData")
load("SOKstoch.RData")
load("MNKstoch.RData")
load("ORKstoch.RData")
load("NRKstoch.RData")
load("GNKstoch.RData")
load("DNKstoch.RData")
load("ALKstoch.RData")
load("TCKstoch.RData")
load("THKstoch.RData")
load("SHKstoch.RData")
load("DMKstoch.RData")
load("TAKstoch.RData")
load("MRKstoch.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.Kstoch.iter <- dim(DP.Nmat2)[1]
DP.Kstoch.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.Kstoch.MRT <- DP.Kstoch.VRT <- DP.Kstoch.stat <- rep(NA,DP.Kstoch.iter)
for(i in 1:DP.Kstoch.iter) {
  DP.Kstoch.MRT[i] <- DP.Kstoch.RT[[i]]$mRT
  DP.Kstoch.VRT[i] <- DP.Kstoch.RT[[i]]$vRT
  DP.Kstoch.stat[i] <- DP.Kstoch.MRT[i] / DP.Kstoch.VRT[i]
}
hist(DP.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.Kstoch.iter <- dim(PA.Nmat2)[1]
PA.Kstoch.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.Kstoch.MRT <- PA.Kstoch.VRT <- PA.Kstoch.stat <- rep(NA,PA.Kstoch.iter)
for(i in 1:PA.Kstoch.iter) {
  PA.Kstoch.MRT[i] <- PA.Kstoch.RT[[i]]$mRT
  PA.Kstoch.VRT[i] <- PA.Kstoch.RT[[i]]$vRT
  PA.Kstoch.stat[i] <- PA.Kstoch.MRT[i] / PA.Kstoch.VRT[i]
}
hist(PA.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.Kstoch.iter <- dim(ZT.Nmat2)[1]
ZT.Kstoch.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.Kstoch.MRT <- ZT.Kstoch.VRT <- ZT.Kstoch.stat <- rep(NA,ZT.Kstoch.iter)
for(i in 1:ZT.Kstoch.iter) {
  ZT.Kstoch.MRT[i] <- ZT.Kstoch.RT[[i]]$mRT
  ZT.Kstoch.VRT[i] <- ZT.Kstoch.RT[[i]]$vRT
  ZT.Kstoch.stat[i] <- ZT.Kstoch.MRT[i] / ZT.Kstoch.VRT[i]
}
hist(ZT.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.Kstoch.iter <- dim(PH.Nmat2)[1]
PH.Kstoch.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.Kstoch.MRT <- PH.Kstoch.VRT <- PH.Kstoch.stat <- rep(NA,PH.Kstoch.iter)
for(i in 1:PH.Kstoch.iter) {
  PH.Kstoch.MRT[i] <- PH.Kstoch.RT[[i]]$mRT
  PH.Kstoch.VRT[i] <- PH.Kstoch.RT[[i]]$vRT
  PH.Kstoch.stat[i] <- PH.Kstoch.MRT[i] / PH.Kstoch.VRT[i]
}
hist(PH.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.Kstoch.iter <- dim(VU.Nmat2)[1]
VU.Kstoch.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.Kstoch.MRT <- VU.Kstoch.VRT <- VU.Kstoch.stat <- rep(NA,VU.Kstoch.iter)
for(i in 1:VU.Kstoch.iter) {
  VU.Kstoch.MRT[i] <- VU.Kstoch.RT[[i]]$mRT
  VU.Kstoch.VRT[i] <- VU.Kstoch.RT[[i]]$vRT
  VU.Kstoch.stat[i] <- VU.Kstoch.MRT[i] / VU.Kstoch.VRT[i]
}
hist(VU.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.Kstoch.iter <- dim(PG.Nmat2)[1]
PG.Kstoch.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.Kstoch.MRT <- PG.Kstoch.VRT <- PG.Kstoch.stat <- rep(NA,PG.Kstoch.iter)
for(i in 1:PG.Kstoch.iter) {
  PG.Kstoch.MRT[i] <- PG.Kstoch.RT[[i]]$mRT
  PG.Kstoch.VRT[i] <- PG.Kstoch.RT[[i]]$vRT
  PG.Kstoch.stat[i] <- PG.Kstoch.MRT[i] / PG.Kstoch.VRT[i]
}
hist(PG.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.Kstoch.iter <- dim(SS.Nmat2)[1]
SS.Kstoch.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.Kstoch.MRT <- SS.Kstoch.VRT <- SS.Kstoch.stat <- rep(NA,SS.Kstoch.iter)
for(i in 1:SS.Kstoch.iter) {
  SS.Kstoch.MRT[i] <- SS.Kstoch.RT[[i]]$mRT
  SS.Kstoch.VRT[i] <- SS.Kstoch.RT[[i]]$vRT
  SS.Kstoch.stat[i] <- SS.Kstoch.MRT[i] / SS.Kstoch.VRT[i]
}
hist(SS.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.Kstoch.iter <- dim(PT.Nmat2)[1]
PT.Kstoch.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.Kstoch.MRT <- PT.Kstoch.VRT <- PT.Kstoch.stat <- rep(NA,PT.Kstoch.iter)
for(i in 1:PT.Kstoch.iter) {
  PT.Kstoch.MRT[i] <- PT.Kstoch.RT[[i]]$mRT
  PT.Kstoch.VRT[i] <- PT.Kstoch.RT[[i]]$vRT
  PT.Kstoch.stat[i] <- PT.Kstoch.MRT[i] / PT.Kstoch.VRT[i]
}
hist(PT.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.Kstoch.iter <- dim(SO.Nmat2)[1]
SO.Kstoch.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.Kstoch.MRT <- SO.Kstoch.VRT <- SO.Kstoch.stat <- rep(NA,SO.Kstoch.iter)
for(i in 1:SO.Kstoch.iter) {
  SO.Kstoch.MRT[i] <- SO.Kstoch.RT[[i]]$mRT
  SO.Kstoch.VRT[i] <- SO.Kstoch.RT[[i]]$vRT
  SO.Kstoch.stat[i] <- SO.Kstoch.MRT[i] / SO.Kstoch.VRT[i]
}
hist(SO.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.Kstoch.iter <- dim(MN.Nmat2)[1]
MN.Kstoch.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.Kstoch.MRT <- MN.Kstoch.VRT <- MN.Kstoch.stat <- rep(NA,MN.Kstoch.iter)
for(i in 1:MN.Kstoch.iter) {
  MN.Kstoch.MRT[i] <- MN.Kstoch.RT[[i]]$mRT
  MN.Kstoch.VRT[i] <- MN.Kstoch.RT[[i]]$vRT
  MN.Kstoch.stat[i] <- MN.Kstoch.MRT[i] / MN.Kstoch.VRT[i]
}
hist(MN.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.Kstoch.iter <- dim(OR.Nmat2)[1]
OR.Kstoch.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.Kstoch.MRT <- OR.Kstoch.VRT <- OR.Kstoch.stat <- rep(NA,OR.Kstoch.iter)
for(i in 1:OR.Kstoch.iter) {
  OR.Kstoch.MRT[i] <- OR.Kstoch.RT[[i]]$mRT
  OR.Kstoch.VRT[i] <- OR.Kstoch.RT[[i]]$vRT
  OR.Kstoch.stat[i] <- OR.Kstoch.MRT[i] / OR.Kstoch.VRT[i]
}
hist(OR.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.Kstoch.iter <- dim(NR.Nmat2)[1]
NR.Kstoch.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.Kstoch.MRT <- NR.Kstoch.VRT <- NR.Kstoch.stat <- rep(NA,NR.Kstoch.iter)
for(i in 1:NR.Kstoch.iter) {
  NR.Kstoch.MRT[i] <- NR.Kstoch.RT[[i]]$mRT
  NR.Kstoch.VRT[i] <- NR.Kstoch.RT[[i]]$vRT
  NR.Kstoch.stat[i] <- NR.Kstoch.MRT[i] / NR.Kstoch.VRT[i]
}
hist(NR.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.Kstoch.iter <- dim(GN.Nmat2)[1]
GN.Kstoch.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.Kstoch.MRT <- GN.Kstoch.VRT <- GN.Kstoch.stat <- rep(NA,GN.Kstoch.iter)
for(i in 1:GN.Kstoch.iter) {
  GN.Kstoch.MRT[i] <- GN.Kstoch.RT[[i]]$mRT
  GN.Kstoch.VRT[i] <- GN.Kstoch.RT[[i]]$vRT
  GN.Kstoch.stat[i] <- GN.Kstoch.MRT[i] / GN.Kstoch.VRT[i]
}
hist(GN.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.Kstoch.iter <- dim(DN.Nmat2)[1]
DN.Kstoch.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.Kstoch.MRT <- DN.Kstoch.VRT <- DN.Kstoch.stat <- rep(NA,DN.Kstoch.iter)
for(i in 1:DN.Kstoch.iter) {
  DN.Kstoch.MRT[i] <- DN.Kstoch.RT[[i]]$mRT
  DN.Kstoch.VRT[i] <- DN.Kstoch.RT[[i]]$vRT
  DN.Kstoch.stat[i] <- DN.Kstoch.MRT[i] / DN.Kstoch.VRT[i]
}
hist(DN.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.Kstoch.iter <- dim(AL.Nmat2)[1]
AL.Kstoch.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.Kstoch.MRT <- AL.Kstoch.VRT <- AL.Kstoch.stat <- rep(NA,AL.Kstoch.iter)
for(i in 1:AL.Kstoch.iter) {
  AL.Kstoch.MRT[i] <- AL.Kstoch.RT[[i]]$mRT
  AL.Kstoch.VRT[i] <- AL.Kstoch.RT[[i]]$vRT
  AL.Kstoch.stat[i] <- AL.Kstoch.MRT[i] / AL.Kstoch.VRT[i]
}
hist(AL.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.Kstoch.iter <- dim(TC.Nmat2)[1]
TC.Kstoch.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.Kstoch.MRT <- TC.Kstoch.VRT <- TC.Kstoch.stat <- rep(NA,TC.Kstoch.iter)
for(i in 1:TC.Kstoch.iter) {
  TC.Kstoch.MRT[i] <- TC.Kstoch.RT[[i]]$mRT
  TC.Kstoch.VRT[i] <- TC.Kstoch.RT[[i]]$vRT
  TC.Kstoch.stat[i] <- TC.Kstoch.MRT[i] / TC.Kstoch.VRT[i]
}
hist(TC.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.Kstoch.iter <- dim(TH.Nmat2)[1]
TH.Kstoch.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.Kstoch.MRT <- TH.Kstoch.VRT <- TH.Kstoch.stat <- rep(NA,TH.Kstoch.iter)
for(i in 1:TH.Kstoch.iter) {
  TH.Kstoch.MRT[i] <- TH.Kstoch.RT[[i]]$mRT
  TH.Kstoch.VRT[i] <- TH.Kstoch.RT[[i]]$vRT
  TH.Kstoch.stat[i] <- TH.Kstoch.MRT[i] / TH.Kstoch.VRT[i]
}
hist(TH.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.Kstoch.iter <- dim(SH.Nmat2)[1]
SH.Kstoch.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.Kstoch.MRT <- SH.Kstoch.VRT <- SH.Kstoch.stat <- rep(NA,SH.Kstoch.iter)
for(i in 1:SH.Kstoch.iter) {
  SH.Kstoch.MRT[i] <- SH.Kstoch.RT[[i]]$mRT
  SH.Kstoch.VRT[i] <- SH.Kstoch.RT[[i]]$vRT
  SH.Kstoch.stat[i] <- SH.Kstoch.MRT[i] / SH.Kstoch.VRT[i]
}
hist(SH.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.Kstoch.iter <- dim(DM.Nmat2)[1]
DM.Kstoch.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.Kstoch.MRT <- DM.Kstoch.VRT <- DM.Kstoch.stat <- rep(NA,DM.Kstoch.iter)
for(i in 1:DM.Kstoch.iter) {
  DM.Kstoch.MRT[i] <- DM.Kstoch.RT[[i]]$mRT
  DM.Kstoch.VRT[i] <- DM.Kstoch.RT[[i]]$vRT
  DM.Kstoch.stat[i] <- DM.Kstoch.MRT[i] / DM.Kstoch.VRT[i]
}
hist(DM.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.Kstoch.iter <- dim(TA.Nmat2)[1]
TA.Kstoch.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.Kstoch.MRT <- TA.Kstoch.VRT <- TA.Kstoch.stat <- rep(NA,TA.Kstoch.iter)
for(i in 1:TA.Kstoch.iter) {
  TA.Kstoch.MRT[i] <- TA.Kstoch.RT[[i]]$mRT
  TA.Kstoch.VRT[i] <- TA.Kstoch.RT[[i]]$vRT
  TA.Kstoch.stat[i] <- TA.Kstoch.MRT[i] / TA.Kstoch.VRT[i]
}
hist(TA.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.Kstoch.iter <- dim(MR.Nmat2)[1]
MR.Kstoch.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.Kstoch.MRT <- MR.Kstoch.VRT <- MR.Kstoch.stat <- rep(NA,MR.Kstoch.iter)
for(i in 1:MR.Kstoch.iter) {
  MR.Kstoch.MRT[i] <- MR.Kstoch.RT[[i]]$mRT
  MR.Kstoch.VRT[i] <- MR.Kstoch.RT[[i]]$vRT
  MR.Kstoch.stat[i] <- MR.Kstoch.MRT[i] / MR.Kstoch.VRT[i]
}
hist(MR.Kstoch.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.Kstoch.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.Kstoch.stat.med <- median(DP.Kstoch.stat, na.rm=T)
PA.Kstoch.stat.med <- median(PA.Kstoch.stat, na.rm=T)
ZT.Kstoch.stat.med <- median(ZT.Kstoch.stat, na.rm=T)
PH.Kstoch.stat.med <- median(PH.Kstoch.stat, na.rm=T)
VU.Kstoch.stat.med <- median(VU.Kstoch.stat, na.rm=T)
PG.Kstoch.stat.med <- median(PG.Kstoch.stat, na.rm=T)
SS.Kstoch.stat.med <- median(SS.Kstoch.stat, na.rm=T)
PT.Kstoch.stat.med <- median(PT.Kstoch.stat, na.rm=T)
SO.Kstoch.stat.med <- median(SO.Kstoch.stat, na.rm=T)
MN.Kstoch.stat.med <- median(MN.Kstoch.stat, na.rm=T)
OR.Kstoch.stat.med <- median(OR.Kstoch.stat, na.rm=T)
NR.Kstoch.stat.med <- median(NR.Kstoch.stat, na.rm=T)
GN.Kstoch.stat.med <- median(GN.Kstoch.stat, na.rm=T)
DN.Kstoch.stat.med <- median(DN.Kstoch.stat, na.rm=T)
AL.Kstoch.stat.med <- median(AL.Kstoch.stat, na.rm=T)
TC.Kstoch.stat.med <- median(TC.Kstoch.stat, na.rm=T)
TH.Kstoch.stat.med <- median(TH.Kstoch.stat, na.rm=T)
SH.Kstoch.stat.med <- median(SH.Kstoch.stat, na.rm=T)
DM.Kstoch.stat.med <- median(DM.Kstoch.stat, na.rm=T)
TA.Kstoch.stat.med <- median(TA.Kstoch.stat, na.rm=T)
MR.Kstoch.stat.med <- median(MR.Kstoch.stat, na.rm=T)

DP.Kstoch.stat.up <- quantile(DP.Kstoch.stat, probs=0.975, na.rm=T)
PA.Kstoch.stat.up <- quantile(PA.Kstoch.stat, probs=0.975, na.rm=T)
ZT.Kstoch.stat.up <- quantile(ZT.Kstoch.stat, probs=0.975, na.rm=T)
PH.Kstoch.stat.up <- quantile(PH.Kstoch.stat, probs=0.975, na.rm=T)
VU.Kstoch.stat.up <- quantile(VU.Kstoch.stat, probs=0.975, na.rm=T)
PG.Kstoch.stat.up <- quantile(PG.Kstoch.stat, probs=0.975, na.rm=T)
SS.Kstoch.stat.up <- quantile(SS.Kstoch.stat, probs=0.975, na.rm=T)
PT.Kstoch.stat.up <- quantile(PT.Kstoch.stat, probs=0.975, na.rm=T)
SO.Kstoch.stat.up <- quantile(SO.Kstoch.stat, probs=0.975, na.rm=T)
MN.Kstoch.stat.up <- quantile(MN.Kstoch.stat, probs=0.975, na.rm=T)
OR.Kstoch.stat.up <- quantile(OR.Kstoch.stat, probs=0.975, na.rm=T)
NR.Kstoch.stat.up <- quantile(NR.Kstoch.stat, probs=0.975, na.rm=T)
GN.Kstoch.stat.up <- quantile(GN.Kstoch.stat, probs=0.975, na.rm=T)
DN.Kstoch.stat.up <- quantile(DN.Kstoch.stat, probs=0.975, na.rm=T)
AL.Kstoch.stat.up <- quantile(AL.Kstoch.stat, probs=0.975, na.rm=T)
TC.Kstoch.stat.up <- quantile(TC.Kstoch.stat, probs=0.975, na.rm=T)
TH.Kstoch.stat.up <- quantile(TH.Kstoch.stat, probs=0.975, na.rm=T)
SH.Kstoch.stat.up <- quantile(SH.Kstoch.stat, probs=0.975, na.rm=T)
DM.Kstoch.stat.up <- quantile(DM.Kstoch.stat, probs=0.975, na.rm=T)
TA.Kstoch.stat.up <- quantile(TA.Kstoch.stat, probs=0.975, na.rm=T)
MR.Kstoch.stat.up <- quantile(MR.Kstoch.stat, probs=0.975, na.rm=T)

DP.Kstoch.stat.lo <- quantile(DP.Kstoch.stat, probs=0.025, na.rm=T)
PA.Kstoch.stat.lo <- quantile(PA.Kstoch.stat, probs=0.025, na.rm=T)
ZT.Kstoch.stat.lo <- quantile(ZT.Kstoch.stat, probs=0.025, na.rm=T)
PH.Kstoch.stat.lo <- quantile(PH.Kstoch.stat, probs=0.025, na.rm=T)
VU.Kstoch.stat.lo <- quantile(VU.Kstoch.stat, probs=0.025, na.rm=T)
PG.Kstoch.stat.lo <- quantile(PG.Kstoch.stat, probs=0.025, na.rm=T)
SS.Kstoch.stat.lo <- quantile(SS.Kstoch.stat, probs=0.025, na.rm=T)
PT.Kstoch.stat.lo <- quantile(PT.Kstoch.stat, probs=0.025, na.rm=T)
SO.Kstoch.stat.lo <- quantile(SO.Kstoch.stat, probs=0.025, na.rm=T)
MN.Kstoch.stat.lo <- quantile(MN.Kstoch.stat, probs=0.025, na.rm=T)
OR.Kstoch.stat.lo <- quantile(OR.Kstoch.stat, probs=0.025, na.rm=T)
NR.Kstoch.stat.lo <- quantile(NR.Kstoch.stat, probs=0.025, na.rm=T)
GN.Kstoch.stat.lo <- quantile(GN.Kstoch.stat, probs=0.025, na.rm=T)
DN.Kstoch.stat.lo <- quantile(DN.Kstoch.stat, probs=0.025, na.rm=T)
AL.Kstoch.stat.lo <- quantile(AL.Kstoch.stat, probs=0.025, na.rm=T)
TC.Kstoch.stat.lo <- quantile(TC.Kstoch.stat, probs=0.025, na.rm=T)
TH.Kstoch.stat.lo <- quantile(TH.Kstoch.stat, probs=0.025, na.rm=T)
SH.Kstoch.stat.lo <- quantile(SH.Kstoch.stat, probs=0.025, na.rm=T)
DM.Kstoch.stat.lo <- quantile(DM.Kstoch.stat, probs=0.025, na.rm=T)
TA.Kstoch.stat.lo <- quantile(TA.Kstoch.stat, probs=0.025, na.rm=T)
MR.Kstoch.stat.lo <- quantile(MR.Kstoch.stat, probs=0.025, na.rm=T)

Kstoch.stat.med <- c(DP.Kstoch.stat.med, PA.Kstoch.stat.med, ZT.Kstoch.stat.med, PH.Kstoch.stat.med, VU.Kstoch.stat.med, PG.Kstoch.stat.med,
                     SS.Kstoch.stat.med, PT.Kstoch.stat.med, SO.Kstoch.stat.med, MN.Kstoch.stat.med, OR.Kstoch.stat.med, NR.Kstoch.stat.med,
                     GN.Kstoch.stat.med, DN.Kstoch.stat.med, AL.Kstoch.stat.med, TC.Kstoch.stat.med, TH.Kstoch.stat.med, SH.Kstoch.stat.med,
                     DM.Kstoch.stat.med, TA.Kstoch.stat.med, MR.Kstoch.stat.med)
Kstoch.stat.up <- c(DP.Kstoch.stat.up, PA.Kstoch.stat.up, ZT.Kstoch.stat.up, PH.Kstoch.stat.up, VU.Kstoch.stat.up, PG.Kstoch.stat.up,
                    SS.Kstoch.stat.up, PT.Kstoch.stat.up, SO.Kstoch.stat.up, MN.Kstoch.stat.up, OR.Kstoch.stat.up, NR.Kstoch.stat.up,
                    GN.Kstoch.stat.up, DN.Kstoch.stat.up, AL.Kstoch.stat.up, TC.Kstoch.stat.up, TH.Kstoch.stat.up, SH.Kstoch.stat.up,
                    DM.Kstoch.stat.up, TA.Kstoch.stat.up, MR.Kstoch.stat.up)
Kstoch.stat.lo <- c(DP.Kstoch.stat.lo, PA.Kstoch.stat.lo, ZT.Kstoch.stat.lo, PH.Kstoch.stat.lo, VU.Kstoch.stat.lo, PG.Kstoch.stat.lo,
                    SS.Kstoch.stat.lo, PT.Kstoch.stat.lo, SO.Kstoch.stat.lo, MN.Kstoch.stat.lo, OR.Kstoch.stat.lo, NR.Kstoch.stat.lo,
                    GN.Kstoch.stat.lo, DN.Kstoch.stat.lo, AL.Kstoch.stat.lo, TC.Kstoch.stat.lo, TH.Kstoch.stat.lo, SH.Kstoch.stat.lo,
                    DM.Kstoch.stat.lo, TA.Kstoch.stat.lo, MR.Kstoch.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), Kstoch.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
Kstoch.dat <- data.frame(spp.mass.vec, Kstoch.stat.med, Kstoch.stat.up, Kstoch.stat.lo)
colnames(Kstoch.dat) <- c("M", "statM", "statUP", "statLO")
rownames(Kstoch.dat) <- labs.vec

p <- ggplot(Kstoch.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPKstochstat <- data.frame(rep("DP",length(DP.Kstoch.stat)), DP.Kstoch.stat)
PAKstochstat <- data.frame(rep("PA",length(PA.Kstoch.stat)), PA.Kstoch.stat)
ZTKstochstat <- data.frame(rep("ZT",length(ZT.Kstoch.stat)), ZT.Kstoch.stat)
PHKstochstat <- data.frame(rep("PH",length(PH.Kstoch.stat)), PH.Kstoch.stat)
VUKstochstat <- data.frame(rep("VU",length(VU.Kstoch.stat)), VU.Kstoch.stat)
PGKstochstat <- data.frame(rep("PG",length(PG.Kstoch.stat)), PG.Kstoch.stat)
SSKstochstat <- data.frame(rep("SS",length(SS.Kstoch.stat)), SS.Kstoch.stat)
PTKstochstat <- data.frame(rep("PT",length(PT.Kstoch.stat)), PT.Kstoch.stat)
SOKstochstat <- data.frame(rep("SO",length(SO.Kstoch.stat)), SO.Kstoch.stat)
MNKstochstat <- data.frame(rep("MN",length(MN.Kstoch.stat)), MN.Kstoch.stat)
ORKstochstat <- data.frame(rep("OR",length(OR.Kstoch.stat)), OR.Kstoch.stat)
NRKstochstat <- data.frame(rep("NR",length(NR.Kstoch.stat)), NR.Kstoch.stat)
GNKstochstat <- data.frame(rep("GN",length(GN.Kstoch.stat)), GN.Kstoch.stat)
DNKstochstat <- data.frame(rep("DN",length(DN.Kstoch.stat)), DN.Kstoch.stat)
ALKstochstat <- data.frame(rep("AL",length(AL.Kstoch.stat)), AL.Kstoch.stat)
TCKstochstat <- data.frame(rep("TC",length(TC.Kstoch.stat)), TC.Kstoch.stat)
THKstochstat <- data.frame(rep("TH",length(TH.Kstoch.stat)), TH.Kstoch.stat)
SHKstochstat <- data.frame(rep("SH",length(SH.Kstoch.stat)), SH.Kstoch.stat)
DMKstochstat <- data.frame(rep("DM",length(DM.Kstoch.stat)), DM.Kstoch.stat)
TAKstochstat <- data.frame(rep("TA",length(TA.Kstoch.stat)), TA.Kstoch.stat)
MRKstochstat <- data.frame(rep("MR",length(MR.Kstoch.stat)), MR.Kstoch.stat)

colnames(DPKstochstat) <- c("SP","mRTvRT")
colnames(PAKstochstat) <- c("SP","mRTvRT")
colnames(ZTKstochstat) <- c("SP","mRTvRT")
colnames(PHKstochstat) <- c("SP","mRTvRT")
colnames(VUKstochstat) <- c("SP","mRTvRT")
colnames(PGKstochstat) <- c("SP","mRTvRT")
colnames(SSKstochstat) <- c("SP","mRTvRT")
colnames(PTKstochstat) <- c("SP","mRTvRT")
colnames(SOKstochstat) <- c("SP","mRTvRT")
colnames(MNKstochstat) <- c("SP","mRTvRT")
colnames(ORKstochstat) <- c("SP","mRTvRT")
colnames(NRKstochstat) <- c("SP","mRTvRT")
colnames(GNKstochstat) <- c("SP","mRTvRT")
colnames(DNKstochstat) <- c("SP","mRTvRT")
colnames(ALKstochstat) <- c("SP","mRTvRT")
colnames(TCKstochstat) <- c("SP","mRTvRT")
colnames(THKstochstat) <- c("SP","mRTvRT")
colnames(SHKstochstat) <- c("SP","mRTvRT")
colnames(DMKstochstat) <- c("SP","mRTvRT")
colnames(TAKstochstat) <- c("SP","mRTvRT")
colnames(MRKstochstat) <- c("SP","mRTvRT")

Kstochstat.dat <- rbind(DPKstochstat, PAKstochstat, ZTKstochstat, PHKstochstat, VUKstochstat, PGKstochstat, SSKstochstat, PTKstochstat,
                        SOKstochstat, MNKstochstat, ORKstochstat, NRKstochstat, GNKstochstat, DNKstochstat, ALKstochstat, TCKstochstat,
                        THKstochstat, SHKstochstat, DMKstochstat, TAKstochstat, MRKstochstat)
Kstochstat.dat$SP <- as.factor(Kstochstat.dat$SP)
nspp <- length(table(Kstochstat.dat$SP))

mu <- Kstochstat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(Kstochstat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(Kstochstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(Kstochstat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(Kstochstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("Kstochstatdistrib.RData")

