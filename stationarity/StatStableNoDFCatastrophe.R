#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

######################################################
##  STABLE FROM CATASTROPHE (NO COMPONENT FEEDBACK) ##
######################################################

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
load("DPstableCat.RData")
load("PAstableCat.RData")
load("ZTstableCat.RData")
load("PHstableCat.RData")
load("VUstableCat.RData")
load("PGstableCat.RData")
load("SSstableCat.RData")
load("PTstableCat.RData")
load("SOstableCat.RData")
load("MNstableCat.RData")
load("ORstableCat.RData")
load("NRstableCat.RData")
load("GNstableCat.RData")
load("DNstableCat.RData")
load("ALstableCat.RData")
load("TCstableCat.RData")
load("THstableCat.RData")
load("SHstableCat.RData")
load("DMstableCat.RData")
load("TAstableCat.RData")
load("MRstableCat.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.stableCat.iter <- dim(DP.Nmat2)[1]
DP.stableCat.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.stableCat.MRT <- DP.stableCat.VRT <- DP.stableCat.stat <- rep(NA,DP.stableCat.iter)
for(i in 1:DP.stableCat.iter) {
  DP.stableCat.MRT[i] <- DP.stableCat.RT[[i]]$mRT
  DP.stableCat.VRT[i] <- DP.stableCat.RT[[i]]$vRT
  DP.stableCat.stat[i] <- DP.stableCat.MRT[i] / DP.stableCat.VRT[i]
}
hist(DP.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.stableCat.iter <- dim(PA.Nmat2)[1]
PA.stableCat.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.stableCat.MRT <- PA.stableCat.VRT <- PA.stableCat.stat <- rep(NA,PA.stableCat.iter)
for(i in 1:PA.stableCat.iter) {
  PA.stableCat.MRT[i] <- PA.stableCat.RT[[i]]$mRT
  PA.stableCat.VRT[i] <- PA.stableCat.RT[[i]]$vRT
  PA.stableCat.stat[i] <- PA.stableCat.MRT[i] / PA.stableCat.VRT[i]
}
hist(PA.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.stableCat.iter <- dim(ZT.Nmat2)[1]
ZT.stableCat.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.stableCat.MRT <- ZT.stableCat.VRT <- ZT.stableCat.stat <- rep(NA,ZT.stableCat.iter)
for(i in 1:ZT.stableCat.iter) {
  ZT.stableCat.MRT[i] <- ZT.stableCat.RT[[i]]$mRT
  ZT.stableCat.VRT[i] <- ZT.stableCat.RT[[i]]$vRT
  ZT.stableCat.stat[i] <- ZT.stableCat.MRT[i] / ZT.stableCat.VRT[i]
}
hist(ZT.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.stableCat.iter <- dim(PH.Nmat2)[1]
PH.stableCat.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.stableCat.MRT <- PH.stableCat.VRT <- PH.stableCat.stat <- rep(NA,PH.stableCat.iter)
for(i in 1:PH.stableCat.iter) {
  PH.stableCat.MRT[i] <- PH.stableCat.RT[[i]]$mRT
  PH.stableCat.VRT[i] <- PH.stableCat.RT[[i]]$vRT
  PH.stableCat.stat[i] <- PH.stableCat.MRT[i] / PH.stableCat.VRT[i]
}
hist(PH.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.stableCat.iter <- dim(VU.Nmat2)[1]
VU.stableCat.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.stableCat.MRT <- VU.stableCat.VRT <- VU.stableCat.stat <- rep(NA,VU.stableCat.iter)
for(i in 1:VU.stableCat.iter) {
  VU.stableCat.MRT[i] <- VU.stableCat.RT[[i]]$mRT
  VU.stableCat.VRT[i] <- VU.stableCat.RT[[i]]$vRT
  VU.stableCat.stat[i] <- VU.stableCat.MRT[i] / VU.stableCat.VRT[i]
}
hist(VU.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.stableCat.iter <- dim(PG.Nmat2)[1]
PG.stableCat.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.stableCat.MRT <- PG.stableCat.VRT <- PG.stableCat.stat <- rep(NA,PG.stableCat.iter)
for(i in 1:PG.stableCat.iter) {
  PG.stableCat.MRT[i] <- PG.stableCat.RT[[i]]$mRT
  PG.stableCat.VRT[i] <- PG.stableCat.RT[[i]]$vRT
  PG.stableCat.stat[i] <- PG.stableCat.MRT[i] / PG.stableCat.VRT[i]
}
hist(PG.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.stableCat.iter <- dim(SS.Nmat2)[1]
SS.stableCat.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.stableCat.MRT <- SS.stableCat.VRT <- SS.stableCat.stat <- rep(NA,SS.stableCat.iter)
for(i in 1:SS.stableCat.iter) {
  SS.stableCat.MRT[i] <- SS.stableCat.RT[[i]]$mRT
  SS.stableCat.VRT[i] <- SS.stableCat.RT[[i]]$vRT
  SS.stableCat.stat[i] <- SS.stableCat.MRT[i] / SS.stableCat.VRT[i]
}
hist(SS.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.stableCat.iter <- dim(PT.Nmat2)[1]
PT.stableCat.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.stableCat.MRT <- PT.stableCat.VRT <- PT.stableCat.stat <- rep(NA,PT.stableCat.iter)
for(i in 1:PT.stableCat.iter) {
  PT.stableCat.MRT[i] <- PT.stableCat.RT[[i]]$mRT
  PT.stableCat.VRT[i] <- PT.stableCat.RT[[i]]$vRT
  PT.stableCat.stat[i] <- PT.stableCat.MRT[i] / PT.stableCat.VRT[i]
}
hist(PT.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.stableCat.iter <- dim(SO.Nmat2)[1]
SO.stableCat.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.stableCat.MRT <- SO.stableCat.VRT <- SO.stableCat.stat <- rep(NA,SO.stableCat.iter)
for(i in 1:SO.stableCat.iter) {
  SO.stableCat.MRT[i] <- SO.stableCat.RT[[i]]$mRT
  SO.stableCat.VRT[i] <- SO.stableCat.RT[[i]]$vRT
  SO.stableCat.stat[i] <- SO.stableCat.MRT[i] / SO.stableCat.VRT[i]
}
hist(SO.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.stableCat.iter <- dim(MN.Nmat2)[1]
MN.stableCat.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.stableCat.MRT <- MN.stableCat.VRT <- MN.stableCat.stat <- rep(NA,MN.stableCat.iter)
for(i in 1:MN.stableCat.iter) {
  MN.stableCat.MRT[i] <- MN.stableCat.RT[[i]]$mRT
  MN.stableCat.VRT[i] <- MN.stableCat.RT[[i]]$vRT
  MN.stableCat.stat[i] <- MN.stableCat.MRT[i] / MN.stableCat.VRT[i]
}
hist(MN.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.stableCat.iter <- dim(OR.Nmat2)[1]
OR.stableCat.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.stableCat.MRT <- OR.stableCat.VRT <- OR.stableCat.stat <- rep(NA,OR.stableCat.iter)
for(i in 1:OR.stableCat.iter) {
  OR.stableCat.MRT[i] <- OR.stableCat.RT[[i]]$mRT
  OR.stableCat.VRT[i] <- OR.stableCat.RT[[i]]$vRT
  OR.stableCat.stat[i] <- OR.stableCat.MRT[i] / OR.stableCat.VRT[i]
}
hist(OR.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.stableCat.iter <- dim(NR.Nmat2)[1]
NR.stableCat.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.stableCat.MRT <- NR.stableCat.VRT <- NR.stableCat.stat <- rep(NA,NR.stableCat.iter)
for(i in 1:NR.stableCat.iter) {
  NR.stableCat.MRT[i] <- NR.stableCat.RT[[i]]$mRT
  NR.stableCat.VRT[i] <- NR.stableCat.RT[[i]]$vRT
  NR.stableCat.stat[i] <- NR.stableCat.MRT[i] / NR.stableCat.VRT[i]
}
hist(NR.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.stableCat.iter <- dim(GN.Nmat2)[1]
GN.stableCat.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.stableCat.MRT <- GN.stableCat.VRT <- GN.stableCat.stat <- rep(NA,GN.stableCat.iter)
for(i in 1:GN.stableCat.iter) {
  GN.stableCat.MRT[i] <- GN.stableCat.RT[[i]]$mRT
  GN.stableCat.VRT[i] <- GN.stableCat.RT[[i]]$vRT
  GN.stableCat.stat[i] <- GN.stableCat.MRT[i] / GN.stableCat.VRT[i]
}
hist(GN.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.stableCat.iter <- dim(DN.Nmat2)[1]
DN.stableCat.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.stableCat.MRT <- DN.stableCat.VRT <- DN.stableCat.stat <- rep(NA,DN.stableCat.iter)
for(i in 1:DN.stableCat.iter) {
  DN.stableCat.MRT[i] <- DN.stableCat.RT[[i]]$mRT
  DN.stableCat.VRT[i] <- DN.stableCat.RT[[i]]$vRT
  DN.stableCat.stat[i] <- DN.stableCat.MRT[i] / DN.stableCat.VRT[i]
}
hist(DN.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.stableCat.iter <- dim(AL.Nmat2)[1]
AL.stableCat.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.stableCat.MRT <- AL.stableCat.VRT <- AL.stableCat.stat <- rep(NA,AL.stableCat.iter)
for(i in 1:AL.stableCat.iter) {
  AL.stableCat.MRT[i] <- AL.stableCat.RT[[i]]$mRT
  AL.stableCat.VRT[i] <- AL.stableCat.RT[[i]]$vRT
  AL.stableCat.stat[i] <- AL.stableCat.MRT[i] / AL.stableCat.VRT[i]
}
hist(AL.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.stableCat.iter <- dim(TC.Nmat2)[1]
TC.stableCat.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.stableCat.MRT <- TC.stableCat.VRT <- TC.stableCat.stat <- rep(NA,TC.stableCat.iter)
for(i in 1:TC.stableCat.iter) {
  TC.stableCat.MRT[i] <- TC.stableCat.RT[[i]]$mRT
  TC.stableCat.VRT[i] <- TC.stableCat.RT[[i]]$vRT
  TC.stableCat.stat[i] <- TC.stableCat.MRT[i] / TC.stableCat.VRT[i]
}
hist(TC.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.stableCat.iter <- dim(TH.Nmat2)[1]
TH.stableCat.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.stableCat.MRT <- TH.stableCat.VRT <- TH.stableCat.stat <- rep(NA,TH.stableCat.iter)
for(i in 1:TH.stableCat.iter) {
  TH.stableCat.MRT[i] <- TH.stableCat.RT[[i]]$mRT
  TH.stableCat.VRT[i] <- TH.stableCat.RT[[i]]$vRT
  TH.stableCat.stat[i] <- TH.stableCat.MRT[i] / TH.stableCat.VRT[i]
}
hist(TH.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.stableCat.iter <- dim(SH.Nmat2)[1]
SH.stableCat.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.stableCat.MRT <- SH.stableCat.VRT <- SH.stableCat.stat <- rep(NA,SH.stableCat.iter)
for(i in 1:SH.stableCat.iter) {
  SH.stableCat.MRT[i] <- SH.stableCat.RT[[i]]$mRT
  SH.stableCat.VRT[i] <- SH.stableCat.RT[[i]]$vRT
  SH.stableCat.stat[i] <- SH.stableCat.MRT[i] / SH.stableCat.VRT[i]
}
hist(SH.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.stableCat.iter <- dim(DM.Nmat2)[1]
DM.stableCat.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.stableCat.MRT <- DM.stableCat.VRT <- DM.stableCat.stat <- rep(NA,DM.stableCat.iter)
for(i in 1:DM.stableCat.iter) {
  DM.stableCat.MRT[i] <- DM.stableCat.RT[[i]]$mRT
  DM.stableCat.VRT[i] <- DM.stableCat.RT[[i]]$vRT
  DM.stableCat.stat[i] <- DM.stableCat.MRT[i] / DM.stableCat.VRT[i]
}
hist(DM.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.stableCat.iter <- dim(TA.Nmat2)[1]
TA.stableCat.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.stableCat.MRT <- TA.stableCat.VRT <- TA.stableCat.stat <- rep(NA,TA.stableCat.iter)
for(i in 1:TA.stableCat.iter) {
  TA.stableCat.MRT[i] <- TA.stableCat.RT[[i]]$mRT
  TA.stableCat.VRT[i] <- TA.stableCat.RT[[i]]$vRT
  TA.stableCat.stat[i] <- TA.stableCat.MRT[i] / TA.stableCat.VRT[i]
}
hist(TA.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.stableCat.iter <- dim(MR.Nmat2)[1]
MR.stableCat.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.stableCat.MRT <- MR.stableCat.VRT <- MR.stableCat.stat <- rep(NA,MR.stableCat.iter)
for(i in 1:MR.stableCat.iter) {
  MR.stableCat.MRT[i] <- MR.stableCat.RT[[i]]$mRT
  MR.stableCat.VRT[i] <- MR.stableCat.RT[[i]]$vRT
  MR.stableCat.stat[i] <- MR.stableCat.MRT[i] / MR.stableCat.VRT[i]
}
hist(MR.stableCat.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.stableCat.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


##############
# summaries ##
##############
DP.stableCat.stat.med <- median(DP.stableCat.stat, na.rm=T)
PA.stableCat.stat.med <- median(PA.stableCat.stat, na.rm=T)
ZT.stableCat.stat.med <- median(ZT.stableCat.stat, na.rm=T)
PH.stableCat.stat.med <- median(PH.stableCat.stat, na.rm=T)
VU.stableCat.stat.med <- median(VU.stableCat.stat, na.rm=T)
PG.stableCat.stat.med <- median(PG.stableCat.stat, na.rm=T)
SS.stableCat.stat.med <- median(SS.stableCat.stat, na.rm=T)
PT.stableCat.stat.med <- median(PT.stableCat.stat, na.rm=T)
SO.stableCat.stat.med <- median(SO.stableCat.stat, na.rm=T)
MN.stableCat.stat.med <- median(MN.stableCat.stat, na.rm=T)
OR.stableCat.stat.med <- median(OR.stableCat.stat, na.rm=T)
NR.stableCat.stat.med <- median(NR.stableCat.stat, na.rm=T)
GN.stableCat.stat.med <- median(GN.stableCat.stat, na.rm=T)
DN.stableCat.stat.med <- median(DN.stableCat.stat, na.rm=T)
AL.stableCat.stat.med <- median(AL.stableCat.stat, na.rm=T)
TC.stableCat.stat.med <- median(TC.stableCat.stat, na.rm=T)
TH.stableCat.stat.med <- median(TH.stableCat.stat, na.rm=T)
SH.stableCat.stat.med <- median(SH.stableCat.stat, na.rm=T)
DM.stableCat.stat.med <- median(DM.stableCat.stat, na.rm=T)
TA.stableCat.stat.med <- median(TA.stableCat.stat, na.rm=T)
MR.stableCat.stat.med <- median(MR.stableCat.stat, na.rm=T)

DP.stableCat.stat.up <- quantile(DP.stableCat.stat, probs=0.975, na.rm=T)
PA.stableCat.stat.up <- quantile(PA.stableCat.stat, probs=0.975, na.rm=T)
ZT.stableCat.stat.up <- quantile(ZT.stableCat.stat, probs=0.975, na.rm=T)
PH.stableCat.stat.up <- quantile(PH.stableCat.stat, probs=0.975, na.rm=T)
VU.stableCat.stat.up <- quantile(VU.stableCat.stat, probs=0.975, na.rm=T)
PG.stableCat.stat.up <- quantile(PG.stableCat.stat, probs=0.975, na.rm=T)
SS.stableCat.stat.up <- quantile(SS.stableCat.stat, probs=0.975, na.rm=T)
PT.stableCat.stat.up <- quantile(PT.stableCat.stat, probs=0.975, na.rm=T)
SO.stableCat.stat.up <- quantile(SO.stableCat.stat, probs=0.975, na.rm=T)
MN.stableCat.stat.up <- quantile(MN.stableCat.stat, probs=0.975, na.rm=T)
OR.stableCat.stat.up <- quantile(OR.stableCat.stat, probs=0.975, na.rm=T)
NR.stableCat.stat.up <- quantile(NR.stableCat.stat, probs=0.975, na.rm=T)
GN.stableCat.stat.up <- quantile(GN.stableCat.stat, probs=0.975, na.rm=T)
DN.stableCat.stat.up <- quantile(DN.stableCat.stat, probs=0.975, na.rm=T)
AL.stableCat.stat.up <- quantile(AL.stableCat.stat, probs=0.975, na.rm=T)
TC.stableCat.stat.up <- quantile(TC.stableCat.stat, probs=0.975, na.rm=T)
TH.stableCat.stat.up <- quantile(TH.stableCat.stat, probs=0.975, na.rm=T)
SH.stableCat.stat.up <- quantile(SH.stableCat.stat, probs=0.975, na.rm=T)
DM.stableCat.stat.up <- quantile(DM.stableCat.stat, probs=0.975, na.rm=T)
TA.stableCat.stat.up <- quantile(TA.stableCat.stat, probs=0.975, na.rm=T)
MR.stableCat.stat.up <- quantile(MR.stableCat.stat, probs=0.975, na.rm=T)

DP.stableCat.stat.lo <- quantile(DP.stableCat.stat, probs=0.025, na.rm=T)
PA.stableCat.stat.lo <- quantile(PA.stableCat.stat, probs=0.025, na.rm=T)
ZT.stableCat.stat.lo <- quantile(ZT.stableCat.stat, probs=0.025, na.rm=T)
PH.stableCat.stat.lo <- quantile(PH.stableCat.stat, probs=0.025, na.rm=T)
VU.stableCat.stat.lo <- quantile(VU.stableCat.stat, probs=0.025, na.rm=T)
PG.stableCat.stat.lo <- quantile(PG.stableCat.stat, probs=0.025, na.rm=T)
SS.stableCat.stat.lo <- quantile(SS.stableCat.stat, probs=0.025, na.rm=T)
PT.stableCat.stat.lo <- quantile(PT.stableCat.stat, probs=0.025, na.rm=T)
SO.stableCat.stat.lo <- quantile(SO.stableCat.stat, probs=0.025, na.rm=T)
MN.stableCat.stat.lo <- quantile(MN.stableCat.stat, probs=0.025, na.rm=T)
OR.stableCat.stat.lo <- quantile(OR.stableCat.stat, probs=0.025, na.rm=T)
NR.stableCat.stat.lo <- quantile(NR.stableCat.stat, probs=0.025, na.rm=T)
GN.stableCat.stat.lo <- quantile(GN.stableCat.stat, probs=0.025, na.rm=T)
DN.stableCat.stat.lo <- quantile(DN.stableCat.stat, probs=0.025, na.rm=T)
AL.stableCat.stat.lo <- quantile(AL.stableCat.stat, probs=0.025, na.rm=T)
TC.stableCat.stat.lo <- quantile(TC.stableCat.stat, probs=0.025, na.rm=T)
TH.stableCat.stat.lo <- quantile(TH.stableCat.stat, probs=0.025, na.rm=T)
SH.stableCat.stat.lo <- quantile(SH.stableCat.stat, probs=0.025, na.rm=T)
DM.stableCat.stat.lo <- quantile(DM.stableCat.stat, probs=0.025, na.rm=T)
TA.stableCat.stat.lo <- quantile(TA.stableCat.stat, probs=0.025, na.rm=T)
MR.stableCat.stat.lo <- quantile(MR.stableCat.stat, probs=0.025, na.rm=T)

stableCat.stat.med <- c(DP.stableCat.stat.med, PA.stableCat.stat.med, ZT.stableCat.stat.med, PH.stableCat.stat.med, VU.stableCat.stat.med, PG.stableCat.stat.med,
                        SS.stableCat.stat.med, PT.stableCat.stat.med, SO.stableCat.stat.med, MN.stableCat.stat.med, OR.stableCat.stat.med, NR.stableCat.stat.med,
                        GN.stableCat.stat.med, DN.stableCat.stat.med, AL.stableCat.stat.med, TC.stableCat.stat.med, TH.stableCat.stat.med, SH.stableCat.stat.med,
                        DM.stableCat.stat.med, TA.stableCat.stat.med, MR.stableCat.stat.med)
stableCat.stat.up <- c(DP.stableCat.stat.up, PA.stableCat.stat.up, ZT.stableCat.stat.up, PH.stableCat.stat.up, VU.stableCat.stat.up, PG.stableCat.stat.up,
                       SS.stableCat.stat.up, PT.stableCat.stat.up, SO.stableCat.stat.up, MN.stableCat.stat.up, OR.stableCat.stat.up, NR.stableCat.stat.up,
                       GN.stableCat.stat.up, DN.stableCat.stat.up, AL.stableCat.stat.up, TC.stableCat.stat.up, TH.stableCat.stat.up, SH.stableCat.stat.up,
                       DM.stableCat.stat.up, TA.stableCat.stat.up, MR.stableCat.stat.up)
stableCat.stat.lo <- c(DP.stableCat.stat.lo, PA.stableCat.stat.lo, ZT.stableCat.stat.lo, PH.stableCat.stat.lo, VU.stableCat.stat.lo, PG.stableCat.stat.lo,
                       SS.stableCat.stat.lo, PT.stableCat.stat.lo, SO.stableCat.stat.lo, MN.stableCat.stat.lo, OR.stableCat.stat.lo, NR.stableCat.stat.lo,
                       GN.stableCat.stat.lo, DN.stableCat.stat.lo, AL.stableCat.stat.lo, TC.stableCat.stat.lo, TH.stableCat.stat.lo, SH.stableCat.stat.lo,
                       DM.stableCat.stat.lo, TA.stableCat.stat.lo, MR.stableCat.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), stableCat.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
stableCat.dat <- data.frame(spp.mass.vec, stableCat.stat.med, stableCat.stat.up, stableCat.stat.lo)
colnames(stableCat.dat) <- c("M", "statM", "statUP", "statLO")
rownames(stableCat.dat) <- labs.vec

p <- ggplot(stableCat.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPstableCatstat <- data.frame(rep("DP",length(DP.stableCat.stat)), DP.stableCat.stat)
PAstableCatstat <- data.frame(rep("PA",length(PA.stableCat.stat)), PA.stableCat.stat)
ZTstableCatstat <- data.frame(rep("ZT",length(ZT.stableCat.stat)), ZT.stableCat.stat)
PHstableCatstat <- data.frame(rep("PH",length(PH.stableCat.stat)), PH.stableCat.stat)
VUstableCatstat <- data.frame(rep("VU",length(VU.stableCat.stat)), VU.stableCat.stat)
PGstableCatstat <- data.frame(rep("PG",length(PG.stableCat.stat)), PG.stableCat.stat)
SSstableCatstat <- data.frame(rep("SS",length(SS.stableCat.stat)), SS.stableCat.stat)
PTstableCatstat <- data.frame(rep("PT",length(PT.stableCat.stat)), PT.stableCat.stat)
SOstableCatstat <- data.frame(rep("SO",length(SO.stableCat.stat)), SO.stableCat.stat)
MNstableCatstat <- data.frame(rep("MN",length(MN.stableCat.stat)), MN.stableCat.stat)
ORstableCatstat <- data.frame(rep("OR",length(OR.stableCat.stat)), OR.stableCat.stat)
NRstableCatstat <- data.frame(rep("NR",length(NR.stableCat.stat)), NR.stableCat.stat)
GNstableCatstat <- data.frame(rep("GN",length(GN.stableCat.stat)), GN.stableCat.stat)
DNstableCatstat <- data.frame(rep("DN",length(DN.stableCat.stat)), DN.stableCat.stat)
ALstableCatstat <- data.frame(rep("AL",length(AL.stableCat.stat)), AL.stableCat.stat)
TCstableCatstat <- data.frame(rep("TC",length(TC.stableCat.stat)), TC.stableCat.stat)
THstableCatstat <- data.frame(rep("TH",length(TH.stableCat.stat)), TH.stableCat.stat)
SHstableCatstat <- data.frame(rep("SH",length(SH.stableCat.stat)), SH.stableCat.stat)
DMstableCatstat <- data.frame(rep("DM",length(DM.stableCat.stat)), DM.stableCat.stat)
TAstableCatstat <- data.frame(rep("TA",length(TA.stableCat.stat)), TA.stableCat.stat)
MRstableCatstat <- data.frame(rep("MR",length(MR.stableCat.stat)), MR.stableCat.stat)

colnames(DPstableCatstat) <- c("SP","mRTvRT")
colnames(PAstableCatstat) <- c("SP","mRTvRT")
colnames(ZTstableCatstat) <- c("SP","mRTvRT")
colnames(PHstableCatstat) <- c("SP","mRTvRT")
colnames(VUstableCatstat) <- c("SP","mRTvRT")
colnames(PGstableCatstat) <- c("SP","mRTvRT")
colnames(SSstableCatstat) <- c("SP","mRTvRT")
colnames(PTstableCatstat) <- c("SP","mRTvRT")
colnames(SOstableCatstat) <- c("SP","mRTvRT")
colnames(MNstableCatstat) <- c("SP","mRTvRT")
colnames(ORstableCatstat) <- c("SP","mRTvRT")
colnames(NRstableCatstat) <- c("SP","mRTvRT")
colnames(GNstableCatstat) <- c("SP","mRTvRT")
colnames(DNstableCatstat) <- c("SP","mRTvRT")
colnames(ALstableCatstat) <- c("SP","mRTvRT")
colnames(TCstableCatstat) <- c("SP","mRTvRT")
colnames(THstableCatstat) <- c("SP","mRTvRT")
colnames(SHstableCatstat) <- c("SP","mRTvRT")
colnames(DMstableCatstat) <- c("SP","mRTvRT")
colnames(TAstableCatstat) <- c("SP","mRTvRT")
colnames(MRstableCatstat) <- c("SP","mRTvRT")

stableCatstat.dat <- rbind(DPstableCatstat, PAstableCatstat, ZTstableCatstat, PHstableCatstat, VUstableCatstat, PGstableCatstat, SSstableCatstat, PTstableCatstat,
                           SOstableCatstat, MNstableCatstat, ORstableCatstat, NRstableCatstat, GNstableCatstat, DNstableCatstat, ALstableCatstat, TCstableCatstat,
                           THstableCatstat, SHstableCatstat, DMstableCatstat, TAstableCatstat, MRstableCatstat)
stableCatstat.dat$SP <- as.factor(stableCatstat.dat$SP)
nspp <- length(table(stableCatstat.dat$SP))

mu <- stableCatstat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(stableCatstat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(stableCatstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(stableCatstat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(stableCatstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("stableCatstatdistrib.RData")

