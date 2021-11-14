#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

##########################################################
## 40 G (STABLE, with catastrophe)
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
load("DPstable40G.RData")
load("PAstable40G.RData")
load("ZTstable40G.RData")
load("PHstable40G.RData")
load("VUstable40G.RData")
load("PGstable40G.RData")
load("SSstable40G.RData")
load("PTstable40G.RData")
load("SOstable40G.RData")
load("MNstable40G.RData")
load("ORstable40G.RData")
load("NRstable40G.RData")
load("GNstable40G.RData")
load("DNstable40G.RData")
load("ALstable40G.RData")
load("TCstable40G.RData")
load("THstable40G.RData")
load("SHstable40G.RData")
load("DMstable40G.RData")
load("TAstable40G.RData")
load("MRstable40G.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.stable40G.iter <- dim(DP.Nmat2)[1]
DP.stable40G.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.stable40G.MRT <- DP.stable40G.VRT <- DP.stable40G.stat <- rep(NA,DP.stable40G.iter)
for(i in 1:DP.stable40G.iter) {
  DP.stable40G.MRT[i] <- DP.stable40G.RT[[i]]$mRT
  DP.stable40G.VRT[i] <- DP.stable40G.RT[[i]]$vRT
  DP.stable40G.stat[i] <- DP.stable40G.MRT[i] / DP.stable40G.VRT[i]
}
hist(DP.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.stable40G.iter <- dim(PA.Nmat2)[1]
PA.stable40G.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.stable40G.MRT <- PA.stable40G.VRT <- PA.stable40G.stat <- rep(NA,PA.stable40G.iter)
for(i in 1:PA.stable40G.iter) {
  PA.stable40G.MRT[i] <- PA.stable40G.RT[[i]]$mRT
  PA.stable40G.VRT[i] <- PA.stable40G.RT[[i]]$vRT
  PA.stable40G.stat[i] <- PA.stable40G.MRT[i] / PA.stable40G.VRT[i]
}
hist(PA.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.stable40G.iter <- dim(ZT.Nmat2)[1]
ZT.stable40G.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.stable40G.MRT <- ZT.stable40G.VRT <- ZT.stable40G.stat <- rep(NA,ZT.stable40G.iter)
for(i in 1:ZT.stable40G.iter) {
  ZT.stable40G.MRT[i] <- ZT.stable40G.RT[[i]]$mRT
  ZT.stable40G.VRT[i] <- ZT.stable40G.RT[[i]]$vRT
  ZT.stable40G.stat[i] <- ZT.stable40G.MRT[i] / ZT.stable40G.VRT[i]
}
hist(ZT.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.stable40G.iter <- dim(PH.Nmat2)[1]
PH.stable40G.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.stable40G.MRT <- PH.stable40G.VRT <- PH.stable40G.stat <- rep(NA,PH.stable40G.iter)
for(i in 1:PH.stable40G.iter) {
  PH.stable40G.MRT[i] <- PH.stable40G.RT[[i]]$mRT
  PH.stable40G.VRT[i] <- PH.stable40G.RT[[i]]$vRT
  PH.stable40G.stat[i] <- PH.stable40G.MRT[i] / PH.stable40G.VRT[i]
}
hist(PH.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.stable40G.iter <- dim(VU.Nmat2)[1]
VU.stable40G.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.stable40G.MRT <- VU.stable40G.VRT <- VU.stable40G.stat <- rep(NA,VU.stable40G.iter)
for(i in 1:VU.stable40G.iter) {
  VU.stable40G.MRT[i] <- VU.stable40G.RT[[i]]$mRT
  VU.stable40G.VRT[i] <- VU.stable40G.RT[[i]]$vRT
  VU.stable40G.stat[i] <- VU.stable40G.MRT[i] / VU.stable40G.VRT[i]
}
hist(VU.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.stable40G.iter <- dim(PG.Nmat2)[1]
PG.stable40G.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.stable40G.MRT <- PG.stable40G.VRT <- PG.stable40G.stat <- rep(NA,PG.stable40G.iter)
for(i in 1:PG.stable40G.iter) {
  PG.stable40G.MRT[i] <- PG.stable40G.RT[[i]]$mRT
  PG.stable40G.VRT[i] <- PG.stable40G.RT[[i]]$vRT
  PG.stable40G.stat[i] <- PG.stable40G.MRT[i] / PG.stable40G.VRT[i]
}
hist(PG.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.stable40G.iter <- dim(SS.Nmat2)[1]
SS.stable40G.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.stable40G.MRT <- SS.stable40G.VRT <- SS.stable40G.stat <- rep(NA,SS.stable40G.iter)
for(i in 1:SS.stable40G.iter) {
  SS.stable40G.MRT[i] <- SS.stable40G.RT[[i]]$mRT
  SS.stable40G.VRT[i] <- SS.stable40G.RT[[i]]$vRT
  SS.stable40G.stat[i] <- SS.stable40G.MRT[i] / SS.stable40G.VRT[i]
}
hist(SS.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.stable40G.iter <- dim(PT.Nmat2)[1]
PT.stable40G.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.stable40G.MRT <- PT.stable40G.VRT <- PT.stable40G.stat <- rep(NA,PT.stable40G.iter)
for(i in 1:PT.stable40G.iter) {
  PT.stable40G.MRT[i] <- PT.stable40G.RT[[i]]$mRT
  PT.stable40G.VRT[i] <- PT.stable40G.RT[[i]]$vRT
  PT.stable40G.stat[i] <- PT.stable40G.MRT[i] / PT.stable40G.VRT[i]
}
hist(PT.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.stable40G.iter <- dim(SO.Nmat2)[1]
SO.stable40G.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.stable40G.MRT <- SO.stable40G.VRT <- SO.stable40G.stat <- rep(NA,SO.stable40G.iter)
for(i in 1:SO.stable40G.iter) {
  SO.stable40G.MRT[i] <- SO.stable40G.RT[[i]]$mRT
  SO.stable40G.VRT[i] <- SO.stable40G.RT[[i]]$vRT
  SO.stable40G.stat[i] <- SO.stable40G.MRT[i] / SO.stable40G.VRT[i]
}
hist(SO.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.stable40G.iter <- dim(MN.Nmat2)[1]
MN.stable40G.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.stable40G.MRT <- MN.stable40G.VRT <- MN.stable40G.stat <- rep(NA,MN.stable40G.iter)
for(i in 1:MN.stable40G.iter) {
  MN.stable40G.MRT[i] <- MN.stable40G.RT[[i]]$mRT
  MN.stable40G.VRT[i] <- MN.stable40G.RT[[i]]$vRT
  MN.stable40G.stat[i] <- MN.stable40G.MRT[i] / MN.stable40G.VRT[i]
}
hist(MN.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.stable40G.iter <- dim(OR.Nmat2)[1]
OR.stable40G.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.stable40G.MRT <- OR.stable40G.VRT <- OR.stable40G.stat <- rep(NA,OR.stable40G.iter)
for(i in 1:OR.stable40G.iter) {
  OR.stable40G.MRT[i] <- OR.stable40G.RT[[i]]$mRT
  OR.stable40G.VRT[i] <- OR.stable40G.RT[[i]]$vRT
  OR.stable40G.stat[i] <- OR.stable40G.MRT[i] / OR.stable40G.VRT[i]
}
hist(OR.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.stable40G.iter <- dim(NR.Nmat2)[1]
NR.stable40G.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.stable40G.MRT <- NR.stable40G.VRT <- NR.stable40G.stat <- rep(NA,NR.stable40G.iter)
for(i in 1:NR.stable40G.iter) {
  NR.stable40G.MRT[i] <- NR.stable40G.RT[[i]]$mRT
  NR.stable40G.VRT[i] <- NR.stable40G.RT[[i]]$vRT
  NR.stable40G.stat[i] <- NR.stable40G.MRT[i] / NR.stable40G.VRT[i]
}
hist(NR.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.stable40G.iter <- dim(GN.Nmat2)[1]
GN.stable40G.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.stable40G.MRT <- GN.stable40G.VRT <- GN.stable40G.stat <- rep(NA,GN.stable40G.iter)
for(i in 1:GN.stable40G.iter) {
  GN.stable40G.MRT[i] <- GN.stable40G.RT[[i]]$mRT
  GN.stable40G.VRT[i] <- GN.stable40G.RT[[i]]$vRT
  GN.stable40G.stat[i] <- GN.stable40G.MRT[i] / GN.stable40G.VRT[i]
}
hist(GN.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.stable40G.iter <- dim(DN.Nmat2)[1]
DN.stable40G.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.stable40G.MRT <- DN.stable40G.VRT <- DN.stable40G.stat <- rep(NA,DN.stable40G.iter)
for(i in 1:DN.stable40G.iter) {
  DN.stable40G.MRT[i] <- DN.stable40G.RT[[i]]$mRT
  DN.stable40G.VRT[i] <- DN.stable40G.RT[[i]]$vRT
  DN.stable40G.stat[i] <- DN.stable40G.MRT[i] / DN.stable40G.VRT[i]
}
hist(DN.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.stable40G.iter <- dim(AL.Nmat2)[1]
AL.stable40G.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.stable40G.MRT <- AL.stable40G.VRT <- AL.stable40G.stat <- rep(NA,AL.stable40G.iter)
for(i in 1:AL.stable40G.iter) {
  AL.stable40G.MRT[i] <- AL.stable40G.RT[[i]]$mRT
  AL.stable40G.VRT[i] <- AL.stable40G.RT[[i]]$vRT
  AL.stable40G.stat[i] <- AL.stable40G.MRT[i] / AL.stable40G.VRT[i]
}
hist(AL.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.stable40G.iter <- dim(TC.Nmat2)[1]
TC.stable40G.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.stable40G.MRT <- TC.stable40G.VRT <- TC.stable40G.stat <- rep(NA,TC.stable40G.iter)
for(i in 1:TC.stable40G.iter) {
  TC.stable40G.MRT[i] <- TC.stable40G.RT[[i]]$mRT
  TC.stable40G.VRT[i] <- TC.stable40G.RT[[i]]$vRT
  TC.stable40G.stat[i] <- TC.stable40G.MRT[i] / TC.stable40G.VRT[i]
}
hist(TC.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.stable40G.iter <- dim(TH.Nmat2)[1]
TH.stable40G.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.stable40G.MRT <- TH.stable40G.VRT <- TH.stable40G.stat <- rep(NA,TH.stable40G.iter)
for(i in 1:TH.stable40G.iter) {
  TH.stable40G.MRT[i] <- TH.stable40G.RT[[i]]$mRT
  TH.stable40G.VRT[i] <- TH.stable40G.RT[[i]]$vRT
  TH.stable40G.stat[i] <- TH.stable40G.MRT[i] / TH.stable40G.VRT[i]
}
hist(TH.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.stable40G.iter <- dim(SH.Nmat2)[1]
SH.stable40G.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.stable40G.MRT <- SH.stable40G.VRT <- SH.stable40G.stat <- rep(NA,SH.stable40G.iter)
for(i in 1:SH.stable40G.iter) {
  SH.stable40G.MRT[i] <- SH.stable40G.RT[[i]]$mRT
  SH.stable40G.VRT[i] <- SH.stable40G.RT[[i]]$vRT
  SH.stable40G.stat[i] <- SH.stable40G.MRT[i] / SH.stable40G.VRT[i]
}
hist(SH.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.stable40G.iter <- dim(DM.Nmat2)[1]
DM.stable40G.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.stable40G.MRT <- DM.stable40G.VRT <- DM.stable40G.stat <- rep(NA,DM.stable40G.iter)
for(i in 1:DM.stable40G.iter) {
  DM.stable40G.MRT[i] <- DM.stable40G.RT[[i]]$mRT
  DM.stable40G.VRT[i] <- DM.stable40G.RT[[i]]$vRT
  DM.stable40G.stat[i] <- DM.stable40G.MRT[i] / DM.stable40G.VRT[i]
}
hist(DM.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.stable40G.iter <- dim(TA.Nmat2)[1]
TA.stable40G.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.stable40G.MRT <- TA.stable40G.VRT <- TA.stable40G.stat <- rep(NA,TA.stable40G.iter)
for(i in 1:TA.stable40G.iter) {
  TA.stable40G.MRT[i] <- TA.stable40G.RT[[i]]$mRT
  TA.stable40G.VRT[i] <- TA.stable40G.RT[[i]]$vRT
  TA.stable40G.stat[i] <- TA.stable40G.MRT[i] / TA.stable40G.VRT[i]
}
hist(TA.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.stable40G.iter <- dim(MR.Nmat2)[1]
MR.stable40G.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.stable40G.MRT <- MR.stable40G.VRT <- MR.stable40G.stat <- rep(NA,MR.stable40G.iter)
for(i in 1:MR.stable40G.iter) {
  MR.stable40G.MRT[i] <- MR.stable40G.RT[[i]]$mRT
  MR.stable40G.VRT[i] <- MR.stable40G.RT[[i]]$vRT
  MR.stable40G.stat[i] <- MR.stable40G.MRT[i] / MR.stable40G.VRT[i]
}
hist(MR.stable40G.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.stable40G.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.stable40G.stat.med <- median(DP.stable40G.stat, na.rm=T)
PA.stable40G.stat.med <- median(PA.stable40G.stat, na.rm=T)
ZT.stable40G.stat.med <- median(ZT.stable40G.stat, na.rm=T)
PH.stable40G.stat.med <- median(PH.stable40G.stat, na.rm=T)
VU.stable40G.stat.med <- median(VU.stable40G.stat, na.rm=T)
PG.stable40G.stat.med <- median(PG.stable40G.stat, na.rm=T)
SS.stable40G.stat.med <- median(SS.stable40G.stat, na.rm=T)
PT.stable40G.stat.med <- median(PT.stable40G.stat, na.rm=T)
SO.stable40G.stat.med <- median(SO.stable40G.stat, na.rm=T)
MN.stable40G.stat.med <- median(MN.stable40G.stat, na.rm=T)
OR.stable40G.stat.med <- median(OR.stable40G.stat, na.rm=T)
NR.stable40G.stat.med <- median(NR.stable40G.stat, na.rm=T)
GN.stable40G.stat.med <- median(GN.stable40G.stat, na.rm=T)
DN.stable40G.stat.med <- median(DN.stable40G.stat, na.rm=T)
AL.stable40G.stat.med <- median(AL.stable40G.stat, na.rm=T)
TC.stable40G.stat.med <- median(TC.stable40G.stat, na.rm=T)
TH.stable40G.stat.med <- median(TH.stable40G.stat, na.rm=T)
SH.stable40G.stat.med <- median(SH.stable40G.stat, na.rm=T)
DM.stable40G.stat.med <- median(DM.stable40G.stat, na.rm=T)
TA.stable40G.stat.med <- median(TA.stable40G.stat, na.rm=T)
MR.stable40G.stat.med <- median(MR.stable40G.stat, na.rm=T)

DP.stable40G.stat.up <- quantile(DP.stable40G.stat, probs=0.975, na.rm=T)
PA.stable40G.stat.up <- quantile(PA.stable40G.stat, probs=0.975, na.rm=T)
ZT.stable40G.stat.up <- quantile(ZT.stable40G.stat, probs=0.975, na.rm=T)
PH.stable40G.stat.up <- quantile(PH.stable40G.stat, probs=0.975, na.rm=T)
VU.stable40G.stat.up <- quantile(VU.stable40G.stat, probs=0.975, na.rm=T)
PG.stable40G.stat.up <- quantile(PG.stable40G.stat, probs=0.975, na.rm=T)
SS.stable40G.stat.up <- quantile(SS.stable40G.stat, probs=0.975, na.rm=T)
PT.stable40G.stat.up <- quantile(PT.stable40G.stat, probs=0.975, na.rm=T)
SO.stable40G.stat.up <- quantile(SO.stable40G.stat, probs=0.975, na.rm=T)
MN.stable40G.stat.up <- quantile(MN.stable40G.stat, probs=0.975, na.rm=T)
OR.stable40G.stat.up <- quantile(OR.stable40G.stat, probs=0.975, na.rm=T)
NR.stable40G.stat.up <- quantile(NR.stable40G.stat, probs=0.975, na.rm=T)
GN.stable40G.stat.up <- quantile(GN.stable40G.stat, probs=0.975, na.rm=T)
DN.stable40G.stat.up <- quantile(DN.stable40G.stat, probs=0.975, na.rm=T)
AL.stable40G.stat.up <- quantile(AL.stable40G.stat, probs=0.975, na.rm=T)
TC.stable40G.stat.up <- quantile(TC.stable40G.stat, probs=0.975, na.rm=T)
TH.stable40G.stat.up <- quantile(TH.stable40G.stat, probs=0.975, na.rm=T)
SH.stable40G.stat.up <- quantile(SH.stable40G.stat, probs=0.975, na.rm=T)
DM.stable40G.stat.up <- quantile(DM.stable40G.stat, probs=0.975, na.rm=T)
TA.stable40G.stat.up <- quantile(TA.stable40G.stat, probs=0.975, na.rm=T)
MR.stable40G.stat.up <- quantile(MR.stable40G.stat, probs=0.975, na.rm=T)

DP.stable40G.stat.lo <- quantile(DP.stable40G.stat, probs=0.025, na.rm=T)
PA.stable40G.stat.lo <- quantile(PA.stable40G.stat, probs=0.025, na.rm=T)
ZT.stable40G.stat.lo <- quantile(ZT.stable40G.stat, probs=0.025, na.rm=T)
PH.stable40G.stat.lo <- quantile(PH.stable40G.stat, probs=0.025, na.rm=T)
VU.stable40G.stat.lo <- quantile(VU.stable40G.stat, probs=0.025, na.rm=T)
PG.stable40G.stat.lo <- quantile(PG.stable40G.stat, probs=0.025, na.rm=T)
SS.stable40G.stat.lo <- quantile(SS.stable40G.stat, probs=0.025, na.rm=T)
PT.stable40G.stat.lo <- quantile(PT.stable40G.stat, probs=0.025, na.rm=T)
SO.stable40G.stat.lo <- quantile(SO.stable40G.stat, probs=0.025, na.rm=T)
MN.stable40G.stat.lo <- quantile(MN.stable40G.stat, probs=0.025, na.rm=T)
OR.stable40G.stat.lo <- quantile(OR.stable40G.stat, probs=0.025, na.rm=T)
NR.stable40G.stat.lo <- quantile(NR.stable40G.stat, probs=0.025, na.rm=T)
GN.stable40G.stat.lo <- quantile(GN.stable40G.stat, probs=0.025, na.rm=T)
DN.stable40G.stat.lo <- quantile(DN.stable40G.stat, probs=0.025, na.rm=T)
AL.stable40G.stat.lo <- quantile(AL.stable40G.stat, probs=0.025, na.rm=T)
TC.stable40G.stat.lo <- quantile(TC.stable40G.stat, probs=0.025, na.rm=T)
TH.stable40G.stat.lo <- quantile(TH.stable40G.stat, probs=0.025, na.rm=T)
SH.stable40G.stat.lo <- quantile(SH.stable40G.stat, probs=0.025, na.rm=T)
DM.stable40G.stat.lo <- quantile(DM.stable40G.stat, probs=0.025, na.rm=T)
TA.stable40G.stat.lo <- quantile(TA.stable40G.stat, probs=0.025, na.rm=T)
MR.stable40G.stat.lo <- quantile(MR.stable40G.stat, probs=0.025, na.rm=T)

stable40G.stat.med <- c(DP.stable40G.stat.med, PA.stable40G.stat.med, ZT.stable40G.stat.med, PH.stable40G.stat.med, VU.stable40G.stat.med, PG.stable40G.stat.med,
                        SS.stable40G.stat.med, PT.stable40G.stat.med, SO.stable40G.stat.med, MN.stable40G.stat.med, OR.stable40G.stat.med, NR.stable40G.stat.med,
                        GN.stable40G.stat.med, DN.stable40G.stat.med, AL.stable40G.stat.med, TC.stable40G.stat.med, TH.stable40G.stat.med, SH.stable40G.stat.med,
                        DM.stable40G.stat.med, TA.stable40G.stat.med, MR.stable40G.stat.med)
stable40G.stat.up <- c(DP.stable40G.stat.up, PA.stable40G.stat.up, ZT.stable40G.stat.up, PH.stable40G.stat.up, VU.stable40G.stat.up, PG.stable40G.stat.up,
                       SS.stable40G.stat.up, PT.stable40G.stat.up, SO.stable40G.stat.up, MN.stable40G.stat.up, OR.stable40G.stat.up, NR.stable40G.stat.up,
                       GN.stable40G.stat.up, DN.stable40G.stat.up, AL.stable40G.stat.up, TC.stable40G.stat.up, TH.stable40G.stat.up, SH.stable40G.stat.up,
                       DM.stable40G.stat.up, TA.stable40G.stat.up, MR.stable40G.stat.up)
stable40G.stat.lo <- c(DP.stable40G.stat.lo, PA.stable40G.stat.lo, ZT.stable40G.stat.lo, PH.stable40G.stat.lo, VU.stable40G.stat.lo, PG.stable40G.stat.lo,
                       SS.stable40G.stat.lo, PT.stable40G.stat.lo, SO.stable40G.stat.lo, MN.stable40G.stat.lo, OR.stable40G.stat.lo, NR.stable40G.stat.lo,
                       GN.stable40G.stat.lo, DN.stable40G.stat.lo, AL.stable40G.stat.lo, TC.stable40G.stat.lo, TH.stable40G.stat.lo, SH.stable40G.stat.lo,
                       DM.stable40G.stat.lo, TA.stable40G.stat.lo, MR.stable40G.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), stable40G.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
stable40G.dat <- data.frame(spp.mass.vec, stable40G.stat.med, stable40G.stat.up, stable40G.stat.lo)
colnames(stable40G.dat) <- c("M", "statM", "statUP", "statLO")
rownames(stable40G.dat) <- labs.vec

p <- ggplot(stable40G.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPstable40Gstat <- data.frame(rep("DP",length(DP.stable40G.stat)), DP.stable40G.stat)
PAstable40Gstat <- data.frame(rep("PA",length(PA.stable40G.stat)), PA.stable40G.stat)
ZTstable40Gstat <- data.frame(rep("ZT",length(ZT.stable40G.stat)), ZT.stable40G.stat)
PHstable40Gstat <- data.frame(rep("PH",length(PH.stable40G.stat)), PH.stable40G.stat)
VUstable40Gstat <- data.frame(rep("VU",length(VU.stable40G.stat)), VU.stable40G.stat)
PGstable40Gstat <- data.frame(rep("PG",length(PG.stable40G.stat)), PG.stable40G.stat)
SSstable40Gstat <- data.frame(rep("SS",length(SS.stable40G.stat)), SS.stable40G.stat)
PTstable40Gstat <- data.frame(rep("PT",length(PT.stable40G.stat)), PT.stable40G.stat)
SOstable40Gstat <- data.frame(rep("SO",length(SO.stable40G.stat)), SO.stable40G.stat)
MNstable40Gstat <- data.frame(rep("MN",length(MN.stable40G.stat)), MN.stable40G.stat)
ORstable40Gstat <- data.frame(rep("OR",length(OR.stable40G.stat)), OR.stable40G.stat)
NRstable40Gstat <- data.frame(rep("NR",length(NR.stable40G.stat)), NR.stable40G.stat)
GNstable40Gstat <- data.frame(rep("GN",length(GN.stable40G.stat)), GN.stable40G.stat)
DNstable40Gstat <- data.frame(rep("DN",length(DN.stable40G.stat)), DN.stable40G.stat)
ALstable40Gstat <- data.frame(rep("AL",length(AL.stable40G.stat)), AL.stable40G.stat)
TCstable40Gstat <- data.frame(rep("TC",length(TC.stable40G.stat)), TC.stable40G.stat)
THstable40Gstat <- data.frame(rep("TH",length(TH.stable40G.stat)), TH.stable40G.stat)
SHstable40Gstat <- data.frame(rep("SH",length(SH.stable40G.stat)), SH.stable40G.stat)
DMstable40Gstat <- data.frame(rep("DM",length(DM.stable40G.stat)), DM.stable40G.stat)
TAstable40Gstat <- data.frame(rep("TA",length(TA.stable40G.stat)), TA.stable40G.stat)
MRstable40Gstat <- data.frame(rep("MR",length(MR.stable40G.stat)), MR.stable40G.stat)

colnames(DPstable40Gstat) <- c("SP","mRTvRT")
colnames(PAstable40Gstat) <- c("SP","mRTvRT")
colnames(ZTstable40Gstat) <- c("SP","mRTvRT")
colnames(PHstable40Gstat) <- c("SP","mRTvRT")
colnames(VUstable40Gstat) <- c("SP","mRTvRT")
colnames(PGstable40Gstat) <- c("SP","mRTvRT")
colnames(SSstable40Gstat) <- c("SP","mRTvRT")
colnames(PTstable40Gstat) <- c("SP","mRTvRT")
colnames(SOstable40Gstat) <- c("SP","mRTvRT")
colnames(MNstable40Gstat) <- c("SP","mRTvRT")
colnames(ORstable40Gstat) <- c("SP","mRTvRT")
colnames(NRstable40Gstat) <- c("SP","mRTvRT")
colnames(GNstable40Gstat) <- c("SP","mRTvRT")
colnames(DNstable40Gstat) <- c("SP","mRTvRT")
colnames(ALstable40Gstat) <- c("SP","mRTvRT")
colnames(TCstable40Gstat) <- c("SP","mRTvRT")
colnames(THstable40Gstat) <- c("SP","mRTvRT")
colnames(SHstable40Gstat) <- c("SP","mRTvRT")
colnames(DMstable40Gstat) <- c("SP","mRTvRT")
colnames(TAstable40Gstat) <- c("SP","mRTvRT")
colnames(MRstable40Gstat) <- c("SP","mRTvRT")

stable40Gstat.dat <- rbind(DPstable40Gstat, PAstable40Gstat, ZTstable40Gstat, PHstable40Gstat, VUstable40Gstat, PGstable40Gstat, SSstable40Gstat, PTstable40Gstat,
                           SOstable40Gstat, MNstable40Gstat, ORstable40Gstat, NRstable40Gstat, GNstable40Gstat, DNstable40Gstat, ALstable40Gstat, TCstable40Gstat,
                           THstable40Gstat, SHstable40Gstat, DMstable40Gstat, TAstable40Gstat, MRstable40Gstat)
stable40Gstat.dat$SP <- as.factor(stable40Gstat.dat$SP)
nspp <- length(table(stable40Gstat.dat$SP))

mu <- stable40Gstat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(stable40Gstat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(stable40Gstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(stable40Gstat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(stable40Gstat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("stable40Gstatdistrib.RData")
