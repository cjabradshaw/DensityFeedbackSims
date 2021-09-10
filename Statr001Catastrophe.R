#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

##########################################################
## r = -0.001
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
load("DPdeclSm.RData")
load("PAdeclSm.RData")
load("ZTdeclSm.RData")
load("PHdeclSm.RData")
load("VUdeclSm.RData")
load("PGdeclSm.RData")
load("SSdeclSm.RData")
load("PTdeclSm.RData")
load("SOdeclSm.RData")
load("MNdeclSm.RData")
load("ORdeclSm.RData")
load("NRdeclSm.RData")
load("GNdeclSm.RData")
load("DNdeclSm.RData")
load("ALdeclSm.RData")
load("TCdeclSm.RData")
load("THdeclSm.RData")
load("SHdeclSm.RData")
load("DMdeclSm.RData")
load("TAdeclSm.RData")
load("MRdeclSm.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.r001.iter <- dim(DP.Nmat2)[1]
DP.r001.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.r001.MRT <- DP.r001.VRT <- DP.r001.stat <- rep(NA,DP.r001.iter)
for(i in 1:DP.r001.iter) {
  DP.r001.MRT[i] <- DP.r001.RT[[i]]$mRT
  DP.r001.VRT[i] <- DP.r001.RT[[i]]$vRT
  DP.r001.stat[i] <- DP.r001.MRT[i] / DP.r001.VRT[i]
}
hist(DP.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.r001.iter <- dim(PA.Nmat2)[1]
PA.r001.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.r001.MRT <- PA.r001.VRT <- PA.r001.stat <- rep(NA,PA.r001.iter)
for(i in 1:PA.r001.iter) {
  PA.r001.MRT[i] <- PA.r001.RT[[i]]$mRT
  PA.r001.VRT[i] <- PA.r001.RT[[i]]$vRT
  PA.r001.stat[i] <- PA.r001.MRT[i] / PA.r001.VRT[i]
}
hist(PA.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.r001.iter <- dim(ZT.Nmat2)[1]
ZT.r001.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.r001.MRT <- ZT.r001.VRT <- ZT.r001.stat <- rep(NA,ZT.r001.iter)
for(i in 1:ZT.r001.iter) {
  ZT.r001.MRT[i] <- ZT.r001.RT[[i]]$mRT
  ZT.r001.VRT[i] <- ZT.r001.RT[[i]]$vRT
  ZT.r001.stat[i] <- ZT.r001.MRT[i] / ZT.r001.VRT[i]
}
hist(ZT.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.r001.iter <- dim(PH.Nmat2)[1]
PH.r001.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.r001.MRT <- PH.r001.VRT <- PH.r001.stat <- rep(NA,PH.r001.iter)
for(i in 1:PH.r001.iter) {
  PH.r001.MRT[i] <- PH.r001.RT[[i]]$mRT
  PH.r001.VRT[i] <- PH.r001.RT[[i]]$vRT
  PH.r001.stat[i] <- PH.r001.MRT[i] / PH.r001.VRT[i]
}
hist(PH.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.r001.iter <- dim(VU.Nmat2)[1]
VU.r001.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.r001.MRT <- VU.r001.VRT <- VU.r001.stat <- rep(NA,VU.r001.iter)
for(i in 1:VU.r001.iter) {
  VU.r001.MRT[i] <- VU.r001.RT[[i]]$mRT
  VU.r001.VRT[i] <- VU.r001.RT[[i]]$vRT
  VU.r001.stat[i] <- VU.r001.MRT[i] / VU.r001.VRT[i]
}
hist(VU.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.r001.iter <- dim(PG.Nmat2)[1]
PG.r001.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.r001.MRT <- PG.r001.VRT <- PG.r001.stat <- rep(NA,PG.r001.iter)
for(i in 1:PG.r001.iter) {
  PG.r001.MRT[i] <- PG.r001.RT[[i]]$mRT
  PG.r001.VRT[i] <- PG.r001.RT[[i]]$vRT
  PG.r001.stat[i] <- PG.r001.MRT[i] / PG.r001.VRT[i]
}
hist(PG.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.r001.iter <- dim(SS.Nmat2)[1]
SS.r001.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.r001.MRT <- SS.r001.VRT <- SS.r001.stat <- rep(NA,SS.r001.iter)
for(i in 1:SS.r001.iter) {
  SS.r001.MRT[i] <- SS.r001.RT[[i]]$mRT
  SS.r001.VRT[i] <- SS.r001.RT[[i]]$vRT
  SS.r001.stat[i] <- SS.r001.MRT[i] / SS.r001.VRT[i]
}
hist(SS.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.r001.iter <- dim(PT.Nmat2)[1]
PT.r001.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.r001.MRT <- PT.r001.VRT <- PT.r001.stat <- rep(NA,PT.r001.iter)
for(i in 1:PT.r001.iter) {
  PT.r001.MRT[i] <- PT.r001.RT[[i]]$mRT
  PT.r001.VRT[i] <- PT.r001.RT[[i]]$vRT
  PT.r001.stat[i] <- PT.r001.MRT[i] / PT.r001.VRT[i]
}
hist(PT.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.r001.iter <- dim(SO.Nmat2)[1]
SO.r001.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.r001.MRT <- SO.r001.VRT <- SO.r001.stat <- rep(NA,SO.r001.iter)
for(i in 1:SO.r001.iter) {
  SO.r001.MRT[i] <- SO.r001.RT[[i]]$mRT
  SO.r001.VRT[i] <- SO.r001.RT[[i]]$vRT
  SO.r001.stat[i] <- SO.r001.MRT[i] / SO.r001.VRT[i]
}
hist(SO.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.r001.iter <- dim(MN.Nmat2)[1]
MN.r001.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.r001.MRT <- MN.r001.VRT <- MN.r001.stat <- rep(NA,MN.r001.iter)
for(i in 1:MN.r001.iter) {
  MN.r001.MRT[i] <- MN.r001.RT[[i]]$mRT
  MN.r001.VRT[i] <- MN.r001.RT[[i]]$vRT
  MN.r001.stat[i] <- MN.r001.MRT[i] / MN.r001.VRT[i]
}
hist(MN.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.r001.iter <- dim(OR.Nmat2)[1]
OR.r001.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.r001.MRT <- OR.r001.VRT <- OR.r001.stat <- rep(NA,OR.r001.iter)
for(i in 1:OR.r001.iter) {
  OR.r001.MRT[i] <- OR.r001.RT[[i]]$mRT
  OR.r001.VRT[i] <- OR.r001.RT[[i]]$vRT
  OR.r001.stat[i] <- OR.r001.MRT[i] / OR.r001.VRT[i]
}
hist(OR.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.r001.iter <- dim(NR.Nmat2)[1]
NR.r001.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.r001.MRT <- NR.r001.VRT <- NR.r001.stat <- rep(NA,NR.r001.iter)
for(i in 1:NR.r001.iter) {
  NR.r001.MRT[i] <- NR.r001.RT[[i]]$mRT
  NR.r001.VRT[i] <- NR.r001.RT[[i]]$vRT
  NR.r001.stat[i] <- NR.r001.MRT[i] / NR.r001.VRT[i]
}
hist(NR.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.r001.iter <- dim(GN.Nmat2)[1]
GN.r001.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.r001.MRT <- GN.r001.VRT <- GN.r001.stat <- rep(NA,GN.r001.iter)
for(i in 1:GN.r001.iter) {
  GN.r001.MRT[i] <- GN.r001.RT[[i]]$mRT
  GN.r001.VRT[i] <- GN.r001.RT[[i]]$vRT
  GN.r001.stat[i] <- GN.r001.MRT[i] / GN.r001.VRT[i]
}
hist(GN.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.r001.iter <- dim(DN.Nmat2)[1]
DN.r001.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.r001.MRT <- DN.r001.VRT <- DN.r001.stat <- rep(NA,DN.r001.iter)
for(i in 1:DN.r001.iter) {
  DN.r001.MRT[i] <- DN.r001.RT[[i]]$mRT
  DN.r001.VRT[i] <- DN.r001.RT[[i]]$vRT
  DN.r001.stat[i] <- DN.r001.MRT[i] / DN.r001.VRT[i]
}
hist(DN.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.r001.iter <- dim(AL.Nmat2)[1]
AL.r001.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.r001.MRT <- AL.r001.VRT <- AL.r001.stat <- rep(NA,AL.r001.iter)
for(i in 1:AL.r001.iter) {
  AL.r001.MRT[i] <- AL.r001.RT[[i]]$mRT
  AL.r001.VRT[i] <- AL.r001.RT[[i]]$vRT
  AL.r001.stat[i] <- AL.r001.MRT[i] / AL.r001.VRT[i]
}
hist(AL.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.r001.iter <- dim(TC.Nmat2)[1]
TC.r001.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.r001.MRT <- TC.r001.VRT <- TC.r001.stat <- rep(NA,TC.r001.iter)
for(i in 1:TC.r001.iter) {
  TC.r001.MRT[i] <- TC.r001.RT[[i]]$mRT
  TC.r001.VRT[i] <- TC.r001.RT[[i]]$vRT
  TC.r001.stat[i] <- TC.r001.MRT[i] / TC.r001.VRT[i]
}
hist(TC.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.r001.iter <- dim(TH.Nmat2)[1]
TH.r001.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.r001.MRT <- TH.r001.VRT <- TH.r001.stat <- rep(NA,TH.r001.iter)
for(i in 1:TH.r001.iter) {
  TH.r001.MRT[i] <- TH.r001.RT[[i]]$mRT
  TH.r001.VRT[i] <- TH.r001.RT[[i]]$vRT
  TH.r001.stat[i] <- TH.r001.MRT[i] / TH.r001.VRT[i]
}
hist(TH.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.r001.iter <- dim(SH.Nmat2)[1]
SH.r001.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.r001.MRT <- SH.r001.VRT <- SH.r001.stat <- rep(NA,SH.r001.iter)
for(i in 1:SH.r001.iter) {
  SH.r001.MRT[i] <- SH.r001.RT[[i]]$mRT
  SH.r001.VRT[i] <- SH.r001.RT[[i]]$vRT
  SH.r001.stat[i] <- SH.r001.MRT[i] / SH.r001.VRT[i]
}
hist(SH.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.r001.iter <- dim(DM.Nmat2)[1]
DM.r001.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.r001.MRT <- DM.r001.VRT <- DM.r001.stat <- rep(NA,DM.r001.iter)
for(i in 1:DM.r001.iter) {
  DM.r001.MRT[i] <- DM.r001.RT[[i]]$mRT
  DM.r001.VRT[i] <- DM.r001.RT[[i]]$vRT
  DM.r001.stat[i] <- DM.r001.MRT[i] / DM.r001.VRT[i]
}
hist(DM.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.r001.iter <- dim(TA.Nmat2)[1]
TA.r001.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.r001.MRT <- TA.r001.VRT <- TA.r001.stat <- rep(NA,TA.r001.iter)
for(i in 1:TA.r001.iter) {
  TA.r001.MRT[i] <- TA.r001.RT[[i]]$mRT
  TA.r001.VRT[i] <- TA.r001.RT[[i]]$vRT
  TA.r001.stat[i] <- TA.r001.MRT[i] / TA.r001.VRT[i]
}
hist(TA.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.r001.iter <- dim(MR.Nmat2)[1]
MR.r001.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.r001.MRT <- MR.r001.VRT <- MR.r001.stat <- rep(NA,MR.r001.iter)
for(i in 1:MR.r001.iter) {
  MR.r001.MRT[i] <- MR.r001.RT[[i]]$mRT
  MR.r001.VRT[i] <- MR.r001.RT[[i]]$vRT
  MR.r001.stat[i] <- MR.r001.MRT[i] / MR.r001.VRT[i]
}
hist(MR.r001.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.r001.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.r001.stat.med <- median(DP.r001.stat, na.rm=T)
PA.r001.stat.med <- median(PA.r001.stat, na.rm=T)
ZT.r001.stat.med <- median(ZT.r001.stat, na.rm=T)
PH.r001.stat.med <- median(PH.r001.stat, na.rm=T)
VU.r001.stat.med <- median(VU.r001.stat, na.rm=T)
PG.r001.stat.med <- median(PG.r001.stat, na.rm=T)
SS.r001.stat.med <- median(SS.r001.stat, na.rm=T)
PT.r001.stat.med <- median(PT.r001.stat, na.rm=T)
SO.r001.stat.med <- median(SO.r001.stat, na.rm=T)
MN.r001.stat.med <- median(MN.r001.stat, na.rm=T)
OR.r001.stat.med <- median(OR.r001.stat, na.rm=T)
NR.r001.stat.med <- median(NR.r001.stat, na.rm=T)
GN.r001.stat.med <- median(GN.r001.stat, na.rm=T)
DN.r001.stat.med <- median(DN.r001.stat, na.rm=T)
AL.r001.stat.med <- median(AL.r001.stat, na.rm=T)
TC.r001.stat.med <- median(TC.r001.stat, na.rm=T)
TH.r001.stat.med <- median(TH.r001.stat, na.rm=T)
SH.r001.stat.med <- median(SH.r001.stat, na.rm=T)
DM.r001.stat.med <- median(DM.r001.stat, na.rm=T)
TA.r001.stat.med <- median(TA.r001.stat, na.rm=T)
MR.r001.stat.med <- median(MR.r001.stat, na.rm=T)

DP.r001.stat.up <- quantile(DP.r001.stat, probs=0.975, na.rm=T)
PA.r001.stat.up <- quantile(PA.r001.stat, probs=0.975, na.rm=T)
ZT.r001.stat.up <- quantile(ZT.r001.stat, probs=0.975, na.rm=T)
PH.r001.stat.up <- quantile(PH.r001.stat, probs=0.975, na.rm=T)
VU.r001.stat.up <- quantile(VU.r001.stat, probs=0.975, na.rm=T)
PG.r001.stat.up <- quantile(PG.r001.stat, probs=0.975, na.rm=T)
SS.r001.stat.up <- quantile(SS.r001.stat, probs=0.975, na.rm=T)
PT.r001.stat.up <- quantile(PT.r001.stat, probs=0.975, na.rm=T)
SO.r001.stat.up <- quantile(SO.r001.stat, probs=0.975, na.rm=T)
MN.r001.stat.up <- quantile(MN.r001.stat, probs=0.975, na.rm=T)
OR.r001.stat.up <- quantile(OR.r001.stat, probs=0.975, na.rm=T)
NR.r001.stat.up <- quantile(NR.r001.stat, probs=0.975, na.rm=T)
GN.r001.stat.up <- quantile(GN.r001.stat, probs=0.975, na.rm=T)
DN.r001.stat.up <- quantile(DN.r001.stat, probs=0.975, na.rm=T)
AL.r001.stat.up <- quantile(AL.r001.stat, probs=0.975, na.rm=T)
TC.r001.stat.up <- quantile(TC.r001.stat, probs=0.975, na.rm=T)
TH.r001.stat.up <- quantile(TH.r001.stat, probs=0.975, na.rm=T)
SH.r001.stat.up <- quantile(SH.r001.stat, probs=0.975, na.rm=T)
DM.r001.stat.up <- quantile(DM.r001.stat, probs=0.975, na.rm=T)
TA.r001.stat.up <- quantile(TA.r001.stat, probs=0.975, na.rm=T)
MR.r001.stat.up <- quantile(MR.r001.stat, probs=0.975, na.rm=T)

DP.r001.stat.lo <- quantile(DP.r001.stat, probs=0.025, na.rm=T)
PA.r001.stat.lo <- quantile(PA.r001.stat, probs=0.025, na.rm=T)
ZT.r001.stat.lo <- quantile(ZT.r001.stat, probs=0.025, na.rm=T)
PH.r001.stat.lo <- quantile(PH.r001.stat, probs=0.025, na.rm=T)
VU.r001.stat.lo <- quantile(VU.r001.stat, probs=0.025, na.rm=T)
PG.r001.stat.lo <- quantile(PG.r001.stat, probs=0.025, na.rm=T)
SS.r001.stat.lo <- quantile(SS.r001.stat, probs=0.025, na.rm=T)
PT.r001.stat.lo <- quantile(PT.r001.stat, probs=0.025, na.rm=T)
SO.r001.stat.lo <- quantile(SO.r001.stat, probs=0.025, na.rm=T)
MN.r001.stat.lo <- quantile(MN.r001.stat, probs=0.025, na.rm=T)
OR.r001.stat.lo <- quantile(OR.r001.stat, probs=0.025, na.rm=T)
NR.r001.stat.lo <- quantile(NR.r001.stat, probs=0.025, na.rm=T)
GN.r001.stat.lo <- quantile(GN.r001.stat, probs=0.025, na.rm=T)
DN.r001.stat.lo <- quantile(DN.r001.stat, probs=0.025, na.rm=T)
AL.r001.stat.lo <- quantile(AL.r001.stat, probs=0.025, na.rm=T)
TC.r001.stat.lo <- quantile(TC.r001.stat, probs=0.025, na.rm=T)
TH.r001.stat.lo <- quantile(TH.r001.stat, probs=0.025, na.rm=T)
SH.r001.stat.lo <- quantile(SH.r001.stat, probs=0.025, na.rm=T)
DM.r001.stat.lo <- quantile(DM.r001.stat, probs=0.025, na.rm=T)
TA.r001.stat.lo <- quantile(TA.r001.stat, probs=0.025, na.rm=T)
MR.r001.stat.lo <- quantile(MR.r001.stat, probs=0.025, na.rm=T)

r001.stat.med <- c(DP.r001.stat.med, PA.r001.stat.med, ZT.r001.stat.med, PH.r001.stat.med, VU.r001.stat.med, PG.r001.stat.med,
                   SS.r001.stat.med, PT.r001.stat.med, SO.r001.stat.med, MN.r001.stat.med, OR.r001.stat.med, NR.r001.stat.med,
                   GN.r001.stat.med, DN.r001.stat.med, AL.r001.stat.med, TC.r001.stat.med, TH.r001.stat.med, SH.r001.stat.med,
                   DM.r001.stat.med, TA.r001.stat.med, MR.r001.stat.med)
r001.stat.up <- c(DP.r001.stat.up, PA.r001.stat.up, ZT.r001.stat.up, PH.r001.stat.up, VU.r001.stat.up, PG.r001.stat.up,
                  SS.r001.stat.up, PT.r001.stat.up, SO.r001.stat.up, MN.r001.stat.up, OR.r001.stat.up, NR.r001.stat.up,
                  GN.r001.stat.up, DN.r001.stat.up, AL.r001.stat.up, TC.r001.stat.up, TH.r001.stat.up, SH.r001.stat.up,
                  DM.r001.stat.up, TA.r001.stat.up, MR.r001.stat.up)
r001.stat.lo <- c(DP.r001.stat.lo, PA.r001.stat.lo, ZT.r001.stat.lo, PH.r001.stat.lo, VU.r001.stat.lo, PG.r001.stat.lo,
                  SS.r001.stat.lo, PT.r001.stat.lo, SO.r001.stat.lo, MN.r001.stat.lo, OR.r001.stat.lo, NR.r001.stat.lo,
                  GN.r001.stat.lo, DN.r001.stat.lo, AL.r001.stat.lo, TC.r001.stat.lo, TH.r001.stat.lo, SH.r001.stat.lo,
                  DM.r001.stat.lo, TA.r001.stat.lo, MR.r001.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), r001.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
r001.dat <- data.frame(spp.mass.vec, r001.stat.med, r001.stat.up, r001.stat.lo)
colnames(r001.dat) <- c("M", "statM", "statUP", "statLO")
rownames(r001.dat) <- labs.vec

p <- ggplot(r001.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPr001stat <- data.frame(rep("DP",length(DP.r001.stat)), DP.r001.stat)
PAr001stat <- data.frame(rep("PA",length(PA.r001.stat)), PA.r001.stat)
ZTr001stat <- data.frame(rep("ZT",length(ZT.r001.stat)), ZT.r001.stat)
PHr001stat <- data.frame(rep("PH",length(PH.r001.stat)), PH.r001.stat)
VUr001stat <- data.frame(rep("VU",length(VU.r001.stat)), VU.r001.stat)
PGr001stat <- data.frame(rep("PG",length(PG.r001.stat)), PG.r001.stat)
SSr001stat <- data.frame(rep("SS",length(SS.r001.stat)), SS.r001.stat)
PTr001stat <- data.frame(rep("PT",length(PT.r001.stat)), PT.r001.stat)
SOr001stat <- data.frame(rep("SO",length(SO.r001.stat)), SO.r001.stat)
MNr001stat <- data.frame(rep("MN",length(MN.r001.stat)), MN.r001.stat)
ORr001stat <- data.frame(rep("OR",length(OR.r001.stat)), OR.r001.stat)
NRr001stat <- data.frame(rep("NR",length(NR.r001.stat)), NR.r001.stat)
GNr001stat <- data.frame(rep("GN",length(GN.r001.stat)), GN.r001.stat)
DNr001stat <- data.frame(rep("DN",length(DN.r001.stat)), DN.r001.stat)
ALr001stat <- data.frame(rep("AL",length(AL.r001.stat)), AL.r001.stat)
TCr001stat <- data.frame(rep("TC",length(TC.r001.stat)), TC.r001.stat)
THr001stat <- data.frame(rep("TH",length(TH.r001.stat)), TH.r001.stat)
SHr001stat <- data.frame(rep("SH",length(SH.r001.stat)), SH.r001.stat)
DMr001stat <- data.frame(rep("DM",length(DM.r001.stat)), DM.r001.stat)
TAr001stat <- data.frame(rep("TA",length(TA.r001.stat)), TA.r001.stat)
MRr001stat <- data.frame(rep("MR",length(MR.r001.stat)), MR.r001.stat)

colnames(DPr001stat) <- c("SP","mRTvRT")
colnames(PAr001stat) <- c("SP","mRTvRT")
colnames(ZTr001stat) <- c("SP","mRTvRT")
colnames(PHr001stat) <- c("SP","mRTvRT")
colnames(VUr001stat) <- c("SP","mRTvRT")
colnames(PGr001stat) <- c("SP","mRTvRT")
colnames(SSr001stat) <- c("SP","mRTvRT")
colnames(PTr001stat) <- c("SP","mRTvRT")
colnames(SOr001stat) <- c("SP","mRTvRT")
colnames(MNr001stat) <- c("SP","mRTvRT")
colnames(ORr001stat) <- c("SP","mRTvRT")
colnames(NRr001stat) <- c("SP","mRTvRT")
colnames(GNr001stat) <- c("SP","mRTvRT")
colnames(DNr001stat) <- c("SP","mRTvRT")
colnames(ALr001stat) <- c("SP","mRTvRT")
colnames(TCr001stat) <- c("SP","mRTvRT")
colnames(THr001stat) <- c("SP","mRTvRT")
colnames(SHr001stat) <- c("SP","mRTvRT")
colnames(DMr001stat) <- c("SP","mRTvRT")
colnames(TAr001stat) <- c("SP","mRTvRT")
colnames(MRr001stat) <- c("SP","mRTvRT")

r001stat.dat <- rbind(DPr001stat, PAr001stat, ZTr001stat, PHr001stat, VUr001stat, PGr001stat, SSr001stat, PTr001stat,
                      SOr001stat, MNr001stat, ORr001stat, NRr001stat, GNr001stat, DNr001stat, ALr001stat, TCr001stat,
                      THr001stat, SHr001stat, DMr001stat, TAr001stat, MRr001stat)
r001stat.dat$SP <- as.factor(r001stat.dat$SP)
nspp <- length(table(r001stat.dat$SP))

mu <- r001stat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(r001stat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(r001stat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(r001stat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(r001stat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("r001statdistrib.RData")

