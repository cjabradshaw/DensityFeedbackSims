#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

##########################################################
## r = -0.01
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
load("DPdecl.RData")
load("PAdecl.RData")
load("ZTdecl.RData")
load("PHdecl.RData")
load("VUdecl.RData")
load("PGdecl.RData")
load("SSdecl.RData")
load("PTdecl.RData")
load("SOdecl.RData")
load("MNdecl.RData")
load("ORdecl.RData")
load("NRdecl.RData")
load("GNdecl.RData")
load("DNdecl.RData")
load("ALdecl.RData")
load("TCdecl.RData")
load("THdecl.RData")
load("SHdecl.RData")
load("DMdecl.RData")
load("TAdecl.RData")
load("MRdecl.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.r01.iter <- dim(DP.Nmat2)[1]
DP.r01.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.r01.MRT <- DP.r01.VRT <- DP.r01.stat <- rep(NA,DP.r01.iter)
for(i in 1:DP.r01.iter) {
  DP.r01.MRT[i] <- DP.r01.RT[[i]]$mRT
  DP.r01.VRT[i] <- DP.r01.RT[[i]]$vRT
  DP.r01.stat[i] <- DP.r01.MRT[i] / DP.r01.VRT[i]
}
hist(DP.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.r01.iter <- dim(PA.Nmat2)[1]
PA.r01.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.r01.MRT <- PA.r01.VRT <- PA.r01.stat <- rep(NA,PA.r01.iter)
for(i in 1:PA.r01.iter) {
  PA.r01.MRT[i] <- PA.r01.RT[[i]]$mRT
  PA.r01.VRT[i] <- PA.r01.RT[[i]]$vRT
  PA.r01.stat[i] <- PA.r01.MRT[i] / PA.r01.VRT[i]
}
hist(PA.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.r01.iter <- dim(ZT.Nmat2)[1]
ZT.r01.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.r01.MRT <- ZT.r01.VRT <- ZT.r01.stat <- rep(NA,ZT.r01.iter)
for(i in 1:ZT.r01.iter) {
  ZT.r01.MRT[i] <- ZT.r01.RT[[i]]$mRT
  ZT.r01.VRT[i] <- ZT.r01.RT[[i]]$vRT
  ZT.r01.stat[i] <- ZT.r01.MRT[i] / ZT.r01.VRT[i]
}
hist(ZT.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.r01.iter <- dim(PH.Nmat2)[1]
PH.r01.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.r01.MRT <- PH.r01.VRT <- PH.r01.stat <- rep(NA,PH.r01.iter)
for(i in 1:PH.r01.iter) {
  PH.r01.MRT[i] <- PH.r01.RT[[i]]$mRT
  PH.r01.VRT[i] <- PH.r01.RT[[i]]$vRT
  PH.r01.stat[i] <- PH.r01.MRT[i] / PH.r01.VRT[i]
}
hist(PH.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.r01.iter <- dim(VU.Nmat2)[1]
VU.r01.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.r01.MRT <- VU.r01.VRT <- VU.r01.stat <- rep(NA,VU.r01.iter)
for(i in 1:VU.r01.iter) {
  VU.r01.MRT[i] <- VU.r01.RT[[i]]$mRT
  VU.r01.VRT[i] <- VU.r01.RT[[i]]$vRT
  VU.r01.stat[i] <- VU.r01.MRT[i] / VU.r01.VRT[i]
}
hist(VU.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.r01.iter <- dim(PG.Nmat2)[1]
PG.r01.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.r01.MRT <- PG.r01.VRT <- PG.r01.stat <- rep(NA,PG.r01.iter)
for(i in 1:PG.r01.iter) {
  PG.r01.MRT[i] <- PG.r01.RT[[i]]$mRT
  PG.r01.VRT[i] <- PG.r01.RT[[i]]$vRT
  PG.r01.stat[i] <- PG.r01.MRT[i] / PG.r01.VRT[i]
}
hist(PG.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.r01.iter <- dim(SS.Nmat2)[1]
SS.r01.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.r01.MRT <- SS.r01.VRT <- SS.r01.stat <- rep(NA,SS.r01.iter)
for(i in 1:SS.r01.iter) {
  SS.r01.MRT[i] <- SS.r01.RT[[i]]$mRT
  SS.r01.VRT[i] <- SS.r01.RT[[i]]$vRT
  SS.r01.stat[i] <- SS.r01.MRT[i] / SS.r01.VRT[i]
}
hist(SS.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.r01.iter <- dim(PT.Nmat2)[1]
PT.r01.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.r01.MRT <- PT.r01.VRT <- PT.r01.stat <- rep(NA,PT.r01.iter)
for(i in 1:PT.r01.iter) {
  PT.r01.MRT[i] <- PT.r01.RT[[i]]$mRT
  PT.r01.VRT[i] <- PT.r01.RT[[i]]$vRT
  PT.r01.stat[i] <- PT.r01.MRT[i] / PT.r01.VRT[i]
}
hist(PT.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.r01.iter <- dim(SO.Nmat2)[1]
SO.r01.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.r01.MRT <- SO.r01.VRT <- SO.r01.stat <- rep(NA,SO.r01.iter)
for(i in 1:SO.r01.iter) {
  SO.r01.MRT[i] <- SO.r01.RT[[i]]$mRT
  SO.r01.VRT[i] <- SO.r01.RT[[i]]$vRT
  SO.r01.stat[i] <- SO.r01.MRT[i] / SO.r01.VRT[i]
}
hist(SO.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.r01.iter <- dim(MN.Nmat2)[1]
MN.r01.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.r01.MRT <- MN.r01.VRT <- MN.r01.stat <- rep(NA,MN.r01.iter)
for(i in 1:MN.r01.iter) {
  MN.r01.MRT[i] <- MN.r01.RT[[i]]$mRT
  MN.r01.VRT[i] <- MN.r01.RT[[i]]$vRT
  MN.r01.stat[i] <- MN.r01.MRT[i] / MN.r01.VRT[i]
}
hist(MN.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.r01.iter <- dim(OR.Nmat2)[1]
OR.r01.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.r01.MRT <- OR.r01.VRT <- OR.r01.stat <- rep(NA,OR.r01.iter)
for(i in 1:OR.r01.iter) {
  OR.r01.MRT[i] <- OR.r01.RT[[i]]$mRT
  OR.r01.VRT[i] <- OR.r01.RT[[i]]$vRT
  OR.r01.stat[i] <- OR.r01.MRT[i] / OR.r01.VRT[i]
}
hist(OR.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.r01.iter <- dim(NR.Nmat2)[1]
NR.r01.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.r01.MRT <- NR.r01.VRT <- NR.r01.stat <- rep(NA,NR.r01.iter)
for(i in 1:NR.r01.iter) {
  NR.r01.MRT[i] <- NR.r01.RT[[i]]$mRT
  NR.r01.VRT[i] <- NR.r01.RT[[i]]$vRT
  NR.r01.stat[i] <- NR.r01.MRT[i] / NR.r01.VRT[i]
}
hist(NR.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.r01.iter <- dim(GN.Nmat2)[1]
GN.r01.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.r01.MRT <- GN.r01.VRT <- GN.r01.stat <- rep(NA,GN.r01.iter)
for(i in 1:GN.r01.iter) {
  GN.r01.MRT[i] <- GN.r01.RT[[i]]$mRT
  GN.r01.VRT[i] <- GN.r01.RT[[i]]$vRT
  GN.r01.stat[i] <- GN.r01.MRT[i] / GN.r01.VRT[i]
}
hist(GN.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.r01.iter <- dim(DN.Nmat2)[1]
DN.r01.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.r01.MRT <- DN.r01.VRT <- DN.r01.stat <- rep(NA,DN.r01.iter)
for(i in 1:DN.r01.iter) {
  DN.r01.MRT[i] <- DN.r01.RT[[i]]$mRT
  DN.r01.VRT[i] <- DN.r01.RT[[i]]$vRT
  DN.r01.stat[i] <- DN.r01.MRT[i] / DN.r01.VRT[i]
}
hist(DN.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.r01.iter <- dim(AL.Nmat2)[1]
AL.r01.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.r01.MRT <- AL.r01.VRT <- AL.r01.stat <- rep(NA,AL.r01.iter)
for(i in 1:AL.r01.iter) {
  AL.r01.MRT[i] <- AL.r01.RT[[i]]$mRT
  AL.r01.VRT[i] <- AL.r01.RT[[i]]$vRT
  AL.r01.stat[i] <- AL.r01.MRT[i] / AL.r01.VRT[i]
}
hist(AL.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.r01.iter <- dim(TC.Nmat2)[1]
TC.r01.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.r01.MRT <- TC.r01.VRT <- TC.r01.stat <- rep(NA,TC.r01.iter)
for(i in 1:TC.r01.iter) {
  TC.r01.MRT[i] <- TC.r01.RT[[i]]$mRT
  TC.r01.VRT[i] <- TC.r01.RT[[i]]$vRT
  TC.r01.stat[i] <- TC.r01.MRT[i] / TC.r01.VRT[i]
}
hist(TC.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.r01.iter <- dim(TH.Nmat2)[1]
TH.r01.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.r01.MRT <- TH.r01.VRT <- TH.r01.stat <- rep(NA,TH.r01.iter)
for(i in 1:TH.r01.iter) {
  TH.r01.MRT[i] <- TH.r01.RT[[i]]$mRT
  TH.r01.VRT[i] <- TH.r01.RT[[i]]$vRT
  TH.r01.stat[i] <- TH.r01.MRT[i] / TH.r01.VRT[i]
}
hist(TH.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.r01.iter <- dim(SH.Nmat2)[1]
SH.r01.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.r01.MRT <- SH.r01.VRT <- SH.r01.stat <- rep(NA,SH.r01.iter)
for(i in 1:SH.r01.iter) {
  SH.r01.MRT[i] <- SH.r01.RT[[i]]$mRT
  SH.r01.VRT[i] <- SH.r01.RT[[i]]$vRT
  SH.r01.stat[i] <- SH.r01.MRT[i] / SH.r01.VRT[i]
}
hist(SH.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.r01.iter <- dim(DM.Nmat2)[1]
DM.r01.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.r01.MRT <- DM.r01.VRT <- DM.r01.stat <- rep(NA,DM.r01.iter)
for(i in 1:DM.r01.iter) {
  DM.r01.MRT[i] <- DM.r01.RT[[i]]$mRT
  DM.r01.VRT[i] <- DM.r01.RT[[i]]$vRT
  DM.r01.stat[i] <- DM.r01.MRT[i] / DM.r01.VRT[i]
}
hist(DM.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.r01.iter <- dim(TA.Nmat2)[1]
TA.r01.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.r01.MRT <- TA.r01.VRT <- TA.r01.stat <- rep(NA,TA.r01.iter)
for(i in 1:TA.r01.iter) {
  TA.r01.MRT[i] <- TA.r01.RT[[i]]$mRT
  TA.r01.VRT[i] <- TA.r01.RT[[i]]$vRT
  TA.r01.stat[i] <- TA.r01.MRT[i] / TA.r01.VRT[i]
}
hist(TA.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.r01.iter <- dim(MR.Nmat2)[1]
MR.r01.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.r01.MRT <- MR.r01.VRT <- MR.r01.stat <- rep(NA,MR.r01.iter)
for(i in 1:MR.r01.iter) {
  MR.r01.MRT[i] <- MR.r01.RT[[i]]$mRT
  MR.r01.VRT[i] <- MR.r01.RT[[i]]$vRT
  MR.r01.stat[i] <- MR.r01.MRT[i] / MR.r01.VRT[i]
}
hist(MR.r01.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.r01.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.r01.stat.med <- median(DP.r01.stat, na.rm=T)
PA.r01.stat.med <- median(PA.r01.stat, na.rm=T)
ZT.r01.stat.med <- median(ZT.r01.stat, na.rm=T)
PH.r01.stat.med <- median(PH.r01.stat, na.rm=T)
VU.r01.stat.med <- median(VU.r01.stat, na.rm=T)
PG.r01.stat.med <- median(PG.r01.stat, na.rm=T)
SS.r01.stat.med <- median(SS.r01.stat, na.rm=T)
PT.r01.stat.med <- median(PT.r01.stat, na.rm=T)
SO.r01.stat.med <- median(SO.r01.stat, na.rm=T)
MN.r01.stat.med <- median(MN.r01.stat, na.rm=T)
OR.r01.stat.med <- median(OR.r01.stat, na.rm=T)
NR.r01.stat.med <- median(NR.r01.stat, na.rm=T)
GN.r01.stat.med <- median(GN.r01.stat, na.rm=T)
DN.r01.stat.med <- median(DN.r01.stat, na.rm=T)
AL.r01.stat.med <- median(AL.r01.stat, na.rm=T)
TC.r01.stat.med <- median(TC.r01.stat, na.rm=T)
TH.r01.stat.med <- median(TH.r01.stat, na.rm=T)
SH.r01.stat.med <- median(SH.r01.stat, na.rm=T)
DM.r01.stat.med <- median(DM.r01.stat, na.rm=T)
TA.r01.stat.med <- median(TA.r01.stat, na.rm=T)
MR.r01.stat.med <- median(MR.r01.stat, na.rm=T)

DP.r01.stat.up <- quantile(DP.r01.stat, probs=0.975, na.rm=T)
PA.r01.stat.up <- quantile(PA.r01.stat, probs=0.975, na.rm=T)
ZT.r01.stat.up <- quantile(ZT.r01.stat, probs=0.975, na.rm=T)
PH.r01.stat.up <- quantile(PH.r01.stat, probs=0.975, na.rm=T)
VU.r01.stat.up <- quantile(VU.r01.stat, probs=0.975, na.rm=T)
PG.r01.stat.up <- quantile(PG.r01.stat, probs=0.975, na.rm=T)
SS.r01.stat.up <- quantile(SS.r01.stat, probs=0.975, na.rm=T)
PT.r01.stat.up <- quantile(PT.r01.stat, probs=0.975, na.rm=T)
SO.r01.stat.up <- quantile(SO.r01.stat, probs=0.975, na.rm=T)
MN.r01.stat.up <- quantile(MN.r01.stat, probs=0.975, na.rm=T)
OR.r01.stat.up <- quantile(OR.r01.stat, probs=0.975, na.rm=T)
NR.r01.stat.up <- quantile(NR.r01.stat, probs=0.975, na.rm=T)
GN.r01.stat.up <- quantile(GN.r01.stat, probs=0.975, na.rm=T)
DN.r01.stat.up <- quantile(DN.r01.stat, probs=0.975, na.rm=T)
AL.r01.stat.up <- quantile(AL.r01.stat, probs=0.975, na.rm=T)
TC.r01.stat.up <- quantile(TC.r01.stat, probs=0.975, na.rm=T)
TH.r01.stat.up <- quantile(TH.r01.stat, probs=0.975, na.rm=T)
SH.r01.stat.up <- quantile(SH.r01.stat, probs=0.975, na.rm=T)
DM.r01.stat.up <- quantile(DM.r01.stat, probs=0.975, na.rm=T)
TA.r01.stat.up <- quantile(TA.r01.stat, probs=0.975, na.rm=T)
MR.r01.stat.up <- quantile(MR.r01.stat, probs=0.975, na.rm=T)

DP.r01.stat.lo <- quantile(DP.r01.stat, probs=0.025, na.rm=T)
PA.r01.stat.lo <- quantile(PA.r01.stat, probs=0.025, na.rm=T)
ZT.r01.stat.lo <- quantile(ZT.r01.stat, probs=0.025, na.rm=T)
PH.r01.stat.lo <- quantile(PH.r01.stat, probs=0.025, na.rm=T)
VU.r01.stat.lo <- quantile(VU.r01.stat, probs=0.025, na.rm=T)
PG.r01.stat.lo <- quantile(PG.r01.stat, probs=0.025, na.rm=T)
SS.r01.stat.lo <- quantile(SS.r01.stat, probs=0.025, na.rm=T)
PT.r01.stat.lo <- quantile(PT.r01.stat, probs=0.025, na.rm=T)
SO.r01.stat.lo <- quantile(SO.r01.stat, probs=0.025, na.rm=T)
MN.r01.stat.lo <- quantile(MN.r01.stat, probs=0.025, na.rm=T)
OR.r01.stat.lo <- quantile(OR.r01.stat, probs=0.025, na.rm=T)
NR.r01.stat.lo <- quantile(NR.r01.stat, probs=0.025, na.rm=T)
GN.r01.stat.lo <- quantile(GN.r01.stat, probs=0.025, na.rm=T)
DN.r01.stat.lo <- quantile(DN.r01.stat, probs=0.025, na.rm=T)
AL.r01.stat.lo <- quantile(AL.r01.stat, probs=0.025, na.rm=T)
TC.r01.stat.lo <- quantile(TC.r01.stat, probs=0.025, na.rm=T)
TH.r01.stat.lo <- quantile(TH.r01.stat, probs=0.025, na.rm=T)
SH.r01.stat.lo <- quantile(SH.r01.stat, probs=0.025, na.rm=T)
DM.r01.stat.lo <- quantile(DM.r01.stat, probs=0.025, na.rm=T)
TA.r01.stat.lo <- quantile(TA.r01.stat, probs=0.025, na.rm=T)
MR.r01.stat.lo <- quantile(MR.r01.stat, probs=0.025, na.rm=T)

r01.stat.med <- c(DP.r01.stat.med, PA.r01.stat.med, ZT.r01.stat.med, PH.r01.stat.med, VU.r01.stat.med, PG.r01.stat.med,
                  SS.r01.stat.med, PT.r01.stat.med, SO.r01.stat.med, MN.r01.stat.med, OR.r01.stat.med, NR.r01.stat.med,
                  GN.r01.stat.med, DN.r01.stat.med, AL.r01.stat.med, TC.r01.stat.med, TH.r01.stat.med, SH.r01.stat.med,
                  DM.r01.stat.med, TA.r01.stat.med, MR.r01.stat.med)
r01.stat.up <- c(DP.r01.stat.up, PA.r01.stat.up, ZT.r01.stat.up, PH.r01.stat.up, VU.r01.stat.up, PG.r01.stat.up,
                 SS.r01.stat.up, PT.r01.stat.up, SO.r01.stat.up, MN.r01.stat.up, OR.r01.stat.up, NR.r01.stat.up,
                 GN.r01.stat.up, DN.r01.stat.up, AL.r01.stat.up, TC.r01.stat.up, TH.r01.stat.up, SH.r01.stat.up,
                 DM.r01.stat.up, TA.r01.stat.up, MR.r01.stat.up)
r01.stat.lo <- c(DP.r01.stat.lo, PA.r01.stat.lo, ZT.r01.stat.lo, PH.r01.stat.lo, VU.r01.stat.lo, PG.r01.stat.lo,
                 SS.r01.stat.lo, PT.r01.stat.lo, SO.r01.stat.lo, MN.r01.stat.lo, OR.r01.stat.lo, NR.r01.stat.lo,
                 GN.r01.stat.lo, DN.r01.stat.lo, AL.r01.stat.lo, TC.r01.stat.lo, TH.r01.stat.lo, SH.r01.stat.lo,
                 DM.r01.stat.lo, TA.r01.stat.lo, MR.r01.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), r01.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
r01.dat <- data.frame(spp.mass.vec, r01.stat.med, r01.stat.up, r01.stat.lo)
colnames(r01.dat) <- c("M", "statM", "statUP", "statLO")
rownames(r01.dat) <- labs.vec

p <- ggplot(r01.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPr01stat <- data.frame(rep("DP",length(DP.r01.stat)), DP.r01.stat)
PAr01stat <- data.frame(rep("PA",length(PA.r01.stat)), PA.r01.stat)
ZTr01stat <- data.frame(rep("ZT",length(ZT.r01.stat)), ZT.r01.stat)
PHr01stat <- data.frame(rep("PH",length(PH.r01.stat)), PH.r01.stat)
VUr01stat <- data.frame(rep("VU",length(VU.r01.stat)), VU.r01.stat)
PGr01stat <- data.frame(rep("PG",length(PG.r01.stat)), PG.r01.stat)
SSr01stat <- data.frame(rep("SS",length(SS.r01.stat)), SS.r01.stat)
PTr01stat <- data.frame(rep("PT",length(PT.r01.stat)), PT.r01.stat)
SOr01stat <- data.frame(rep("SO",length(SO.r01.stat)), SO.r01.stat)
MNr01stat <- data.frame(rep("MN",length(MN.r01.stat)), MN.r01.stat)
ORr01stat <- data.frame(rep("OR",length(OR.r01.stat)), OR.r01.stat)
NRr01stat <- data.frame(rep("NR",length(NR.r01.stat)), NR.r01.stat)
GNr01stat <- data.frame(rep("GN",length(GN.r01.stat)), GN.r01.stat)
DNr01stat <- data.frame(rep("DN",length(DN.r01.stat)), DN.r01.stat)
ALr01stat <- data.frame(rep("AL",length(AL.r01.stat)), AL.r01.stat)
TCr01stat <- data.frame(rep("TC",length(TC.r01.stat)), TC.r01.stat)
THr01stat <- data.frame(rep("TH",length(TH.r01.stat)), TH.r01.stat)
SHr01stat <- data.frame(rep("SH",length(SH.r01.stat)), SH.r01.stat)
DMr01stat <- data.frame(rep("DM",length(DM.r01.stat)), DM.r01.stat)
TAr01stat <- data.frame(rep("TA",length(TA.r01.stat)), TA.r01.stat)
MRr01stat <- data.frame(rep("MR",length(MR.r01.stat)), MR.r01.stat)

colnames(DPr01stat) <- c("SP","mRTvRT")
colnames(PAr01stat) <- c("SP","mRTvRT")
colnames(ZTr01stat) <- c("SP","mRTvRT")
colnames(PHr01stat) <- c("SP","mRTvRT")
colnames(VUr01stat) <- c("SP","mRTvRT")
colnames(PGr01stat) <- c("SP","mRTvRT")
colnames(SSr01stat) <- c("SP","mRTvRT")
colnames(PTr01stat) <- c("SP","mRTvRT")
colnames(SOr01stat) <- c("SP","mRTvRT")
colnames(MNr01stat) <- c("SP","mRTvRT")
colnames(ORr01stat) <- c("SP","mRTvRT")
colnames(NRr01stat) <- c("SP","mRTvRT")
colnames(GNr01stat) <- c("SP","mRTvRT")
colnames(DNr01stat) <- c("SP","mRTvRT")
colnames(ALr01stat) <- c("SP","mRTvRT")
colnames(TCr01stat) <- c("SP","mRTvRT")
colnames(THr01stat) <- c("SP","mRTvRT")
colnames(SHr01stat) <- c("SP","mRTvRT")
colnames(DMr01stat) <- c("SP","mRTvRT")
colnames(TAr01stat) <- c("SP","mRTvRT")
colnames(MRr01stat) <- c("SP","mRTvRT")

r01stat.dat <- rbind(DPr01stat, PAr01stat, ZTr01stat, PHr01stat, VUr01stat, PGr01stat, SSr01stat, PTr01stat,
                     SOr01stat, MNr01stat, ORr01stat, NRr01stat, GNr01stat, DNr01stat, ALr01stat, TCr01stat,
                     THr01stat, SHr01stat, DMr01stat, TAr01stat, MRr01stat)
r01stat.dat$SP <- as.factor(r01stat.dat$SP)
nspp <- length(table(r01stat.dat$SP))

mu <- r01stat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(r01stat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(r01stat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(r01stat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(r01stat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("r01statdistrib.RData")
