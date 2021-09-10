#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

##########################################################
## stationarity measurements
##########################################################

#####################################
##  PULSE DISTURBANCE (90% at 20G) ##
#####################################

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
load("DPpulse.RData")
load("PApulse.RData")
load("ZTpulse.RData")
load("PHpulse.RData")
load("VUpulse.RData")
load("PGpulse.RData")
load("SSpulse.RData")
load("PTpulse.RData")
load("SOpulse.RData")
load("MNpulse.RData")
load("ORpulse.RData")
load("NRpulse.RData")
load("GNpulse.RData")
load("DNpulse.RData")
load("ALpulse.RData")
load("TCpulse.RData")
load("THpulse.RData")
load("SHpulse.RData")
load("DMpulse.RData")
load("TApulse.RData")
load("MRpulse.RData")

# DP
DP.Nmat1 <- DP.n.sums.matb; DP.zeroFix <- which(rowSums(DP.Nmat1)==0)
if (length(DP.zeroFix)==0) {
  DP.Nmat2 <- DP.Nmat1
} else {
  DP.Nmat2 <- DP.Nmat1[-DP.zeroFix,]
}
DP.pulse.iter <- dim(DP.Nmat2)[1]
DP.pulse.RT <- apply(DP.Nmat2, MARGIN=1, return.time)
DP.pulse.MRT <- DP.pulse.VRT <- DP.pulse.stat <- rep(NA,DP.pulse.iter)
for(i in 1:DP.pulse.iter) {
  DP.pulse.MRT[i] <- DP.pulse.RT[[i]]$mRT
  DP.pulse.VRT[i] <- DP.pulse.RT[[i]]$vRT
  DP.pulse.stat[i] <- DP.pulse.MRT[i] / DP.pulse.VRT[i]
}
hist(DP.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DP.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PA
PA.Nmat1 <- PA.n.sums.matb; PA.zeroFix <- which(rowSums(PA.Nmat1)==0)
if (length(PA.zeroFix)==0) {
  PA.Nmat2 <- PA.Nmat1
} else {
  PA.Nmat2 <- PA.Nmat1[-PA.zeroFix,]
}
PA.pulse.iter <- dim(PA.Nmat2)[1]
PA.pulse.RT <- apply(PA.Nmat2, MARGIN=1, return.time)
PA.pulse.MRT <- PA.pulse.VRT <- PA.pulse.stat <- rep(NA,PA.pulse.iter)
for(i in 1:PA.pulse.iter) {
  PA.pulse.MRT[i] <- PA.pulse.RT[[i]]$mRT
  PA.pulse.VRT[i] <- PA.pulse.RT[[i]]$vRT
  PA.pulse.stat[i] <- PA.pulse.MRT[i] / PA.pulse.VRT[i]
}
hist(PA.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PA.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# ZT
ZT.Nmat1 <- ZT.n.sums.matb; ZT.zeroFix <- which(rowSums(ZT.Nmat1)==0)
if (length(ZT.zeroFix)==0) {
  ZT.Nmat2 <- ZT.Nmat1
} else {
  ZT.Nmat2 <- ZT.Nmat1[-ZT.zeroFix,]
}
ZT.pulse.iter <- dim(ZT.Nmat2)[1]
ZT.pulse.RT <- apply(ZT.Nmat2, MARGIN=1, return.time)
ZT.pulse.MRT <- ZT.pulse.VRT <- ZT.pulse.stat <- rep(NA,ZT.pulse.iter)
for(i in 1:ZT.pulse.iter) {
  ZT.pulse.MRT[i] <- ZT.pulse.RT[[i]]$mRT
  ZT.pulse.VRT[i] <- ZT.pulse.RT[[i]]$vRT
  ZT.pulse.stat[i] <- ZT.pulse.MRT[i] / ZT.pulse.VRT[i]
}
hist(ZT.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(ZT.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PH
PH.Nmat1 <- PH.n.sums.matb; PH.zeroFix <- which(rowSums(PH.Nmat1)==0)
if (length(PH.zeroFix)==0) {
  PH.Nmat2 <- PH.Nmat1
} else {
  PH.Nmat2 <- PH.Nmat1[-PH.zeroFix,]
}
PH.pulse.iter <- dim(PH.Nmat2)[1]
PH.pulse.RT <- apply(PH.Nmat2, MARGIN=1, return.time)
PH.pulse.MRT <- PH.pulse.VRT <- PH.pulse.stat <- rep(NA,PH.pulse.iter)
for(i in 1:PH.pulse.iter) {
  PH.pulse.MRT[i] <- PH.pulse.RT[[i]]$mRT
  PH.pulse.VRT[i] <- PH.pulse.RT[[i]]$vRT
  PH.pulse.stat[i] <- PH.pulse.MRT[i] / PH.pulse.VRT[i]
}
hist(PH.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PH.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# VU
VU.Nmat1 <- VU.n.sums.matb; VU.zeroFix <- which(rowSums(VU.Nmat1)==0)
if (length(VU.zeroFix)==0) {
  VU.Nmat2 <- VU.Nmat1
} else {
  VU.Nmat2 <- VU.Nmat1[-VU.zeroFix,]
}
VU.pulse.iter <- dim(VU.Nmat2)[1]
VU.pulse.RT <- apply(VU.Nmat2, MARGIN=1, return.time)
VU.pulse.MRT <- VU.pulse.VRT <- VU.pulse.stat <- rep(NA,VU.pulse.iter)
for(i in 1:VU.pulse.iter) {
  VU.pulse.MRT[i] <- VU.pulse.RT[[i]]$mRT
  VU.pulse.VRT[i] <- VU.pulse.RT[[i]]$vRT
  VU.pulse.stat[i] <- VU.pulse.MRT[i] / VU.pulse.VRT[i]
}
hist(VU.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(VU.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PG
PG.Nmat1 <- PG.n.sums.matb; PG.zeroFix <- which(rowSums(PG.Nmat1)==0)
if (length(PG.zeroFix)==0) {
  PG.Nmat2 <- PG.Nmat1
} else {
  PG.Nmat2 <- PG.Nmat1[-PG.zeroFix,]
}
PG.pulse.iter <- dim(PG.Nmat2)[1]
PG.pulse.RT <- apply(PG.Nmat2, MARGIN=1, return.time)
PG.pulse.MRT <- PG.pulse.VRT <- PG.pulse.stat <- rep(NA,PG.pulse.iter)
for(i in 1:PG.pulse.iter) {
  PG.pulse.MRT[i] <- PG.pulse.RT[[i]]$mRT
  PG.pulse.VRT[i] <- PG.pulse.RT[[i]]$vRT
  PG.pulse.stat[i] <- PG.pulse.MRT[i] / PG.pulse.VRT[i]
}
hist(PG.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PG.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SS
SS.Nmat1 <- SS.n.sums.matb; SS.zeroFix <- which(rowSums(SS.Nmat1)==0)
if (length(SS.zeroFix)==0) {
  SS.Nmat2 <- SS.Nmat1
} else {
  SS.Nmat2 <- SS.Nmat1[-SS.zeroFix,]
}
SS.pulse.iter <- dim(SS.Nmat2)[1]
SS.pulse.RT <- apply(SS.Nmat2, MARGIN=1, return.time)
SS.pulse.MRT <- SS.pulse.VRT <- SS.pulse.stat <- rep(NA,SS.pulse.iter)
for(i in 1:SS.pulse.iter) {
  SS.pulse.MRT[i] <- SS.pulse.RT[[i]]$mRT
  SS.pulse.VRT[i] <- SS.pulse.RT[[i]]$vRT
  SS.pulse.stat[i] <- SS.pulse.MRT[i] / SS.pulse.VRT[i]
}
hist(SS.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SS.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# PT
PT.Nmat1 <- PT.n.sums.matb; PT.zeroFix <- which(rowSums(PT.Nmat1)==0)
if (length(PT.zeroFix)==0) {
  PT.Nmat2 <- PT.Nmat1
} else {
  PT.Nmat2 <- PT.Nmat1[-PT.zeroFix,]
}
PT.pulse.iter <- dim(PT.Nmat2)[1]
PT.pulse.RT <- apply(PT.Nmat2, MARGIN=1, return.time)
PT.pulse.MRT <- PT.pulse.VRT <- PT.pulse.stat <- rep(NA,PT.pulse.iter)
for(i in 1:PT.pulse.iter) {
  PT.pulse.MRT[i] <- PT.pulse.RT[[i]]$mRT
  PT.pulse.VRT[i] <- PT.pulse.RT[[i]]$vRT
  PT.pulse.stat[i] <- PT.pulse.MRT[i] / PT.pulse.VRT[i]
}
hist(PT.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(PT.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SO
SO.Nmat1 <- SO.n.sums.matb; SO.zeroFix <- which(rowSums(SO.Nmat1)==0)
if (length(SO.zeroFix)==0) {
  SO.Nmat2 <- SO.Nmat1
} else {
  SO.Nmat2 <- SO.Nmat1[-SO.zeroFix,]
}
SO.pulse.iter <- dim(SO.Nmat2)[1]
SO.pulse.RT <- apply(SO.Nmat2, MARGIN=1, return.time)
SO.pulse.MRT <- SO.pulse.VRT <- SO.pulse.stat <- rep(NA,SO.pulse.iter)
for(i in 1:SO.pulse.iter) {
  SO.pulse.MRT[i] <- SO.pulse.RT[[i]]$mRT
  SO.pulse.VRT[i] <- SO.pulse.RT[[i]]$vRT
  SO.pulse.stat[i] <- SO.pulse.MRT[i] / SO.pulse.VRT[i]
}
hist(SO.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SO.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MN
MN.Nmat1 <- MN.n.sums.matb; MN.zeroFix <- which(rowSums(MN.Nmat1)==0)
if (length(MN.zeroFix)==0) {
  MN.Nmat2 <- MN.Nmat1
} else {
  MN.Nmat2 <- MN.Nmat1[-MN.zeroFix,]
}
MN.pulse.iter <- dim(MN.Nmat2)[1]
MN.pulse.RT <- apply(MN.Nmat2, MARGIN=1, return.time)
MN.pulse.MRT <- MN.pulse.VRT <- MN.pulse.stat <- rep(NA,MN.pulse.iter)
for(i in 1:MN.pulse.iter) {
  MN.pulse.MRT[i] <- MN.pulse.RT[[i]]$mRT
  MN.pulse.VRT[i] <- MN.pulse.RT[[i]]$vRT
  MN.pulse.stat[i] <- MN.pulse.MRT[i] / MN.pulse.VRT[i]
}
hist(MN.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MN.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# OR
OR.Nmat1 <- OR.n.sums.matb; OR.zeroFix <- which(rowSums(OR.Nmat1)==0)
if (length(OR.zeroFix)==0) {
  OR.Nmat2 <- OR.Nmat1
} else {
  OR.Nmat2 <- OR.Nmat1[-OR.zeroFix,]
}
OR.pulse.iter <- dim(OR.Nmat2)[1]
OR.pulse.RT <- apply(OR.Nmat2, MARGIN=1, return.time)
OR.pulse.MRT <- OR.pulse.VRT <- OR.pulse.stat <- rep(NA,OR.pulse.iter)
for(i in 1:OR.pulse.iter) {
  OR.pulse.MRT[i] <- OR.pulse.RT[[i]]$mRT
  OR.pulse.VRT[i] <- OR.pulse.RT[[i]]$vRT
  OR.pulse.stat[i] <- OR.pulse.MRT[i] / OR.pulse.VRT[i]
}
hist(OR.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(OR.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# NR
NR.Nmat1 <- NR.n.sums.matb; NR.zeroFix <- which(rowSums(NR.Nmat1)==0)
if (length(NR.zeroFix)==0) {
  NR.Nmat2 <- NR.Nmat1
} else {
  NR.Nmat2 <- NR.Nmat1[-NR.zeroFix,]
}
NR.pulse.iter <- dim(NR.Nmat2)[1]
NR.pulse.RT <- apply(NR.Nmat2, MARGIN=1, return.time)
NR.pulse.MRT <- NR.pulse.VRT <- NR.pulse.stat <- rep(NA,NR.pulse.iter)
for(i in 1:NR.pulse.iter) {
  NR.pulse.MRT[i] <- NR.pulse.RT[[i]]$mRT
  NR.pulse.VRT[i] <- NR.pulse.RT[[i]]$vRT
  NR.pulse.stat[i] <- NR.pulse.MRT[i] / NR.pulse.VRT[i]
}
hist(NR.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(NR.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# GN
GN.Nmat1 <- GN.n.sums.matb; GN.zeroFix <- which(rowSums(GN.Nmat1)==0)
if (length(GN.zeroFix)==0) {
  GN.Nmat2 <- GN.Nmat1
} else {
  GN.Nmat2 <- GN.Nmat1[-GN.zeroFix,]
}
GN.pulse.iter <- dim(GN.Nmat2)[1]
GN.pulse.RT <- apply(GN.Nmat2, MARGIN=1, return.time)
GN.pulse.MRT <- GN.pulse.VRT <- GN.pulse.stat <- rep(NA,GN.pulse.iter)
for(i in 1:GN.pulse.iter) {
  GN.pulse.MRT[i] <- GN.pulse.RT[[i]]$mRT
  GN.pulse.VRT[i] <- GN.pulse.RT[[i]]$vRT
  GN.pulse.stat[i] <- GN.pulse.MRT[i] / GN.pulse.VRT[i]
}
hist(GN.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(GN.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DN
DN.Nmat1 <- DN.n.sums.matb; DN.zeroFix <- which(rowSums(DN.Nmat1)==0)
if (length(DN.zeroFix)==0) {
  DN.Nmat2 <- DN.Nmat1
} else {
  DN.Nmat2 <- DN.Nmat1[-DN.zeroFix,]
}
DN.pulse.iter <- dim(DN.Nmat2)[1]
DN.pulse.RT <- apply(DN.Nmat2, MARGIN=1, return.time)
DN.pulse.MRT <- DN.pulse.VRT <- DN.pulse.stat <- rep(NA,DN.pulse.iter)
for(i in 1:DN.pulse.iter) {
  DN.pulse.MRT[i] <- DN.pulse.RT[[i]]$mRT
  DN.pulse.VRT[i] <- DN.pulse.RT[[i]]$vRT
  DN.pulse.stat[i] <- DN.pulse.MRT[i] / DN.pulse.VRT[i]
}
hist(DN.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DN.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# AL
AL.Nmat1 <- AL.n.sums.matb; AL.zeroFix <- which(rowSums(AL.Nmat1)==0)
if (length(AL.zeroFix)==0) {
  AL.Nmat2 <- AL.Nmat1
} else {
  AL.Nmat2 <- AL.Nmat1[-AL.zeroFix,]
}
AL.pulse.iter <- dim(AL.Nmat2)[1]
AL.pulse.RT <- apply(AL.Nmat2, MARGIN=1, return.time)
AL.pulse.MRT <- AL.pulse.VRT <- AL.pulse.stat <- rep(NA,AL.pulse.iter)
for(i in 1:AL.pulse.iter) {
  AL.pulse.MRT[i] <- AL.pulse.RT[[i]]$mRT
  AL.pulse.VRT[i] <- AL.pulse.RT[[i]]$vRT
  AL.pulse.stat[i] <- AL.pulse.MRT[i] / AL.pulse.VRT[i]
}
hist(AL.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(AL.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TC
TC.Nmat1 <- TC.n.sums.matb; TC.zeroFix <- which(rowSums(TC.Nmat1)==0)
if (length(TC.zeroFix)==0) {
  TC.Nmat2 <- TC.Nmat1
} else {
  TC.Nmat2 <- TC.Nmat1[-TC.zeroFix,]
}
TC.pulse.iter <- dim(TC.Nmat2)[1]
TC.pulse.RT <- apply(TC.Nmat2, MARGIN=1, return.time)
TC.pulse.MRT <- TC.pulse.VRT <- TC.pulse.stat <- rep(NA,TC.pulse.iter)
for(i in 1:TC.pulse.iter) {
  TC.pulse.MRT[i] <- TC.pulse.RT[[i]]$mRT
  TC.pulse.VRT[i] <- TC.pulse.RT[[i]]$vRT
  TC.pulse.stat[i] <- TC.pulse.MRT[i] / TC.pulse.VRT[i]
}
hist(TC.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TC.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TH
TH.Nmat1 <- TH.n.sums.matb; TH.zeroFix <- which(rowSums(TH.Nmat1)==0)
if (length(TH.zeroFix)==0) {
  TH.Nmat2 <- TH.Nmat1
} else {
  TH.Nmat2 <- TH.Nmat1[-TH.zeroFix,]
}
TH.pulse.iter <- dim(TH.Nmat2)[1]
TH.pulse.RT <- apply(TH.Nmat2, MARGIN=1, return.time)
TH.pulse.MRT <- TH.pulse.VRT <- TH.pulse.stat <- rep(NA,TH.pulse.iter)
for(i in 1:TH.pulse.iter) {
  TH.pulse.MRT[i] <- TH.pulse.RT[[i]]$mRT
  TH.pulse.VRT[i] <- TH.pulse.RT[[i]]$vRT
  TH.pulse.stat[i] <- TH.pulse.MRT[i] / TH.pulse.VRT[i]
}
hist(TH.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TH.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# SH
SH.Nmat1 <- SH.n.sums.matb; SH.zeroFix <- which(rowSums(SH.Nmat1)==0)
if (length(SH.zeroFix)==0) {
  SH.Nmat2 <- SH.Nmat1
} else {
  SH.Nmat2 <- SH.Nmat1[-SH.zeroFix,]
}
SH.pulse.iter <- dim(SH.Nmat2)[1]
SH.pulse.RT <- apply(SH.Nmat2, MARGIN=1, return.time)
SH.pulse.MRT <- SH.pulse.VRT <- SH.pulse.stat <- rep(NA,SH.pulse.iter)
for(i in 1:SH.pulse.iter) {
  SH.pulse.MRT[i] <- SH.pulse.RT[[i]]$mRT
  SH.pulse.VRT[i] <- SH.pulse.RT[[i]]$vRT
  SH.pulse.stat[i] <- SH.pulse.MRT[i] / SH.pulse.VRT[i]
}
hist(SH.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(SH.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# DM
DM.Nmat1 <- DM.n.sums.matb; DM.zeroFix <- which(rowSums(DM.Nmat1)==0)
if (length(DM.zeroFix)==0) {
  DM.Nmat2 <- DM.Nmat1
} else {
  DM.Nmat2 <- DM.Nmat1[-DM.zeroFix,]
}
DM.pulse.iter <- dim(DM.Nmat2)[1]
DM.pulse.RT <- apply(DM.Nmat2, MARGIN=1, return.time)
DM.pulse.MRT <- DM.pulse.VRT <- DM.pulse.stat <- rep(NA,DM.pulse.iter)
for(i in 1:DM.pulse.iter) {
  DM.pulse.MRT[i] <- DM.pulse.RT[[i]]$mRT
  DM.pulse.VRT[i] <- DM.pulse.RT[[i]]$vRT
  DM.pulse.stat[i] <- DM.pulse.MRT[i] / DM.pulse.VRT[i]
}
hist(DM.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(DM.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# TA
TA.Nmat1 <- TA.n.sums.matb; TA.zeroFix <- which(rowSums(TA.Nmat1)==0)
if (length(TA.zeroFix)==0) {
  TA.Nmat2 <- TA.Nmat1
} else {
  TA.Nmat2 <- TA.Nmat1[-TA.zeroFix,]
}
TA.pulse.iter <- dim(TA.Nmat2)[1]
TA.pulse.RT <- apply(TA.Nmat2, MARGIN=1, return.time)
TA.pulse.MRT <- TA.pulse.VRT <- TA.pulse.stat <- rep(NA,TA.pulse.iter)
for(i in 1:TA.pulse.iter) {
  TA.pulse.MRT[i] <- TA.pulse.RT[[i]]$mRT
  TA.pulse.VRT[i] <- TA.pulse.RT[[i]]$vRT
  TA.pulse.stat[i] <- TA.pulse.MRT[i] / TA.pulse.VRT[i]
}
hist(TA.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(TA.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")

# MR
MR.Nmat1 <- MR.n.sums.matb; MR.zeroFix <- which(rowSums(MR.Nmat1)==0)
if (length(MR.zeroFix)==0) {
  MR.Nmat2 <- MR.Nmat1
} else {
  MR.Nmat2 <- MR.Nmat1[-MR.zeroFix,]
}
MR.pulse.iter <- dim(MR.Nmat2)[1]
MR.pulse.RT <- apply(MR.Nmat2, MARGIN=1, return.time)
MR.pulse.MRT <- MR.pulse.VRT <- MR.pulse.stat <- rep(NA,MR.pulse.iter)
for(i in 1:MR.pulse.iter) {
  MR.pulse.MRT[i] <- MR.pulse.RT[[i]]$mRT
  MR.pulse.VRT[i] <- MR.pulse.RT[[i]]$vRT
  MR.pulse.stat[i] <- MR.pulse.MRT[i] / MR.pulse.VRT[i]
}
hist(MR.pulse.MRT, main="", xlab="mean return time", ylab="")
abline(v=2, lty=2, lwd=2, col="red")
hist(MR.pulse.stat, main="", xlab="mRT/vRT", ylab="")
abline(v=1, lty=2, lwd=2, col="red")


# summaries
DP.pulse.stat.med <- median(DP.pulse.stat, na.rm=T)
PA.pulse.stat.med <- median(PA.pulse.stat, na.rm=T)
ZT.pulse.stat.med <- median(ZT.pulse.stat, na.rm=T)
PH.pulse.stat.med <- median(PH.pulse.stat, na.rm=T)
VU.pulse.stat.med <- median(VU.pulse.stat, na.rm=T)
PG.pulse.stat.med <- median(PG.pulse.stat, na.rm=T)
SS.pulse.stat.med <- median(SS.pulse.stat, na.rm=T)
PT.pulse.stat.med <- median(PT.pulse.stat, na.rm=T)
SO.pulse.stat.med <- median(SO.pulse.stat, na.rm=T)
MN.pulse.stat.med <- median(MN.pulse.stat, na.rm=T)
OR.pulse.stat.med <- median(OR.pulse.stat, na.rm=T)
NR.pulse.stat.med <- median(NR.pulse.stat, na.rm=T)
GN.pulse.stat.med <- median(GN.pulse.stat, na.rm=T)
DN.pulse.stat.med <- median(DN.pulse.stat, na.rm=T)
AL.pulse.stat.med <- median(AL.pulse.stat, na.rm=T)
TC.pulse.stat.med <- median(TC.pulse.stat, na.rm=T)
TH.pulse.stat.med <- median(TH.pulse.stat, na.rm=T)
SH.pulse.stat.med <- median(SH.pulse.stat, na.rm=T)
DM.pulse.stat.med <- median(DM.pulse.stat, na.rm=T)
TA.pulse.stat.med <- median(TA.pulse.stat, na.rm=T)
MR.pulse.stat.med <- median(MR.pulse.stat, na.rm=T)

DP.pulse.stat.up <- quantile(DP.pulse.stat, probs=0.975, na.rm=T)
PA.pulse.stat.up <- quantile(PA.pulse.stat, probs=0.975, na.rm=T)
ZT.pulse.stat.up <- quantile(ZT.pulse.stat, probs=0.975, na.rm=T)
PH.pulse.stat.up <- quantile(PH.pulse.stat, probs=0.975, na.rm=T)
VU.pulse.stat.up <- quantile(VU.pulse.stat, probs=0.975, na.rm=T)
PG.pulse.stat.up <- quantile(PG.pulse.stat, probs=0.975, na.rm=T)
SS.pulse.stat.up <- quantile(SS.pulse.stat, probs=0.975, na.rm=T)
PT.pulse.stat.up <- quantile(PT.pulse.stat, probs=0.975, na.rm=T)
SO.pulse.stat.up <- quantile(SO.pulse.stat, probs=0.975, na.rm=T)
MN.pulse.stat.up <- quantile(MN.pulse.stat, probs=0.975, na.rm=T)
OR.pulse.stat.up <- quantile(OR.pulse.stat, probs=0.975, na.rm=T)
NR.pulse.stat.up <- quantile(NR.pulse.stat, probs=0.975, na.rm=T)
GN.pulse.stat.up <- quantile(GN.pulse.stat, probs=0.975, na.rm=T)
DN.pulse.stat.up <- quantile(DN.pulse.stat, probs=0.975, na.rm=T)
AL.pulse.stat.up <- quantile(AL.pulse.stat, probs=0.975, na.rm=T)
TC.pulse.stat.up <- quantile(TC.pulse.stat, probs=0.975, na.rm=T)
TH.pulse.stat.up <- quantile(TH.pulse.stat, probs=0.975, na.rm=T)
SH.pulse.stat.up <- quantile(SH.pulse.stat, probs=0.975, na.rm=T)
DM.pulse.stat.up <- quantile(DM.pulse.stat, probs=0.975, na.rm=T)
TA.pulse.stat.up <- quantile(TA.pulse.stat, probs=0.975, na.rm=T)
MR.pulse.stat.up <- quantile(MR.pulse.stat, probs=0.975, na.rm=T)

DP.pulse.stat.lo <- quantile(DP.pulse.stat, probs=0.025, na.rm=T)
PA.pulse.stat.lo <- quantile(PA.pulse.stat, probs=0.025, na.rm=T)
ZT.pulse.stat.lo <- quantile(ZT.pulse.stat, probs=0.025, na.rm=T)
PH.pulse.stat.lo <- quantile(PH.pulse.stat, probs=0.025, na.rm=T)
VU.pulse.stat.lo <- quantile(VU.pulse.stat, probs=0.025, na.rm=T)
PG.pulse.stat.lo <- quantile(PG.pulse.stat, probs=0.025, na.rm=T)
SS.pulse.stat.lo <- quantile(SS.pulse.stat, probs=0.025, na.rm=T)
PT.pulse.stat.lo <- quantile(PT.pulse.stat, probs=0.025, na.rm=T)
SO.pulse.stat.lo <- quantile(SO.pulse.stat, probs=0.025, na.rm=T)
MN.pulse.stat.lo <- quantile(MN.pulse.stat, probs=0.025, na.rm=T)
OR.pulse.stat.lo <- quantile(OR.pulse.stat, probs=0.025, na.rm=T)
NR.pulse.stat.lo <- quantile(NR.pulse.stat, probs=0.025, na.rm=T)
GN.pulse.stat.lo <- quantile(GN.pulse.stat, probs=0.025, na.rm=T)
DN.pulse.stat.lo <- quantile(DN.pulse.stat, probs=0.025, na.rm=T)
AL.pulse.stat.lo <- quantile(AL.pulse.stat, probs=0.025, na.rm=T)
TC.pulse.stat.lo <- quantile(TC.pulse.stat, probs=0.025, na.rm=T)
TH.pulse.stat.lo <- quantile(TH.pulse.stat, probs=0.025, na.rm=T)
SH.pulse.stat.lo <- quantile(SH.pulse.stat, probs=0.025, na.rm=T)
DM.pulse.stat.lo <- quantile(DM.pulse.stat, probs=0.025, na.rm=T)
TA.pulse.stat.lo <- quantile(TA.pulse.stat, probs=0.025, na.rm=T)
MR.pulse.stat.lo <- quantile(MR.pulse.stat, probs=0.025, na.rm=T)

pulse.stat.med <- c(DP.pulse.stat.med, PA.pulse.stat.med, ZT.pulse.stat.med, PH.pulse.stat.med, VU.pulse.stat.med, PG.pulse.stat.med,
                    SS.pulse.stat.med, PT.pulse.stat.med, SO.pulse.stat.med, MN.pulse.stat.med, OR.pulse.stat.med, NR.pulse.stat.med,
                    GN.pulse.stat.med, DN.pulse.stat.med, AL.pulse.stat.med, TC.pulse.stat.med, TH.pulse.stat.med, SH.pulse.stat.med,
                    DM.pulse.stat.med, TA.pulse.stat.med, MR.pulse.stat.med)
pulse.stat.up <- c(DP.pulse.stat.up, PA.pulse.stat.up, ZT.pulse.stat.up, PH.pulse.stat.up, VU.pulse.stat.up, PG.pulse.stat.up,
                   SS.pulse.stat.up, PT.pulse.stat.up, SO.pulse.stat.up, MN.pulse.stat.up, OR.pulse.stat.up, NR.pulse.stat.up,
                   GN.pulse.stat.up, DN.pulse.stat.up, AL.pulse.stat.up, TC.pulse.stat.up, TH.pulse.stat.up, SH.pulse.stat.up,
                   DM.pulse.stat.up, TA.pulse.stat.up, MR.pulse.stat.up)
pulse.stat.lo <- c(DP.pulse.stat.lo, PA.pulse.stat.lo, ZT.pulse.stat.lo, PH.pulse.stat.lo, VU.pulse.stat.lo, PG.pulse.stat.lo,
                   SS.pulse.stat.lo, PT.pulse.stat.lo, SO.pulse.stat.lo, MN.pulse.stat.lo, OR.pulse.stat.lo, NR.pulse.stat.lo,
                   GN.pulse.stat.lo, DN.pulse.stat.lo, AL.pulse.stat.lo, TC.pulse.stat.lo, TH.pulse.stat.lo, SH.pulse.stat.lo,
                   DM.pulse.stat.lo, TA.pulse.stat.lo, MR.pulse.stat.lo)

spp.mass.vec <- c(DP.mass,PA.mass,ZT.mass,PH.mass,VU.mass,PG.mass,SS.mass,PT.mass,SO.mass,MN.mass,OR.mass,NR.mass,GN.mass,DN.mass,AL.mass,TC.mass,TH.mass,SH.mass,DM.mass,TA.mass,MR.mass)
labs.vec <- c("DP","PA","ZT","PH","VU","PG","SS","PT","SO","MN","OR","NR","GN","DN","AL","TC","TH","SH","DM","TA","MR")

plot(log10(spp.mass.vec), pulse.stat.med, pch=19, xlab="log10 mass (kg)", ylab="mRT/vRT")
pulse.dat <- data.frame(spp.mass.vec, pulse.stat.med, pulse.stat.up, pulse.stat.lo)
colnames(pulse.dat) <- c("M", "statM", "statUP", "statLO")
rownames(pulse.dat) <- labs.vec

p <- ggplot(pulse.dat, aes(x=log10(M), y=statM)) + 
  geom_point() +
  geom_errorbar(aes(ymin=statLO, ymax=statUP), width=.2)
p + labs(x="log10 mass (kg)", y ="mRT/vRT")+
  theme_classic()

## overlapping histograms themselves
# combine data frames
DPpulsestat <- data.frame(rep("DP",length(DP.pulse.stat)), DP.pulse.stat)
PApulsestat <- data.frame(rep("PA",length(PA.pulse.stat)), PA.pulse.stat)
ZTpulsestat <- data.frame(rep("ZT",length(ZT.pulse.stat)), ZT.pulse.stat)
PHpulsestat <- data.frame(rep("PH",length(PH.pulse.stat)), PH.pulse.stat)
VUpulsestat <- data.frame(rep("VU",length(VU.pulse.stat)), VU.pulse.stat)
PGpulsestat <- data.frame(rep("PG",length(PG.pulse.stat)), PG.pulse.stat)
SSpulsestat <- data.frame(rep("SS",length(SS.pulse.stat)), SS.pulse.stat)
PTpulsestat <- data.frame(rep("PT",length(PT.pulse.stat)), PT.pulse.stat)
SOpulsestat <- data.frame(rep("SO",length(SO.pulse.stat)), SO.pulse.stat)
MNpulsestat <- data.frame(rep("MN",length(MN.pulse.stat)), MN.pulse.stat)
ORpulsestat <- data.frame(rep("OR",length(OR.pulse.stat)), OR.pulse.stat)
NRpulsestat <- data.frame(rep("NR",length(NR.pulse.stat)), NR.pulse.stat)
GNpulsestat <- data.frame(rep("GN",length(GN.pulse.stat)), GN.pulse.stat)
DNpulsestat <- data.frame(rep("DN",length(DN.pulse.stat)), DN.pulse.stat)
ALpulsestat <- data.frame(rep("AL",length(AL.pulse.stat)), AL.pulse.stat)
TCpulsestat <- data.frame(rep("TC",length(TC.pulse.stat)), TC.pulse.stat)
THpulsestat <- data.frame(rep("TH",length(TH.pulse.stat)), TH.pulse.stat)
SHpulsestat <- data.frame(rep("SH",length(SH.pulse.stat)), SH.pulse.stat)
DMpulsestat <- data.frame(rep("DM",length(DM.pulse.stat)), DM.pulse.stat)
TApulsestat <- data.frame(rep("TA",length(TA.pulse.stat)), TA.pulse.stat)
MRpulsestat <- data.frame(rep("MR",length(MR.pulse.stat)), MR.pulse.stat)

colnames(DPpulsestat) <- c("SP","mRTvRT")
colnames(PApulsestat) <- c("SP","mRTvRT")
colnames(ZTpulsestat) <- c("SP","mRTvRT")
colnames(PHpulsestat) <- c("SP","mRTvRT")
colnames(VUpulsestat) <- c("SP","mRTvRT")
colnames(PGpulsestat) <- c("SP","mRTvRT")
colnames(SSpulsestat) <- c("SP","mRTvRT")
colnames(PTpulsestat) <- c("SP","mRTvRT")
colnames(SOpulsestat) <- c("SP","mRTvRT")
colnames(MNpulsestat) <- c("SP","mRTvRT")
colnames(ORpulsestat) <- c("SP","mRTvRT")
colnames(NRpulsestat) <- c("SP","mRTvRT")
colnames(GNpulsestat) <- c("SP","mRTvRT")
colnames(DNpulsestat) <- c("SP","mRTvRT")
colnames(ALpulsestat) <- c("SP","mRTvRT")
colnames(TCpulsestat) <- c("SP","mRTvRT")
colnames(THpulsestat) <- c("SP","mRTvRT")
colnames(SHpulsestat) <- c("SP","mRTvRT")
colnames(DMpulsestat) <- c("SP","mRTvRT")
colnames(TApulsestat) <- c("SP","mRTvRT")
colnames(MRpulsestat) <- c("SP","mRTvRT")

pulsestat.dat <- rbind(DPpulsestat, PApulsestat, ZTpulsestat, PHpulsestat, VUpulsestat, PGpulsestat, SSpulsestat, PTpulsestat,
                       SOpulsestat, MNpulsestat, ORpulsestat, NRpulsestat, GNpulsestat, DNpulsestat, ALpulsestat, TCpulsestat,
                       THpulsestat, SHpulsestat, DMpulsestat, TApulsestat, MRpulsestat)
pulsestat.dat$SP <- as.factor(pulsestat.dat$SP)
nspp <- length(table(pulsestat.dat$SP))

mu <- pulsestat.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(mRTvRT))
mu

theme_set(theme_ridges())
ggplot(pulsestat.dat, aes(x = mRTvRT, y = SP, show.legend=F)) +
  xlim(0, max(pulsestat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(pulsestat.dat, x = "mRTvRT", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, max(pulsestat.dat$mRTvRT)) +
  xlab("mRT/vRT") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

save.image("pulsestatdistrib.RData")

