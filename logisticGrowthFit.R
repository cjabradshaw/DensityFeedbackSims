#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

################################
## logistic-growth model fits ##
################################

## Remove everything
rm(list = ls())
library(dplyr)
library(plotly)
library(expss)
library(car)
library(Hmisc)
library(cluster)
library(bootstrap)
library(data.table)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(biostat)
library(reshape2)

# Set functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc.glm <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]}
  return(AICc.vec)}

k.glm <- function(x) {
  if (length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff)+1)}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
chdev.glm <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]))/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object

RmatCalcFunc <- function(Nmat) {
  lenSeries <- dim(Nmat)[2]
  lenIt <- dim(Nmat)[1]
  Rmat <- matrix(NA, nrow=lenIt, ncol=(lenSeries-1))
  for (i in 1:lenIt) {
    Rmat[i,] <- log(Nmat[i,2:lenSeries] / Nmat[i,1:(lenSeries-1)])
  }
  return(Rmat)
}

# functions
RWfunc <- function(nSeries) {
  which0 <- which(nSeries==0)
  if (length(which0)==0) {
    nSeriesP <- nSeries
  } else {
    nSeriesP <- nSeries[-(which0[1]:length(nSeries))]
  }
  lenSeries <- length(nSeriesP)
  rSeries <- log(nSeriesP[2:lenSeries] / nSeriesP[1:(lenSeries-1)])
  q <- length(rSeries)
  sigma.rw <- sqrt((q*(sum((rSeries^2))/q))/(q-1))
  a.rw <- 0
  LL.rw <- (-q/2)*(log(2*pi*((sigma.rw^2)*(q-1)/q))+1)
  param.rw <- 1
  AICc.rw <- (-2*LL.rw) + ((2*1*q)/(q-1-1))
  return(list(LL=LL.rw, AICc=AICc.rw))
}

EXfunc <- function(nSeries) {
  which0 <- which(nSeries==0)
  if (length(which0)==0) {
    nSeriesP <- nSeries
  } else {
    nSeriesP <- nSeries[-(which0[1]:length(nSeries))]
  }
  lenSeries <- length(nSeriesP)
  rSeries <- log(nSeriesP[2:lenSeries] / nSeriesP[1:(lenSeries-1)])
  q <- length(rSeries)
  r.mean <- mean(rSeries,na.rm=T)
  r.var <- (sd(rSeries,na.rm=T))^2
  sigma.ex <- sqrt(r.var)
  a.ex <- r.mean
  LL.ex <- (-q/2)*(log(2*pi*((sigma.ex^2)*(q-1)/q))+1)
  AICc.ex <- (-2*LL.ex) + ((2*2*q)/(q-2-1))
  return(list(a=a.ex, LL=LL.ex, AICc=AICc.ex))
}

RLfunc <- function(nSeries) {
  which0 <- which(nSeries==0)
  if (length(which0)==0) {
    nSeriesP <- nSeries
  } else {
    nSeriesP <- nSeries[-(which0[1]:length(nSeries))]
  }
  lenSeries <- length(nSeriesP)
  rSeries <- log(nSeriesP[2:lenSeries] / nSeriesP[1:(lenSeries-1)])
  q <- length(rSeries)
  rl.fit <- lm(rSeries ~ nSeriesP[1:(lenSeries-1)])
  a.rl <- rl.fit$coeff[1]
  b.rl <- rl.fit$coeff[2]
  sigma.rl <- sqrt((q*(sum((rSeries-(a.rl+(b.rl*nSeriesP[1:(lenSeries-1)])))^2)/q))/(q-1))
  k.rl <- -a.rl/b.rl
  LL.rl <- (-q/2)*(log(2*pi*((sigma.rl^2)*(q-1)/q))+1)
  AICc.rl <- (-2*LL.rl) + ((2*3*q)/(q-3-1))
  if (a.rl < 0 || b.rl > 0) AICc.rl <- 1000000*abs(EXfunc(nSeries)$AICc)
  return(list(a=as.numeric(a.rl), b=as.numeric(b.rl), K=as.numeric(k.rl), LL=LL.rl, AICc=AICc.rl))
}

GLfunc <- function(nSeries) {
  which0 <- which(nSeries==0)
  if (length(which0)==0) {
    nSeriesP <- nSeries
  } else {
    nSeriesP <- nSeries[-(which0[1]:length(nSeries))]
  }
  lenSeries <- length(nSeriesP)
  rSeries <- log(nSeriesP[2:lenSeries] / nSeriesP[1:(lenSeries-1)])
  q <- length(rSeries)
  gl.fit <- lm(rSeries ~ log(nSeriesP[1:(lenSeries-1)]))
  a.gl <- gl.fit$coeff[1]
  b.gl <- gl.fit$coeff[2]
  sigma.gl <- sqrt((q*(sum((rSeries-(a.gl+(b.gl*log(nSeriesP[1:(lenSeries-1)]))))^2)/q))/(q-1))
  k.gl <- exp(-a.gl/b.gl)
  LL.gl <- (-q/2)*(log(2*pi*((sigma.gl^2)*(q-1)/q))+1)
  AICc.gl <- (-2*LL.gl) + ((2*3*q)/(q-3-1))
  return(list(a=as.numeric(a.gl), b=as.numeric(b.gl), K=as.numeric(k.gl), LL=LL.gl, AICc=AICc.gl))
}

## source
source("matrixOperators.r") 


#####################################################################
## choose which scenario (load only species for intended scenario) ##
#####################################################################

# stable 40 G with catastrophe
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


##  pulse
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


##  r ~ -0.001
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

##  r ~ -0.01
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

##  K stoch
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

##  K stoch increasing variance
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

##  K stoch declining
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
genL <- DP.gen.l; Nmed <- DP.n.md; Nmat1 <- DP.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- DP.pred.red.mn; Ninit <- DP.pop.found; Sb <- DP.b.lp; popmat <- DP.popmat.orig; mass <- DP.mass; prim <- DP.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
DP.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
DP.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
DP.DDw.dat <- data.frame(rep("DP",length(DP.DDwAICvec)), DP.DDwAICvec)
colnames(DP.DDw.dat) <- c("SP","DDwAIC")

DP.DDwAICmed <- median(DP.DDwAICvec, na.rm=T)
DP.DDwAICUp <- as.numeric(quantile(DP.DDwAICvec, probs=0.975, na.rm=T))
DP.DDwAICLo <- as.numeric(quantile(DP.DDwAICvec, probs=0.025, na.rm=T))

DP.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
DP.GLbMed <- -(median(DP.GLbs, na.rm=T))
DP.GLbLo <- -(as.numeric(quantile(DP.GLbs, probs=0.975, na.rm=T)))
DP.GLbUp <- -(as.numeric(quantile(DP.GLbs, probs=0.025, na.rm=T)))

DP.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
DP.RLbMed <- -(median(DP.RLbs, na.rm=T))
DP.RLbLo <- -(as.numeric(quantile(DP.RLbs, probs=0.975, na.rm=T)))
DP.RLbUp <- -(as.numeric(quantile(DP.RLbs, probs=0.025, na.rm=T)))

DP.predRedMed <- 1 - median(predRed, na.rm=T)
DP.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
DP.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
DP.predRedErr <- mean(c((DP.predRedMed-DP.predRedLo),(DP.predRedUp-DP.predRedMed)))

DP.q <- lenSeries - 1
DP.DDsStrength <- 1/(as.numeric(Sb/Ninit))
DP.longev <- round(dim(popmat)[2] - 1, 0)
DP.genL <- round(G.val(popmat, DP.longev), 2)
DP.M <- round(mass,0)
DP.alpha <- round(prim,0)
DP.out <- data.frame("DP",DP.predRedMed, DP.predRedErr, DP.GLbMed, DP.GLbUp, DP.GLbLo, DP.DDwAICmed, DP.DDwAICUp, DP.DDwAICLo, DP.q, DP.longev, DP.genL, DP.M, DP.alpha)
colnames(DP.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# PA
genL <- PA.gen.l; Nmed <- PA.n.md; Nmat1 <- PA.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- PA.pred.red.mn; Ninit <- PA.pop.found; Sb <- PA.b.lp; popmat <- PA.popmat.orig; mass <- PA.mass; prim <- PA.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
PA.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
PA.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
PA.DDw.dat <- data.frame(rep("PA",length(PA.DDwAICvec)), PA.DDwAICvec)
colnames(PA.DDw.dat) <- c("SP","DDwAIC")

PA.DDwAICmed <- median(PA.DDwAICvec, na.rm=T)
PA.DDwAICUp <- as.numeric(quantile(PA.DDwAICvec, probs=0.975, na.rm=T))
PA.DDwAICLo <- as.numeric(quantile(PA.DDwAICvec, probs=0.025, na.rm=T))

PA.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
PA.GLbMed <- -(median(PA.GLbs, na.rm=T))
PA.GLbLo <- -(as.numeric(quantile(PA.GLbs, probs=0.975, na.rm=T)))
PA.GLbUp <- -(as.numeric(quantile(PA.GLbs, probs=0.025, na.rm=T)))

PA.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
PA.RLbMed <- -(median(PA.RLbs, na.rm=T))
PA.RLbLo <- -(as.numeric(quantile(PA.RLbs, probs=0.975, na.rm=T)))
PA.RLbUp <- -(as.numeric(quantile(PA.RLbs, probs=0.025, na.rm=T)))

PA.predRedMed <- 1 - median(predRed, na.rm=T)
PA.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
PA.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
PA.predRedErr <- mean(c((PA.predRedMed-PA.predRedLo),(PA.predRedUp-PA.predRedMed)))

PA.q <- lenSeries - 1
PA.DDsStrength <- 1/(as.numeric(Sb/Ninit))
PA.longev <- round(dim(popmat)[2] - 1, 0)
PA.genL <- round(G.val(popmat, PA.longev), 2)
PA.M <- round(mass,0)
PA.alpha <- round(prim,0)
PA.out <- data.frame("PA",PA.predRedMed, PA.predRedErr, PA.GLbMed, PA.GLbUp, PA.GLbLo, PA.DDwAICmed, PA.DDwAICUp, PA.DDwAICLo, PA.q, PA.longev, PA.genL, PA.M, PA.alpha)
colnames(PA.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# ZT
genL <- ZT.gen.l; Nmed <- ZT.n.md; Nmat1 <- ZT.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- ZT.pred.red.mn; Ninit <- ZT.pop.found; Sb <- ZT.b.lp; popmat <- ZT.popmat.orig; mass <- ZT.mass; prim <- ZT.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
ZT.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
ZT.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
ZT.DDw.dat <- data.frame(rep("ZT",length(ZT.DDwAICvec)), ZT.DDwAICvec)
colnames(ZT.DDw.dat) <- c("SP","DDwAIC")

ZT.DDwAICmed <- median(ZT.DDwAICvec, na.rm=T)
ZT.DDwAICUp <- as.numeric(quantile(ZT.DDwAICvec, probs=0.975, na.rm=T))
ZT.DDwAICLo <- as.numeric(quantile(ZT.DDwAICvec, probs=0.025, na.rm=T))

ZT.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
ZT.GLbMed <- -(median(ZT.GLbs, na.rm=T))
ZT.GLbLo <- -(as.numeric(quantile(ZT.GLbs, probs=0.975, na.rm=T)))
ZT.GLbUp <- -(as.numeric(quantile(ZT.GLbs, probs=0.025, na.rm=T)))

ZT.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
ZT.RLbMed <- -(median(ZT.RLbs, na.rm=T))
ZT.RLbLo <- -(as.numeric(quantile(ZT.RLbs, probs=0.975, na.rm=T)))
ZT.RLbUp <- -(as.numeric(quantile(ZT.RLbs, probs=0.025, na.rm=T)))

ZT.predRedMed <- 1 - median(predRed, na.rm=T)
ZT.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
ZT.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
ZT.predRedErr <- mean(c((ZT.predRedMed-ZT.predRedLo),(ZT.predRedUp-ZT.predRedMed)))

ZT.q <- lenSeries - 1
ZT.DDsStrength <- 1/(as.numeric(Sb/Ninit))
ZT.longev <- round(dim(popmat)[2] - 1, 0)
ZT.genL <- round(G.val(popmat, ZT.longev), 2)
ZT.M <- round(mass,0)
ZT.alpha <- round(prim,0)
ZT.out <- data.frame("ZT",ZT.predRedMed, ZT.predRedErr, ZT.GLbMed, ZT.GLbUp, ZT.GLbLo, ZT.DDwAICmed, ZT.DDwAICUp, ZT.DDwAICLo, ZT.q, ZT.longev, ZT.genL, ZT.M, ZT.alpha)
colnames(ZT.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# PH
genL <- PH.gen.l; Nmed <- PH.n.md; Nmat1 <- PH.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- PH.pred.red.mn; Ninit <- PH.pop.found; Sb <- PH.b.lp; popmat <- PH.popmat.orig; mass <- PH.mass; prim <- PH.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
PH.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
PH.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
PH.DDw.dat <- data.frame(rep("PH",length(PH.DDwAICvec)), PH.DDwAICvec)
colnames(PH.DDw.dat) <- c("SP","DDwAIC")

PH.DDwAICmed <- median(PH.DDwAICvec, na.rm=T)
PH.DDwAICUp <- as.numeric(quantile(PH.DDwAICvec, probs=0.975, na.rm=T))
PH.DDwAICLo <- as.numeric(quantile(PH.DDwAICvec, probs=0.025, na.rm=T))

PH.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
PH.GLbMed <- -(median(PH.GLbs, na.rm=T))
PH.GLbLo <- -(as.numeric(quantile(PH.GLbs, probs=0.975, na.rm=T)))
PH.GLbUp <- -(as.numeric(quantile(PH.GLbs, probs=0.025, na.rm=T)))

PH.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
PH.RLbMed <- -(median(PH.RLbs, na.rm=T))
PH.RLbLo <- -(as.numeric(quantile(PH.RLbs, probs=0.975, na.rm=T)))
PH.RLbUp <- -(as.numeric(quantile(PH.RLbs, probs=0.025, na.rm=T)))

PH.predRedMed <- 1 - median(predRed, na.rm=T)
PH.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
PH.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
PH.predRedErr <- mean(c((PH.predRedMed-PH.predRedLo),(PH.predRedUp-PH.predRedMed)))

PH.q <- lenSeries - 1
PH.DDsStrength <- 1/(as.numeric(Sb/Ninit))
PH.longev <- round(dim(popmat)[2] - 1, 0)
PH.genL <- round(G.val(popmat, PH.longev), 2)
PH.M <- round(mass,0)
PH.alpha <- round(prim,0)
PH.out <- data.frame("PH",PH.predRedMed, PH.predRedErr, PH.GLbMed, PH.GLbUp, PH.GLbLo, PH.DDwAICmed, PH.DDwAICUp, PH.DDwAICLo, PH.q, PH.longev, PH.genL, PH.M, PH.alpha)
colnames(PH.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# VU
genL <- VU.gen.l; Nmed <- VU.n.md; Nmat1 <- VU.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- VU.pred.red.mn; Ninit <- VU.pop.found; Sb <- VU.b.lp; popmat <- VU.popmat.orig; mass <- VU.mass; prim <- VU.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
VU.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
VU.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
VU.DDw.dat <- data.frame(rep("VU",length(VU.DDwAICvec)), VU.DDwAICvec)
colnames(VU.DDw.dat) <- c("SP","DDwAIC")

VU.DDwAICmed <- median(VU.DDwAICvec, na.rm=T)
VU.DDwAICUp <- as.numeric(quantile(VU.DDwAICvec, probs=0.975, na.rm=T))
VU.DDwAICLo <- as.numeric(quantile(VU.DDwAICvec, probs=0.025, na.rm=T))

VU.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
VU.GLbMed <- -(median(VU.GLbs, na.rm=T))
VU.GLbLo <- -(as.numeric(quantile(VU.GLbs, probs=0.975, na.rm=T)))
VU.GLbUp <- -(as.numeric(quantile(VU.GLbs, probs=0.025, na.rm=T)))

VU.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
VU.RLbMed <- -(median(VU.RLbs, na.rm=T))
VU.RLbLo <- -(as.numeric(quantile(VU.RLbs, probs=0.975, na.rm=T)))
VU.RLbUp <- -(as.numeric(quantile(VU.RLbs, probs=0.025, na.rm=T)))

VU.predRedMed <- 1 - median(predRed, na.rm=T)
VU.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
VU.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
VU.predRedErr <- mean(c((VU.predRedMed-VU.predRedLo),(VU.predRedUp-VU.predRedMed)))

VU.q <- lenSeries - 1
VU.DDsStrength <- 1/(as.numeric(Sb/Ninit))
VU.longev <- round(dim(popmat)[2] - 1, 0)
VU.genL <- round(G.val(popmat, VU.longev), 2)
VU.M <- round(mass,0)
VU.alpha <- round(prim,0)
VU.out <- data.frame("VU",VU.predRedMed, VU.predRedErr, VU.GLbMed, VU.GLbUp, VU.GLbLo, VU.DDwAICmed, VU.DDwAICUp, VU.DDwAICLo, VU.q, VU.longev, VU.genL, VU.M, VU.alpha)
colnames(VU.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# PG
genL <- PG.gen.l; Nmed <- PG.n.md; Nmat1 <- PG.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- PG.pred.red.mn; Ninit <- PG.pop.found; Sb <- PG.b.lp; popmat <- PG.popmat.orig; mass <- PG.mass; prim <- PG.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
PG.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
PG.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
PG.DDw.dat <- data.frame(rep("PG",length(PG.DDwAICvec)), PG.DDwAICvec)
colnames(PG.DDw.dat) <- c("SP","DDwAIC")

PG.DDwAICmed <- median(PG.DDwAICvec, na.rm=T)
PG.DDwAICUp <- as.numeric(quantile(PG.DDwAICvec, probs=0.975, na.rm=T))
PG.DDwAICLo <- as.numeric(quantile(PG.DDwAICvec, probs=0.025, na.rm=T))

PG.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
PG.GLbMed <- -(median(PG.GLbs, na.rm=T))
PG.GLbLo <- -(as.numeric(quantile(PG.GLbs, probs=0.975, na.rm=T)))
PG.GLbUp <- -(as.numeric(quantile(PG.GLbs, probs=0.025, na.rm=T)))

PG.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
PG.RLbMed <- -(median(PG.RLbs, na.rm=T))
PG.RLbLo <- -(as.numeric(quantile(PG.RLbs, probs=0.975, na.rm=T)))
PG.RLbUp <- -(as.numeric(quantile(PG.RLbs, probs=0.025, na.rm=T)))

PG.predRedMed <- 1 - median(predRed, na.rm=T)
PG.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
PG.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
PG.predRedErr <- mean(c((PG.predRedMed-PG.predRedLo),(PG.predRedUp-PG.predRedMed)))

PG.q <- lenSeries - 1
PG.DDsStrength <- 1/(as.numeric(Sb/Ninit))
PG.longev <- round(dim(popmat)[2] - 1, 0)
PG.genL <- round(G.val(popmat, PG.longev), 2)
PG.M <- round(mass,0)
PG.alpha <- round(prim,0)
PG.out <- data.frame("PG",PG.predRedMed, PG.predRedErr, PG.GLbMed, PG.GLbUp, PG.GLbLo, PG.DDwAICmed, PG.DDwAICUp, PG.DDwAICLo, PG.q, PG.longev, PG.genL, PG.M, PG.alpha)
colnames(PG.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# SS
genL <- SS.gen.l; Nmed <- SS.n.md; Nmat1 <- SS.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- SS.pred.red.mn; Ninit <- SS.pop.found; Sb <- SS.b.lp; popmat <- SS.popmat.orig; mass <- SS.mass; prim <- SS.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
SS.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
SS.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
SS.DDw.dat <- data.frame(rep("SS",length(SS.DDwAICvec)), SS.DDwAICvec)
colnames(SS.DDw.dat) <- c("SP","DDwAIC")

SS.DDwAICmed <- median(SS.DDwAICvec, na.rm=T)
SS.DDwAICUp <- as.numeric(quantile(SS.DDwAICvec, probs=0.975, na.rm=T))
SS.DDwAICLo <- as.numeric(quantile(SS.DDwAICvec, probs=0.025, na.rm=T))

SS.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
SS.GLbMed <- -(median(SS.GLbs, na.rm=T))
SS.GLbLo <- -(as.numeric(quantile(SS.GLbs, probs=0.975, na.rm=T)))
SS.GLbUp <- -(as.numeric(quantile(SS.GLbs, probs=0.025, na.rm=T)))

SS.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
SS.RLbMed <- -(median(SS.RLbs, na.rm=T))
SS.RLbLo <- -(as.numeric(quantile(SS.RLbs, probs=0.975, na.rm=T)))
SS.RLbUp <- -(as.numeric(quantile(SS.RLbs, probs=0.025, na.rm=T)))

SS.predRedMed <- 1 - median(predRed, na.rm=T)
SS.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
SS.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
SS.predRedErr <- mean(c((SS.predRedMed-SS.predRedLo),(SS.predRedUp-SS.predRedMed)))

SS.q <- lenSeries - 1
SS.DDsStrength <- 1/(as.numeric(Sb/Ninit))
SS.longev <- round(dim(popmat)[2] - 1, 0)
SS.genL <- round(G.val(popmat, SS.longev), 2)
SS.M <- round(mass,0)
SS.alpha <- round(prim,0)
SS.out <- data.frame("SS",SS.predRedMed, SS.predRedErr, SS.GLbMed, SS.GLbUp, SS.GLbLo, SS.DDwAICmed, SS.DDwAICUp, SS.DDwAICLo, SS.q, SS.longev, SS.genL, SS.M, SS.alpha)
colnames(SS.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# PT
genL <- PT.gen.l; Nmed <- PT.n.md; Nmat1 <- PT.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- PT.pred.red.mn; Ninit <- PT.pop.found; Sb <- PT.b.lp; popmat <- PT.popmat.orig; mass <- PT.mass; prim <- PT.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
PT.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
PT.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
PT.DDw.dat <- data.frame(rep("PT",length(PT.DDwAICvec)), PT.DDwAICvec)
colnames(PT.DDw.dat) <- c("SP","DDwAIC")

PT.DDwAICmed <- median(PT.DDwAICvec, na.rm=T)
PT.DDwAICUp <- as.numeric(quantile(PT.DDwAICvec, probs=0.975, na.rm=T))
PT.DDwAICLo <- as.numeric(quantile(PT.DDwAICvec, probs=0.025, na.rm=T))

PT.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
PT.GLbMed <- -(median(PT.GLbs, na.rm=T))
PT.GLbLo <- -(as.numeric(quantile(PT.GLbs, probs=0.975, na.rm=T)))
PT.GLbUp <- -(as.numeric(quantile(PT.GLbs, probs=0.025, na.rm=T)))

PT.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
PT.RLbMed <- -(median(PT.RLbs, na.rm=T))
PT.RLbLo <- -(as.numeric(quantile(PT.RLbs, probs=0.975, na.rm=T)))
PT.RLbUp <- -(as.numeric(quantile(PT.RLbs, probs=0.025, na.rm=T)))

PT.predRedMed <- 1 - median(predRed, na.rm=T)
PT.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
PT.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
PT.predRedErr <- mean(c((PT.predRedMed-PT.predRedLo),(PT.predRedUp-PT.predRedMed)))

PT.q <- lenSeries - 1
PT.DDsStrength <- 1/(as.numeric(Sb/Ninit))
PT.longev <- round(dim(popmat)[2] - 1, 0)
PT.genL <- round(G.val(popmat, PT.longev), 2)
PT.M <- round(mass,0)
PT.alpha <- round(prim,0)
PT.out <- data.frame("PT",PT.predRedMed, PT.predRedErr, PT.GLbMed, PT.GLbUp, PT.GLbLo, PT.DDwAICmed, PT.DDwAICUp, PT.DDwAICLo, PT.q, PT.longev, PT.genL, PT.M, PT.alpha)
colnames(PT.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# SO
genL <- SO.gen.l; Nmed <- SO.n.md; Nmat1 <- SO.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- SO.pred.red.mn; Ninit <- SO.pop.found; Sb <- SO.b.lp; popmat <- SO.popmat.orig; mass <- SO.mass; prim <- SO.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
SO.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
SO.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
SO.DDw.dat <- data.frame(rep("SO",length(SO.DDwAICvec)), SO.DDwAICvec)
colnames(SO.DDw.dat) <- c("SP","DDwAIC")

SO.DDwAICmed <- median(SO.DDwAICvec, na.rm=T)
SO.DDwAICUp <- as.numeric(quantile(SO.DDwAICvec, probs=0.975, na.rm=T))
SO.DDwAICLo <- as.numeric(quantile(SO.DDwAICvec, probs=0.025, na.rm=T))

SO.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
SO.GLbMed <- -(median(SO.GLbs, na.rm=T))
SO.GLbLo <- -(as.numeric(quantile(SO.GLbs, probs=0.975, na.rm=T)))
SO.GLbUp <- -(as.numeric(quantile(SO.GLbs, probs=0.025, na.rm=T)))

SO.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
SO.RLbMed <- -(median(SO.RLbs, na.rm=T))
SO.RLbLo <- -(as.numeric(quantile(SO.RLbs, probs=0.975, na.rm=T)))
SO.RLbUp <- -(as.numeric(quantile(SO.RLbs, probs=0.025, na.rm=T)))

SO.predRedMed <- 1 - median(predRed, na.rm=T)
SO.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
SO.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
SO.predRedErr <- mean(c((SO.predRedMed-SO.predRedLo),(SO.predRedUp-SO.predRedMed)))

SO.q <- lenSeries - 1
SO.DDsStrength <- 1/(as.numeric(Sb/Ninit))
SO.longev <- round(dim(popmat)[2] - 1, 0)
SO.genL <- round(G.val(popmat, SO.longev), 2)
SO.M <- round(mass,0)
SO.alpha <- round(prim,0)
SO.out <- data.frame("SO",SO.predRedMed, SO.predRedErr, SO.GLbMed, SO.GLbUp, SO.GLbLo, SO.DDwAICmed, SO.DDwAICUp, SO.DDwAICLo, SO.q, SO.longev, SO.genL, SO.M, SO.alpha)
colnames(SO.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# MN
genL <- MN.gen.l; Nmed <- MN.n.md; Nmat1 <- MN.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- MN.pred.red.mn; Ninit <- MN.pop.found; Sb <- MN.b.lp; popmat <- MN.popmat.orig; mass <- MN.mass; prim <- MN.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
MN.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
MN.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
MN.DDw.dat <- data.frame(rep("MN",length(MN.DDwAICvec)), MN.DDwAICvec)
colnames(MN.DDw.dat) <- c("SP","DDwAIC")

MN.DDwAICmed <- median(MN.DDwAICvec, na.rm=T)
MN.DDwAICUp <- as.numeric(quantile(MN.DDwAICvec, probs=0.975, na.rm=T))
MN.DDwAICLo <- as.numeric(quantile(MN.DDwAICvec, probs=0.025, na.rm=T))

MN.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
MN.GLbMed <- -(median(MN.GLbs, na.rm=T))
MN.GLbLo <- -(as.numeric(quantile(MN.GLbs, probs=0.975, na.rm=T)))
MN.GLbUp <- -(as.numeric(quantile(MN.GLbs, probs=0.025, na.rm=T)))

MN.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
MN.RLbMed <- -(median(MN.RLbs, na.rm=T))
MN.RLbLo <- -(as.numeric(quantile(MN.RLbs, probs=0.975, na.rm=T)))
MN.RLbUp <- -(as.numeric(quantile(MN.RLbs, probs=0.025, na.rm=T)))

MN.predRedMed <- 1 - median(predRed, na.rm=T)
MN.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
MN.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
MN.predRedErr <- mean(c((MN.predRedMed-MN.predRedLo),(MN.predRedUp-MN.predRedMed)))

MN.q <- lenSeries - 1
MN.DDsStrength <- 1/(as.numeric(Sb/Ninit))
MN.longev <- round(dim(popmat)[2] - 1, 0)
MN.genL <- round(G.val(popmat, MN.longev), 2)
MN.M <- round(mass,0)
MN.alpha <- round(prim,0)
MN.out <- data.frame("MN",MN.predRedMed, MN.predRedErr, MN.GLbMed, MN.GLbUp, MN.GLbLo, MN.DDwAICmed, MN.DDwAICUp, MN.DDwAICLo, MN.q, MN.longev, MN.genL, MN.M, MN.alpha)
colnames(MN.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# OR
genL <- OR.gen.l; Nmed <- OR.n.md; Nmat1 <- OR.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- OR.pred.red.mn; Ninit <- OR.pop.found; Sb <- OR.b.lp; popmat <- OR.popmat.orig; mass <- OR.mass; prim <- OR.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
OR.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
OR.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
OR.DDw.dat <- data.frame(rep("OR",length(OR.DDwAICvec)), OR.DDwAICvec)
colnames(OR.DDw.dat) <- c("SP","DDwAIC")

OR.DDwAICmed <- median(OR.DDwAICvec, na.rm=T)
OR.DDwAICUp <- as.numeric(quantile(OR.DDwAICvec, probs=0.975, na.rm=T))
OR.DDwAICLo <- as.numeric(quantile(OR.DDwAICvec, probs=0.025, na.rm=T))

OR.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
OR.GLbMed <- -(median(OR.GLbs, na.rm=T))
OR.GLbLo <- -(as.numeric(quantile(OR.GLbs, probs=0.975, na.rm=T)))
OR.GLbUp <- -(as.numeric(quantile(OR.GLbs, probs=0.025, na.rm=T)))

OR.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
OR.RLbMed <- -(median(OR.RLbs, na.rm=T))
OR.RLbLo <- -(as.numeric(quantile(OR.RLbs, probs=0.975, na.rm=T)))
OR.RLbUp <- -(as.numeric(quantile(OR.RLbs, probs=0.025, na.rm=T)))

OR.predRedMed <- 1 - median(predRed, na.rm=T)
OR.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
OR.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
OR.predRedErr <- mean(c((OR.predRedMed-OR.predRedLo),(OR.predRedUp-OR.predRedMed)))

OR.q <- lenSeries - 1
OR.DDsStrength <- 1/(as.numeric(Sb/Ninit))
OR.longev <- round(dim(popmat)[2] - 1, 0)
OR.genL <- round(G.val(popmat, OR.longev), 2)
OR.M <- round(mass,0)
OR.alpha <- round(prim,0)
OR.out <- data.frame("OR",OR.predRedMed, OR.predRedErr, OR.GLbMed, OR.GLbUp, OR.GLbLo, OR.DDwAICmed, OR.DDwAICUp, OR.DDwAICLo, OR.q, OR.longev, OR.genL, OR.M, OR.alpha)
colnames(OR.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# NR
genL <- NR.gen.l; Nmed <- NR.n.md; Nmat1 <- NR.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- NR.pred.red.mn; Ninit <- NR.pop.found; Sb <- NR.b.lp; popmat <- NR.popmat.orig; mass <- NR.mass; prim <- NR.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
NR.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
NR.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
NR.DDw.dat <- data.frame(rep("NR",length(NR.DDwAICvec)), NR.DDwAICvec)
colnames(NR.DDw.dat) <- c("SP","DDwAIC")

NR.DDwAICmed <- median(NR.DDwAICvec, na.rm=T)
NR.DDwAICUp <- as.numeric(quantile(NR.DDwAICvec, probs=0.975, na.rm=T))
NR.DDwAICLo <- as.numeric(quantile(NR.DDwAICvec, probs=0.025, na.rm=T))

NR.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
NR.GLbMed <- -(median(NR.GLbs, na.rm=T))
NR.GLbLo <- -(as.numeric(quantile(NR.GLbs, probs=0.975, na.rm=T)))
NR.GLbUp <- -(as.numeric(quantile(NR.GLbs, probs=0.025, na.rm=T)))

NR.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
NR.RLbMed <- -(median(NR.RLbs, na.rm=T))
NR.RLbLo <- -(as.numeric(quantile(NR.RLbs, probs=0.975, na.rm=T)))
NR.RLbUp <- -(as.numeric(quantile(NR.RLbs, probs=0.025, na.rm=T)))

NR.predRedMed <- 1 - median(predRed, na.rm=T)
NR.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
NR.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
NR.predRedErr <- mean(c((NR.predRedMed-NR.predRedLo),(NR.predRedUp-NR.predRedMed)))

NR.q <- lenSeries - 1
NR.DDsStrength <- 1/(as.numeric(Sb/Ninit))
NR.longev <- round(dim(popmat)[2] - 1, 0)
NR.genL <- round(G.val(popmat, NR.longev), 2)
NR.M <- round(mass,0)
NR.alpha <- round(prim,0)
NR.out <- data.frame("NR",NR.predRedMed, NR.predRedErr, NR.GLbMed, NR.GLbUp, NR.GLbLo, NR.DDwAICmed, NR.DDwAICUp, NR.DDwAICLo, NR.q, NR.longev, NR.genL, NR.M, NR.alpha)
colnames(NR.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# GN
genL <- GN.gen.l; Nmed <- GN.n.md; Nmat1 <- GN.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- GN.pred.red.mn; Ninit <- GN.pop.found; Sb <- GN.b.lp; popmat <- GN.popmat.orig; mass <- GN.mass; prim <- GN.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
GN.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
GN.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
GN.DDw.dat <- data.frame(rep("GN",length(GN.DDwAICvec)), GN.DDwAICvec)
colnames(GN.DDw.dat) <- c("SP","DDwAIC")

GN.DDwAICmed <- median(GN.DDwAICvec, na.rm=T)
GN.DDwAICUp <- as.numeric(quantile(GN.DDwAICvec, probs=0.975, na.rm=T))
GN.DDwAICLo <- as.numeric(quantile(GN.DDwAICvec, probs=0.025, na.rm=T))

GN.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
GN.GLbMed <- -(median(GN.GLbs, na.rm=T))
GN.GLbLo <- -(as.numeric(quantile(GN.GLbs, probs=0.975, na.rm=T)))
GN.GLbUp <- -(as.numeric(quantile(GN.GLbs, probs=0.025, na.rm=T)))

GN.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
GN.RLbMed <- -(median(GN.RLbs, na.rm=T))
GN.RLbLo <- -(as.numeric(quantile(GN.RLbs, probs=0.975, na.rm=T)))
GN.RLbUp <- -(as.numeric(quantile(GN.RLbs, probs=0.025, na.rm=T)))

GN.predRedMed <- 1 - median(predRed, na.rm=T)
GN.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
GN.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
GN.predRedErr <- mean(c((GN.predRedMed-GN.predRedLo),(GN.predRedUp-GN.predRedMed)))

GN.q <- lenSeries - 1
GN.DDsStrength <- 1/(as.numeric(Sb/Ninit))
GN.longev <- round(dim(popmat)[2] - 1, 0)
GN.genL <- round(G.val(popmat, GN.longev), 2)
GN.M <- round(mass,0)
GN.alpha <- round(prim,0)
GN.out <- data.frame("GN",GN.predRedMed, GN.predRedErr, GN.GLbMed, GN.GLbUp, GN.GLbLo, GN.DDwAICmed, GN.DDwAICUp, GN.DDwAICLo, GN.q, GN.longev, GN.genL, GN.M, GN.alpha)
colnames(GN.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# DN
genL <- DN.gen.l; Nmed <- DN.n.md; Nmat1 <- DN.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- DN.pred.red.mn; Ninit <- DN.pop.found; Sb <- DN.b.lp; popmat <- DN.popmat.orig; mass <- DN.mass; prim <- DN.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
DN.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
DN.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
DN.DDw.dat <- data.frame(rep("DN",length(DN.DDwAICvec)), DN.DDwAICvec)
colnames(DN.DDw.dat) <- c("SP","DDwAIC")

DN.DDwAICmed <- median(DN.DDwAICvec, na.rm=T)
DN.DDwAICUp <- as.numeric(quantile(DN.DDwAICvec, probs=0.975, na.rm=T))
DN.DDwAICLo <- as.numeric(quantile(DN.DDwAICvec, probs=0.025, na.rm=T))

DN.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
DN.GLbMed <- -(median(DN.GLbs, na.rm=T))
DN.GLbLo <- -(as.numeric(quantile(DN.GLbs, probs=0.975, na.rm=T)))
DN.GLbUp <- -(as.numeric(quantile(DN.GLbs, probs=0.025, na.rm=T)))

DN.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
DN.RLbMed <- -(median(DN.RLbs, na.rm=T))
DN.RLbLo <- -(as.numeric(quantile(DN.RLbs, probs=0.975, na.rm=T)))
DN.RLbUp <- -(as.numeric(quantile(DN.RLbs, probs=0.025, na.rm=T)))

DN.predRedMed <- 1 - median(predRed, na.rm=T)
DN.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
DN.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
DN.predRedErr <- mean(c((DN.predRedMed-DN.predRedLo),(DN.predRedUp-DN.predRedMed)))

DN.q <- lenSeries - 1
DN.DDsStrength <- 1/(as.numeric(Sb/Ninit))
DN.longev <- round(dim(popmat)[2] - 1, 0)
DN.genL <- round(G.val(popmat, DN.longev), 2)
DN.M <- round(mass,0)
DN.alpha <- round(prim,0)
DN.out <- data.frame("DN",DN.predRedMed, DN.predRedErr, DN.GLbMed, DN.GLbUp, DN.GLbLo, DN.DDwAICmed, DN.DDwAICUp, DN.DDwAICLo, DN.q, DN.longev, DN.genL, DN.M, DN.alpha)
colnames(DN.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# AL
genL <- AL.gen.l; Nmed <- AL.n.md; Nmat1 <- AL.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- AL.pred.red.mn; Ninit <- AL.pop.found; Sb <- AL.b.lp; popmat <- AL.popmat.orig; mass <- AL.mass; prim <- AL.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
AL.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
AL.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
AL.DDw.dat <- data.frame(rep("AL",length(AL.DDwAICvec)), AL.DDwAICvec)
colnames(AL.DDw.dat) <- c("SP","DDwAIC")

AL.DDwAICmed <- median(AL.DDwAICvec, na.rm=T)
AL.DDwAICUp <- as.numeric(quantile(AL.DDwAICvec, probs=0.975, na.rm=T))
AL.DDwAICLo <- as.numeric(quantile(AL.DDwAICvec, probs=0.025, na.rm=T))

AL.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
AL.GLbMed <- -(median(AL.GLbs, na.rm=T))
AL.GLbLo <- -(as.numeric(quantile(AL.GLbs, probs=0.975, na.rm=T)))
AL.GLbUp <- -(as.numeric(quantile(AL.GLbs, probs=0.025, na.rm=T)))

AL.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
AL.RLbMed <- -(median(AL.RLbs, na.rm=T))
AL.RLbLo <- -(as.numeric(quantile(AL.RLbs, probs=0.975, na.rm=T)))
AL.RLbUp <- -(as.numeric(quantile(AL.RLbs, probs=0.025, na.rm=T)))

AL.predRedMed <- 1 - median(predRed, na.rm=T)
AL.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
AL.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
AL.predRedErr <- mean(c((AL.predRedMed-AL.predRedLo),(AL.predRedUp-AL.predRedMed)))

AL.q <- lenSeries - 1
AL.DDsStrength <- 1/(as.numeric(Sb/Ninit))
AL.longev <- round(dim(popmat)[2] - 1, 0)
AL.genL <- round(G.val(popmat, AL.longev), 2)
AL.M <- round(mass,0)
AL.alpha <- round(prim,0)
AL.out <- data.frame("AL",AL.predRedMed, AL.predRedErr, AL.GLbMed, AL.GLbUp, AL.GLbLo, AL.DDwAICmed, AL.DDwAICUp, AL.DDwAICLo, AL.q, AL.longev, AL.genL, AL.M, AL.alpha)
colnames(AL.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# TC
genL <- TC.gen.l; Nmed <- TC.n.md; Nmat1 <- TC.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- TC.pred.red.mn; Ninit <- TC.pop.found; Sb <- TC.b.lp; popmat <- TC.popmat.orig; mass <- TC.mass; prim <- TC.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
TC.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
TC.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
TC.DDw.dat <- data.frame(rep("TC",length(TC.DDwAICvec)), TC.DDwAICvec)
colnames(TC.DDw.dat) <- c("SP","DDwAIC")

TC.DDwAICmed <- median(TC.DDwAICvec, na.rm=T)
TC.DDwAICUp <- as.numeric(quantile(TC.DDwAICvec, probs=0.975, na.rm=T))
TC.DDwAICLo <- as.numeric(quantile(TC.DDwAICvec, probs=0.025, na.rm=T))

TC.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
TC.GLbMed <- -(median(TC.GLbs, na.rm=T))
TC.GLbLo <- -(as.numeric(quantile(TC.GLbs, probs=0.975, na.rm=T)))
TC.GLbUp <- -(as.numeric(quantile(TC.GLbs, probs=0.025, na.rm=T)))

TC.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
TC.RLbMed <- -(median(TC.RLbs, na.rm=T))
TC.RLbLo <- -(as.numeric(quantile(TC.RLbs, probs=0.975, na.rm=T)))
TC.RLbUp <- -(as.numeric(quantile(TC.RLbs, probs=0.025, na.rm=T)))

TC.predRedMed <- 1 - median(predRed, na.rm=T)
TC.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
TC.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
TC.predRedErr <- mean(c((TC.predRedMed-TC.predRedLo),(TC.predRedUp-TC.predRedMed)))

TC.q <- lenSeries - 1
TC.DDsStrength <- 1/(as.numeric(Sb/Ninit))
TC.longev <- round(dim(popmat)[2] - 1, 0)
TC.genL <- round(G.val(popmat, TC.longev), 2)
TC.M <- round(mass,0)
TC.alpha <- round(prim,0)
TC.out <- data.frame("TC",TC.predRedMed, TC.predRedErr, TC.GLbMed, TC.GLbUp, TC.GLbLo, TC.DDwAICmed, TC.DDwAICUp, TC.DDwAICLo, TC.q, TC.longev, TC.genL, TC.M, TC.alpha)
colnames(TC.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")


# TH
genL <- TH.gen.l; Nmed <- TH.n.md; Nmat1 <- TH.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- TH.pred.red.mn; Ninit <- TH.pop.found; Sb <- TH.b.lp; popmat <- TH.popmat.orig; mass <- TH.mass; prim <- TH.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
TH.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
TH.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
TH.DDw.dat <- data.frame(rep("TH",length(TH.DDwAICvec)), TH.DDwAICvec)
colnames(TH.DDw.dat) <- c("SP","DDwAIC")

TH.DDwAICmed <- median(TH.DDwAICvec, na.rm=T)
TH.DDwAICUp <- as.numeric(quantile(TH.DDwAICvec, probs=0.975, na.rm=T))
TH.DDwAICLo <- as.numeric(quantile(TH.DDwAICvec, probs=0.025, na.rm=T))

TH.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
TH.GLbMed <- -(median(TH.GLbs, na.rm=T))
TH.GLbLo <- -(as.numeric(quantile(TH.GLbs, probs=0.975, na.rm=T)))
TH.GLbUp <- -(as.numeric(quantile(TH.GLbs, probs=0.025, na.rm=T)))

TH.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
TH.RLbMed <- -(median(TH.RLbs, na.rm=T))
TH.RLbLo <- -(as.numeric(quantile(TH.RLbs, probs=0.975, na.rm=T)))
TH.RLbUp <- -(as.numeric(quantile(TH.RLbs, probs=0.025, na.rm=T)))

TH.predRedMed <- 1 - median(predRed, na.rm=T)
TH.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
TH.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
TH.predRedErr <- mean(c((TH.predRedMed-TH.predRedLo),(TH.predRedUp-TH.predRedMed)))

TH.q <- lenSeries - 1
TH.DDsStrength <- 1/(as.numeric(Sb/Ninit))
TH.longev <- round(dim(popmat)[2] - 1, 0)
TH.genL <- round(G.val(popmat, TH.longev), 2)
TH.M <- round(mass,0)
TH.alpha <- round(prim,0)
TH.out <- data.frame("TH",TH.predRedMed, TH.predRedErr, TH.GLbMed, TH.GLbUp, TH.GLbLo, TH.DDwAICmed, TH.DDwAICUp, TH.DDwAICLo, TH.q, TH.longev, TH.genL, TH.M, TH.alpha)
colnames(TH.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# SH
genL <- SH.gen.l; Nmed <- SH.n.md; Nmat1 <- SH.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- SH.pred.red.mn; Ninit <- SH.pop.found; Sb <- SH.b.lp; popmat <- SH.popmat.orig; mass <- SH.mass; prim <- SH.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
SH.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
SH.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
SH.DDw.dat <- data.frame(rep("SH",length(SH.DDwAICvec)), SH.DDwAICvec)
colnames(SH.DDw.dat) <- c("SP","DDwAIC")

SH.DDwAICmed <- median(SH.DDwAICvec, na.rm=T)
SH.DDwAICUp <- as.numeric(quantile(SH.DDwAICvec, probs=0.975, na.rm=T))
SH.DDwAICLo <- as.numeric(quantile(SH.DDwAICvec, probs=0.025, na.rm=T))

SH.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
SH.GLbMed <- -(median(SH.GLbs, na.rm=T))
SH.GLbLo <- -(as.numeric(quantile(SH.GLbs, probs=0.975, na.rm=T)))
SH.GLbUp <- -(as.numeric(quantile(SH.GLbs, probs=0.025, na.rm=T)))

SH.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
SH.RLbMed <- -(median(SH.RLbs, na.rm=T))
SH.RLbLo <- -(as.numeric(quantile(SH.RLbs, probs=0.975, na.rm=T)))
SH.RLbUp <- -(as.numeric(quantile(SH.RLbs, probs=0.025, na.rm=T)))

SH.predRedMed <- 1 - median(predRed, na.rm=T)
SH.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
SH.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
SH.predRedErr <- mean(c((SH.predRedMed-SH.predRedLo),(SH.predRedUp-SH.predRedMed)))

SH.q <- lenSeries - 1
SH.DDsStrength <- 1/(as.numeric(Sb/Ninit))
SH.longev <- round(dim(popmat)[2] - 1, 0)
SH.genL <- round(G.val(popmat, SH.longev), 2)
SH.M <- round(mass,0)
SH.alpha <- round(prim,0)
SH.out <- data.frame("SH",SH.predRedMed, SH.predRedErr, SH.GLbMed, SH.GLbUp, SH.GLbLo, SH.DDwAICmed, SH.DDwAICUp, SH.DDwAICLo, SH.q, SH.longev, SH.genL, SH.M, SH.alpha)
colnames(SH.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# DM
genL <- DM.gen.l; Nmed <- DM.n.md; Nmat1 <- DM.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- DM.pred.red.mn; Ninit <- DM.pop.found; Sb <- DM.b.lp; popmat <- DM.popmat.orig; mass <- DM.mass; prim <- DM.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
DM.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
DM.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
DM.DDw.dat <- data.frame(rep("DM",length(DM.DDwAICvec)), DM.DDwAICvec)
colnames(DM.DDw.dat) <- c("SP","DDwAIC")

DM.DDwAICmed <- median(DM.DDwAICvec, na.rm=T)
DM.DDwAICUp <- as.numeric(quantile(DM.DDwAICvec, probs=0.975, na.rm=T))
DM.DDwAICLo <- as.numeric(quantile(DM.DDwAICvec, probs=0.025, na.rm=T))

DM.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
DM.GLbMed <- -(median(DM.GLbs, na.rm=T))
DM.GLbLo <- -(as.numeric(quantile(DM.GLbs, probs=0.975, na.rm=T)))
DM.GLbUp <- -(as.numeric(quantile(DM.GLbs, probs=0.025, na.rm=T)))

DM.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
DM.RLbMed <- -(median(DM.RLbs, na.rm=T))
DM.RLbLo <- -(as.numeric(quantile(DM.RLbs, probs=0.975, na.rm=T)))
DM.RLbUp <- -(as.numeric(quantile(DM.RLbs, probs=0.025, na.rm=T)))

DM.predRedMed <- 1 - median(predRed, na.rm=T)
DM.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
DM.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
DM.predRedErr <- mean(c((DM.predRedMed-DM.predRedLo),(DM.predRedUp-DM.predRedMed)))

DM.q <- lenSeries - 1
DM.DDsStrength <- 1/(as.numeric(Sb/Ninit))
DM.longev <- round(dim(popmat)[2] - 1, 0)
DM.genL <- round(G.val(popmat, DM.longev), 2)
DM.M <- round(mass,0)
DM.alpha <- round(prim,0)
DM.out <- data.frame("DM",DM.predRedMed, DM.predRedErr, DM.GLbMed, DM.GLbUp, DM.GLbLo, DM.DDwAICmed, DM.DDwAICUp, DM.DDwAICLo, DM.q, DM.longev, DM.genL, DM.M, DM.alpha)
colnames(DM.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# TA
genL <- TA.gen.l; Nmed <- TA.n.md; Nmat1 <- TA.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- TA.pred.red.mn; Ninit <- TA.pop.found; Sb <- TA.b.lp; popmat <- TA.popmat.orig; mass <- TA.mass; prim <- TA.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
TA.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
TA.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
TA.DDw.dat <- data.frame(rep("TA",length(TA.DDwAICvec)), TA.DDwAICvec)
colnames(TA.DDw.dat) <- c("SP","DDwAIC")

TA.DDwAICmed <- median(TA.DDwAICvec, na.rm=T)
TA.DDwAICUp <- as.numeric(quantile(TA.DDwAICvec, probs=0.975, na.rm=T))
TA.DDwAICLo <- as.numeric(quantile(TA.DDwAICvec, probs=0.025, na.rm=T))

TA.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
TA.GLbMed <- -(median(TA.GLbs, na.rm=T))
TA.GLbLo <- -(as.numeric(quantile(TA.GLbs, probs=0.975, na.rm=T)))
TA.GLbUp <- -(as.numeric(quantile(TA.GLbs, probs=0.025, na.rm=T)))

TA.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
TA.RLbMed <- -(median(TA.RLbs, na.rm=T))
TA.RLbLo <- -(as.numeric(quantile(TA.RLbs, probs=0.975, na.rm=T)))
TA.RLbUp <- -(as.numeric(quantile(TA.RLbs, probs=0.025, na.rm=T)))

TA.predRedMed <- 1 - median(predRed, na.rm=T)
TA.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
TA.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
TA.predRedErr <- mean(c((TA.predRedMed-TA.predRedLo),(TA.predRedUp-TA.predRedMed)))

TA.q <- lenSeries - 1
TA.DDsStrength <- 1/(as.numeric(Sb/Ninit))
TA.longev <- round(dim(popmat)[2] - 1, 0)
TA.genL <- round(G.val(popmat, TA.longev), 2)
TA.M <- round(mass,0)
TA.alpha <- round(prim,0)
TA.out <- data.frame("TA",TA.predRedMed, TA.predRedErr, TA.GLbMed, TA.GLbUp, TA.GLbLo, TA.DDwAICmed, TA.DDwAICUp, TA.DDwAICLo, TA.q, TA.longev, TA.genL, TA.M, TA.alpha)
colnames(TA.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

# MR
genL <- MR.gen.l; Nmed <- MR.n.md; Nmat1 <- MR.n.sums.matb; zeroFix <- which(rowSums(Nmat1)==0)
if (length(zeroFix)==0) {
  Nmat2 <- Nmat1
} else {
  Nmat2 <- Nmat1[-zeroFix,]
}
lenvec <- apply(Nmat2, MARGIN=1, count_if, criterion=gt(0))
lenvecthree <- which(lenvec < 3)
if (length(lenvecthree)==0) {
  Nmat <- Nmat2
} else {
  Nmat <- Nmat2[-lenvecthree,]
}
predRed <- MR.pred.red.mn; Ninit <- MR.pop.found; Sb <- MR.b.lp; popmat <- MR.popmat.orig; mass <- MR.mass; prim <- MR.alpha
lenSeries <- dim(Nmat)[2]

RWout <- unlist(apply(Nmat, MARGIN=1, RWfunc))
EXout <- unlist(apply(Nmat, MARGIN=1, EXfunc))
RLout <- unlist(apply(Nmat, MARGIN=1, RLfunc))
GLout <- unlist(apply(Nmat, MARGIN=1, GLfunc))
RWAICs <- as.numeric(RWout[which(attr(RWout, "names") == "AICc")])
EXAICs <- as.numeric(EXout[which(attr(EXout, "names") == "AICc")])
RLAICs <- as.numeric(RLout[which(attr(RLout, "names") == "AICc")])
GLAICs <- as.numeric(GLout[which(attr(GLout, "names") == "AICc")])
AICtable <- data.frame(RWAICs,EXAICs,RLAICs,GLAICs)
dAICtable <- as.data.frame(t(apply(AICtable, MARGIN=1, delta.IC)))
wAICtable <- as.data.frame(t(apply(dAICtable, MARGIN=1, weight.IC)))
MR.DDwAICvec <- apply(wAICtable[,3:4], MARGIN=1, sum)  
MR.DIwAICvec <- apply(wAICtable[,1:2], MARGIN=1, sum)  
MR.DDw.dat <- data.frame(rep("MR",length(MR.DDwAICvec)), MR.DDwAICvec)
colnames(MR.DDw.dat) <- c("SP","DDwAIC")

MR.DDwAICmed <- median(MR.DDwAICvec, na.rm=T)
MR.DDwAICUp <- as.numeric(quantile(MR.DDwAICvec, probs=0.975, na.rm=T))
MR.DDwAICLo <- as.numeric(quantile(MR.DDwAICvec, probs=0.025, na.rm=T))

MR.GLbs <- as.numeric(GLout[which(attr(GLout, "names") == "b")])
MR.GLbMed <- -(median(MR.GLbs, na.rm=T))
MR.GLbLo <- -(as.numeric(quantile(MR.GLbs, probs=0.975, na.rm=T)))
MR.GLbUp <- -(as.numeric(quantile(MR.GLbs, probs=0.025, na.rm=T)))

MR.RLbs <- as.numeric(RLout[which(attr(RLout, "names") == "b")])
MR.RLbMed <- -(median(MR.RLbs, na.rm=T))
MR.RLbLo <- -(as.numeric(quantile(MR.RLbs, probs=0.975, na.rm=T)))
MR.RLbUp <- -(as.numeric(quantile(MR.RLbs, probs=0.025, na.rm=T)))

MR.predRedMed <- 1 - median(predRed, na.rm=T)
MR.predRedUp <- 1 - quantile(predRed, probs=0.025, na.rm=T)
MR.predRedLo <- 1 - quantile(predRed, probs=0.975, na.rm=T)
MR.predRedErr <- mean(c((MR.predRedMed-MR.predRedLo),(MR.predRedUp-MR.predRedMed)))

MR.q <- lenSeries - 1
MR.DDsStrength <- 1/(as.numeric(Sb/Ninit))
MR.longev <- round(dim(popmat)[2] - 1, 0)
MR.genL <- round(G.val(popmat, MR.longev), 2)
MR.M <- round(mass,0)
MR.alpha <- round(prim,0)
MR.out <- data.frame("MR",MR.predRedMed, MR.predRedErr, MR.GLbMed, MR.GLbUp, MR.GLbLo, MR.DDwAICmed, MR.DDwAICUp, MR.DDwAICLo, MR.q, MR.longev, MR.genL, MR.M, MR.alpha)
colnames(MR.out) <- c("SP","predRedMd","predRedE","GLbMd","GLbUp","GLbLo","DDwAICmd","DDwAICup","DDwAIClo","q","longev","G","M","alpha")

## create output datafile for species-level summary data
summ.dat <- rbind(DP.out, PA.out, ZT.out, PH.out, VU.out, PG.out, SS.out, PT.out,
                  SO.out, MN.out, OR.out, NR.out, GN.out, DN.out, AL.out, TC.out,
                  TH.out, SH.out, DM.out, TA.out, MR.out)

## build histograms themselves
# combine data frames
DDw.dat <- rbind(DP.DDw.dat, PA.DDw.dat, ZT.DDw.dat, PH.DDw.dat, VU.DDw.dat, PG.DDw.dat, SS.DDw.dat, PT.DDw.dat,
                 SO.DDw.dat, MN.DDw.dat, OR.DDw.dat, NR.DDw.dat, GN.DDw.dat, DN.DDw.dat, AL.DDw.dat, TC.DDw.dat,
                 TH.DDw.dat, SH.DDw.dat, DM.DDw.dat, TA.DDw.dat, MR.DDw.dat)
DDw.dat$SP <- as.factor(DDw.dat$SP)
nspp <- length(table(DDw.dat$SP))

mu <- DDw.dat %>% 
  group_by(SP) %>%
  summarise(grp.mean = mean(DDwAIC))
mu

theme_set(theme_ridges())
ggplot(DDw.dat, aes(x = DDwAIC, y = SP, show.legend=F)) +
  xlim(0.95, 1.0) +
  xlab("wAICc") + ylab("") +
  geom_density_ridges(aes(fill = SP), alpha=0.6, show.legend = FALSE) +
  scale_fill_manual(values = rep("blue", nspp))

ggdensity(DDw.dat, x = "DDwAIC", y="..ndensity..", add = "none", rug = TRUE, color=NA, fill = "SP", alpha=0.3) +
  xlim(0, 1.0) +
  xlab("density-feedback wAICc") + ylab("") +
  scale_fill_manual(values = rep("blue", nspp)) +
  theme_bw() +
  theme(legend.position="none")

