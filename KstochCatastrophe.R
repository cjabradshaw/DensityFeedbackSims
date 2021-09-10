#############################################################################################################################
## density-feedback simulations
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

#########################################################
## CARRYING CAPACITY STABLE, STOCHASTICALLY RESAMPLED, ##
## VARIANCE CONSTANT                                   ##
#########################################################

## DIPROTODON (DP)
load("DP.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*DP.gen.l, 0) - 1)
DP.b.lp.vec2 <- rnorm(n=t, mean=DP.b.lp, sd=0.05*DP.b.lp)
DP.b.lp.vec <- ifelse(DP.b.lp.vec2 < 10, 10, DP.b.lp.vec2)
plot(DP.b.lp.vec, type="l", lty=2)
mean(log(DP.b.lp.vec[2:length(DP.b.lp.vec)] / DP.b.lp.vec[1:(length(DP.b.lp.vec)-1)]))

## set storage matrices & vectors
DP.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DP.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  DP.popmat <- DP.popmat.orig
  
  DP.n.mat <- matrix(0, nrow=DP.age.max+1,ncol=(t+1))
  DP.n.mat[,1] <- DP.init.vec
  DP.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    DP.s.alpha <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$alpha
    DP.s.beta <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$beta
    DP.s.stoch <- rbeta(length(DP.s.alpha), DP.s.alpha, DP.s.beta)
    
    if (rbinom(1, 1, 0.14/DP.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      DP.s.stoch <- DP.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    DP.fert.stch <- rnorm(length(DP.popmat[,1]), DP.pred.p.mm, DP.m.sd.vec)
    
    DP.totN.i <- sum(DP.n.mat[,i], na.rm=T)
    DP.pred.red.vec[i] <- DP.a.lp/(1+(DP.totN.i/DP.b.lp.vec[i])^DP.c.lp)
    
    diag(DP.popmat[2:(DP.age.max+1),]) <- (DP.s.stoch[-(DP.age.max+1)])*DP.pred.red.vec[i]
    DP.popmat[DP.age.max+1,DP.age.max+1] <- (DP.s.stoch[DP.age.max+1])*DP.pred.red.vec[i]
    DP.popmat[1,] <- ifelse(DP.fert.stch < 0, 0, DP.fert.stch)
    DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]
    
  } # end i loop
  
  DP.n.sums.mat[e,] <- (as.vector(colSums(DP.n.mat)))
  DP.pred.red.mn[e] <- mean(DP.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
DP.n.sums.matb <- DP.n.sums.mat[, -(1:round(DP.gen.l*2, 0))]

# total N
DP.n.md <- apply(DP.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
DP.n.up <- apply(DP.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DP.n.lo <- apply(DP.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(DP.n.sums.matb)[2])
plot(yrs, DP.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(DP.n.lo), max(DP.n.up)))
lines(yrs,DP.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,DP.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(DP.n.md[2:length(yrs)] / DP.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("DPKstoch.RData")
rm(list = ls())



## PALORCHESTES (PA)
load("PA.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*PA.gen.l, 0) - 1)
PA.b.lp.vec2 <- rnorm(n=t, mean=PA.b.lp, sd=0.05*PA.b.lp)
PA.b.lp.vec <- ifelse(PA.b.lp.vec2 < 10, 10, PA.b.lp.vec2)
plot(PA.b.lp.vec, type="l", lty=2)
mean(log(PA.b.lp.vec[2:length(PA.b.lp.vec)] / PA.b.lp.vec[1:(length(PA.b.lp.vec)-1)]))

## set storage matrices & vectors
PA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PA.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  PA.popmat <- PA.popmat.orig
  
  PA.n.mat <- matrix(0, nrow=PA.age.max+1,ncol=(t+1))
  PA.n.mat[,1] <- PA.init.vec
  PA.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    PA.s.alpha <- estBetaParams(PA.Sx, PA.s.sd.vec^2)$alpha
    PA.s.beta <- estBetaParams(PA.Sx, PA.s.sd.vec^2)$beta
    PA.s.stoch <- rbeta(length(PA.s.alpha), PA.s.alpha, PA.s.beta)
    
    if (rbinom(1, 1, 0.14/PA.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      PA.s.stoch <- PA.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    PA.fert.stch <- rnorm(length(PA.popmat[,1]), PA.pred.p.mm, PA.m.sd.vec)
    
    PA.totN.i <- sum(PA.n.mat[,i], na.rm=T)
    PA.pred.red.vec[i] <- PA.a.lp/(1+(PA.totN.i/PA.b.lp.vec[i])^PA.c.lp)
    
    diag(PA.popmat[2:(PA.age.max+1),]) <- (PA.s.stoch[-(PA.age.max+1)])*PA.pred.red.vec[i]
    PA.popmat[PA.age.max+1,PA.age.max+1] <- (PA.s.stoch[PA.age.max+1])*PA.pred.red.vec[i]
    PA.popmat[1,] <- ifelse(PA.fert.stch < 0, 0, PA.fert.stch)
    PA.n.mat[,i+1] <- PA.popmat %*% PA.n.mat[,i]
    
  } # end i loop
  
  PA.n.sums.mat[e,] <- (as.vector(colSums(PA.n.mat)))
  PA.pred.red.mn[e] <- mean(PA.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
PA.n.sums.matb <- PA.n.sums.mat[, -(1:round(PA.gen.l*2, 0))]

# total N
PA.n.md <- apply(PA.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
PA.n.up <- apply(PA.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PA.n.lo <- apply(PA.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(PA.n.sums.matb)[2])
plot(yrs, PA.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(PA.n.lo), max(PA.n.up)))
lines(yrs,PA.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,PA.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(PA.n.md[2:length(yrs)] / PA.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("PAKstoch.RData")
rm(list = ls())



## ZYGOMATURUS (ZT)
load("ZT.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*ZT.gen.l, 0) - 1)
ZT.b.lp.vec2 <- rnorm(n=t, mean=ZT.b.lp, sd=0.05*ZT.b.lp)
ZT.b.lp.vec <- ifelse(ZT.b.lp.vec2 < 10, 10, ZT.b.lp.vec2)
plot(ZT.b.lp.vec, type="l", lty=2)
mean(log(ZT.b.lp.vec[2:length(ZT.b.lp.vec)] / ZT.b.lp.vec[1:(length(ZT.b.lp.vec)-1)]))

## set storage matrices & vectors
ZT.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
ZT.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  ZT.popmat <- ZT.popmat.orig
  
  ZT.n.mat <- matrix(0, nrow=ZT.age.max+1,ncol=(t+1))
  ZT.n.mat[,1] <- ZT.init.vec
  ZT.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    ZT.s.alpha <- estBetaParams(ZT.Sx, ZT.s.sd.vec^2)$alpha
    ZT.s.beta <- estBetaParams(ZT.Sx, ZT.s.sd.vec^2)$beta
    ZT.s.stoch <- rbeta(length(ZT.s.alpha), ZT.s.alpha, ZT.s.beta)
    
    if (rbinom(1, 1, 0.14/ZT.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      ZT.s.stoch <- ZT.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    ZT.fert.stch <- rnorm(length(ZT.popmat[,1]), ZT.pred.p.mm, ZT.m.sd.vec)
    
    ZT.totN.i <- sum(ZT.n.mat[,i], na.rm=T)
    ZT.pred.red.vec[i] <- ZT.a.lp/(1+(ZT.totN.i/ZT.b.lp.vec[i])^ZT.c.lp)
    
    diag(ZT.popmat[2:(ZT.age.max+1),]) <- (ZT.s.stoch[-(ZT.age.max+1)])*ZT.pred.red.vec[i]
    ZT.popmat[ZT.age.max+1,ZT.age.max+1] <- (ZT.s.stoch[ZT.age.max+1])*ZT.pred.red.vec[i]
    ZT.popmat[1,] <- ifelse(ZT.fert.stch < 0, 0, ZT.fert.stch)
    ZT.n.mat[,i+1] <- ZT.popmat %*% ZT.n.mat[,i]
    
  } # end i loop
  
  ZT.n.sums.mat[e,] <- (as.vector(colSums(ZT.n.mat)))
  ZT.pred.red.mn[e] <- mean(ZT.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
ZT.n.sums.matb <- ZT.n.sums.mat[, -(1:round(ZT.gen.l*2, 0))]

# total N
ZT.n.md <- apply(ZT.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
ZT.n.up <- apply(ZT.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
ZT.n.lo <- apply(ZT.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(ZT.n.sums.matb)[2])
plot(yrs, ZT.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(ZT.n.lo), max(ZT.n.up)))
lines(yrs,ZT.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,ZT.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(ZT.n.md[2:length(yrs)] / ZT.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("ZTKstoch.RData")
rm(list = ls())



## PHASCOLONUS (PH)
load("PH.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*PH.gen.l, 0) - 1)
PH.b.lp.vec2 <- rnorm(n=t, mean=PH.b.lp, sd=0.05*PH.b.lp)
PH.b.lp.vec <- ifelse(PH.b.lp.vec2 < 10, 10, PH.b.lp.vec2)
plot(PH.b.lp.vec, type="l", lty=2)
mean(log(PH.b.lp.vec[2:length(PH.b.lp.vec)] / PH.b.lp.vec[1:(length(PH.b.lp.vec)-1)]))

## set storage matrices & vectors
PH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PH.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  PH.popmat <- PH.popmat.orig
  
  PH.n.mat <- matrix(0, nrow=PH.age.max+1,ncol=(t+1))
  PH.n.mat[,1] <- PH.init.vec
  PH.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    PH.s.alpha <- estBetaParams(PH.Sx, PH.s.sd.vec^2)$alpha
    PH.s.beta <- estBetaParams(PH.Sx, PH.s.sd.vec^2)$beta
    PH.s.stoch <- rbeta(length(PH.s.alpha), PH.s.alpha, PH.s.beta)
    
    if (rbinom(1, 1, 0.14/PH.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      PH.s.stoch <- PH.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    PH.fert.stch <- rnorm(length(PH.popmat[,1]), PH.pred.p.mm, PH.m.sd.vec)
    
    PH.totN.i <- sum(PH.n.mat[,i], na.rm=T)
    PH.pred.red.vec[i] <- PH.a.lp/(1+(PH.totN.i/PH.b.lp.vec[i])^PH.c.lp)
    
    diag(PH.popmat[2:(PH.age.max+1),]) <- (PH.s.stoch[-(PH.age.max+1)])*PH.pred.red.vec[i]
    PH.popmat[PH.age.max+1,PH.age.max+1] <- (PH.s.stoch[PH.age.max+1])*PH.pred.red.vec[i]
    PH.popmat[1,] <- ifelse(PH.fert.stch < 0, 0, PH.fert.stch)
    PH.n.mat[,i+1] <- PH.popmat %*% PH.n.mat[,i]
    
  } # end i loop
  
  PH.n.sums.mat[e,] <- (as.vector(colSums(PH.n.mat)))
  PH.pred.red.mn[e] <- mean(PH.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
PH.n.sums.matb <- PH.n.sums.mat[, -(1:round(PH.gen.l*2, 0))]

# total N
PH.n.md <- apply(PH.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
PH.n.up <- apply(PH.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PH.n.lo <- apply(PH.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(PH.n.sums.matb)[2])
plot(yrs, PH.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(PH.n.lo), max(PH.n.up)))
lines(yrs,PH.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,PH.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(PH.n.md[2:length(yrs)] / PH.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("PHKstoch.RData")
rm(list = ls())



## VOMBATUS (VU)
load("VU.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*VU.gen.l, 0) - 1)
VU.b.lp.vec2 <- rnorm(n=t, mean=VU.b.lp, sd=0.05*VU.b.lp)
VU.b.lp.vec <- ifelse(VU.b.lp.vec2 < 10, 10, VU.b.lp.vec2)
plot(VU.b.lp.vec, type="l", lty=2)
mean(log(VU.b.lp.vec[2:length(VU.b.lp.vec)] / VU.b.lp.vec[1:(length(VU.b.lp.vec)-1)]))

## set storage matrices & vectors
VU.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
VU.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  VU.popmat <- VU.popmat.orig
  
  VU.n.mat <- matrix(0, nrow=VU.age.max+1,ncol=(t+1))
  VU.n.mat[,1] <- VU.init.vec
  VU.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    VU.s.alpha <- estBetaParams(VU.Sx, VU.s.sd.vec^2)$alpha
    VU.s.beta <- estBetaParams(VU.Sx, VU.s.sd.vec^2)$beta
    VU.s.stoch <- rbeta(length(VU.s.alpha), VU.s.alpha, VU.s.beta)
    
    if (rbinom(1, 1, 0.14/VU.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      VU.s.stoch <- VU.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    VU.fert.stch <- rnorm(length(VU.popmat[,1]), VU.pred.p.mm, VU.m.sd.vec)
    
    VU.totN.i <- sum(VU.n.mat[,i], na.rm=T)
    VU.pred.red.vec[i] <- VU.a.lp/(1+(VU.totN.i/VU.b.lp.vec[i])^VU.c.lp)
    
    diag(VU.popmat[2:(VU.age.max+1),]) <- (VU.s.stoch[-(VU.age.max+1)])*VU.pred.red.vec[i]
    VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
    VU.popmat[1,] <- ifelse(VU.fert.stch < 0, 0, VU.fert.stch)
    VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]
    
  } # end i loop
  
  VU.n.sums.mat[e,] <- (as.vector(colSums(VU.n.mat)))
  VU.pred.red.mn[e] <- mean(VU.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
VU.n.sums.matb <- VU.n.sums.mat[, -(1:round(VU.gen.l*2, 0))]

# total N
VU.n.md <- apply(VU.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
VU.n.up <- apply(VU.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
VU.n.lo <- apply(VU.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(VU.n.sums.matb)[2])
plot(yrs, VU.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(VU.n.lo), max(VU.n.up)))
lines(yrs,VU.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,VU.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(VU.n.md[2:length(yrs)] / VU.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("VUKstoch.RData")
rm(list = ls())



## PROCOPTODON (PG)
load("PG.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*PG.gen.l, 0) - 1)
PG.b.lp.vec2 <- rnorm(n=t, mean=PG.b.lp, sd=0.05*PG.b.lp)
PG.b.lp.vec <- ifelse(PG.b.lp.vec2 < 10, 10, PG.b.lp.vec2)
plot(PG.b.lp.vec, type="l", lty=2)
mean(log(PG.b.lp.vec[2:length(PG.b.lp.vec)] / PG.b.lp.vec[1:(length(PG.b.lp.vec)-1)]))

## set storage matrices & vectors
PG.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PG.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  PG.popmat <- PG.popmat.orig
  
  PG.n.mat <- matrix(0, nrow=PG.age.max+1,ncol=(t+1))
  PG.n.mat[,1] <- PG.init.vec
  PG.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    PG.s.alpha <- estBetaParams(PG.Sx, PG.s.sd.vec^2)$alpha
    PG.s.beta <- estBetaParams(PG.Sx, PG.s.sd.vec^2)$beta
    PG.s.stoch <- rbeta(length(PG.s.alpha), PG.s.alpha, PG.s.beta)
    
    if (rbinom(1, 1, 0.14/PG.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      PG.s.stoch <- PG.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    PG.fert.stch <- rnorm(length(PG.popmat[,1]), PG.pred.p.mm, PG.m.sd.vec)
    
    PG.totN.i <- sum(PG.n.mat[,i], na.rm=T)
    PG.pred.red.vec[i] <- PG.a.lp/(1+(PG.totN.i/PG.b.lp.vec[i])^PG.c.lp)
    
    diag(PG.popmat[2:(PG.age.max+1),]) <- (PG.s.stoch[-(PG.age.max+1)])*PG.pred.red.vec[i]
    PG.popmat[PG.age.max+1,PG.age.max+1] <- (PG.s.stoch[PG.age.max+1])*PG.pred.red.vec[i]
    PG.popmat[1,] <- ifelse(PG.fert.stch < 0, 0, PG.fert.stch)
    PG.n.mat[,i+1] <- PG.popmat %*% PG.n.mat[,i]
    
  } # end i loop
  
  PG.n.sums.mat[e,] <- (as.vector(colSums(PG.n.mat)))
  PG.pred.red.mn[e] <- mean(PG.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
PG.n.sums.matb <- PG.n.sums.mat[, -(1:round(PG.gen.l*2, 0))]

# total N
PG.n.md <- apply(PG.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
PG.n.up <- apply(PG.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PG.n.lo <- apply(PG.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(PG.n.sums.matb)[2])
plot(yrs, PG.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(PG.n.lo), max(PG.n.up)))
lines(yrs,PG.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,PG.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(PG.n.md[2:length(yrs)] / PG.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("PGKstoch.RData")
rm(list = ls())



## STHENURUS (SS)
load("SS.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*SS.gen.l, 0) - 1)
SS.b.lp.vec2 <- rnorm(n=t, mean=SS.b.lp, sd=0.05*SS.b.lp)
SS.b.lp.vec <- ifelse(SS.b.lp.vec2 < 10, 10, SS.b.lp.vec2)
plot(SS.b.lp.vec, type="l", lty=2)
mean(log(SS.b.lp.vec[2:length(SS.b.lp.vec)] / SS.b.lp.vec[1:(length(SS.b.lp.vec)-1)]))

## set storage matrices & vectors
SS.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
SS.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  SS.popmat <- SS.popmat.orig
  
  SS.n.mat <- matrix(0, nrow=SS.age.max+1,ncol=(t+1))
  SS.n.mat[,1] <- SS.init.vec
  SS.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    SS.s.alpha <- estBetaParams(SS.Sx, SS.s.sd.vec^2)$alpha
    SS.s.beta <- estBetaParams(SS.Sx, SS.s.sd.vec^2)$beta
    SS.s.stoch <- rbeta(length(SS.s.alpha), SS.s.alpha, SS.s.beta)
    
    if (rbinom(1, 1, 0.14/SS.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      SS.s.stoch <- SS.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    SS.fert.stch <- rnorm(length(SS.popmat[,1]), SS.pred.p.mm, SS.m.sd.vec)
    
    SS.totN.i <- sum(SS.n.mat[,i], na.rm=T)
    SS.pred.red.vec[i] <- SS.a.lp/(1+(SS.totN.i/SS.b.lp.vec[i])^SS.c.lp)
    
    diag(SS.popmat[2:(SS.age.max+1),]) <- (SS.s.stoch[-(SS.age.max+1)])*SS.pred.red.vec[i]
    SS.popmat[SS.age.max+1,SS.age.max+1] <- (SS.s.stoch[SS.age.max+1])*SS.pred.red.vec[i]
    SS.popmat[1,] <- ifelse(SS.fert.stch < 0, 0, SS.fert.stch)
    SS.n.mat[,i+1] <- SS.popmat %*% SS.n.mat[,i]
    
  } # end i loop
  
  SS.n.sums.mat[e,] <- (as.vector(colSums(SS.n.mat)))
  SS.pred.red.mn[e] <- mean(SS.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
SS.n.sums.matb <- SS.n.sums.mat[, -(1:round(SS.gen.l*2, 0))]

# total N
SS.n.md <- apply(SS.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
SS.n.up <- apply(SS.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
SS.n.lo <- apply(SS.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(SS.n.sums.matb)[2])
plot(yrs, SS.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(SS.n.lo), max(SS.n.up)))
lines(yrs,SS.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,SS.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(SS.n.md[2:length(yrs)] / SS.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("SSKstoch.RData")
rm(list = ls())



## PROTEMNODON (PT)
load("PT.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*PT.gen.l, 0) - 1)
PT.b.lp.vec2 <- rnorm(n=t, mean=PT.b.lp, sd=0.05*PT.b.lp)
PT.b.lp.vec <- ifelse(PT.b.lp.vec2 < 10, 10, PT.b.lp.vec2)
plot(PT.b.lp.vec, type="l", lty=2)
mean(log(PT.b.lp.vec[2:length(PT.b.lp.vec)] / PT.b.lp.vec[1:(length(PT.b.lp.vec)-1)]))

## set storage matrices & vectors
PT.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PT.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  PT.popmat <- PT.popmat.orig
  
  PT.n.mat <- matrix(0, nrow=PT.age.max+1,ncol=(t+1))
  PT.n.mat[,1] <- PT.init.vec
  PT.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    PT.s.alpha <- estBetaParams(PT.Sx, PT.s.sd.vec^2)$alpha
    PT.s.beta <- estBetaParams(PT.Sx, PT.s.sd.vec^2)$beta
    PT.s.stoch <- rbeta(length(PT.s.alpha), PT.s.alpha, PT.s.beta)
    
    if (rbinom(1, 1, 0.14/PT.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      PT.s.stoch <- PT.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    PT.fert.stch <- rnorm(length(PT.popmat[,1]), PT.pred.p.mm, PT.m.sd.vec)
    
    PT.totN.i <- sum(PT.n.mat[,i], na.rm=T)
    PT.pred.red.vec[i] <- PT.a.lp/(1+(PT.totN.i/PT.b.lp.vec[i])^PT.c.lp)
    
    diag(PT.popmat[2:(PT.age.max+1),]) <- (PT.s.stoch[-(PT.age.max+1)])*PT.pred.red.vec[i]
    PT.popmat[PT.age.max+1,PT.age.max+1] <- (PT.s.stoch[PT.age.max+1])*PT.pred.red.vec[i]
    PT.popmat[1,] <- ifelse(PT.fert.stch < 0, 0, PT.fert.stch)
    PT.n.mat[,i+1] <- PT.popmat %*% PT.n.mat[,i]
    
  } # end i loop
  
  PT.n.sums.mat[e,] <- (as.vector(colSums(PT.n.mat)))
  PT.pred.red.mn[e] <- mean(PT.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
PT.n.sums.matb <- PT.n.sums.mat[, -(1:round(PT.gen.l*2, 0))]

# total N
PT.n.md <- apply(PT.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
PT.n.up <- apply(PT.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PT.n.lo <- apply(PT.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(PT.n.sums.matb)[2])
plot(yrs, PT.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(PT.n.lo), max(PT.n.up)))
lines(yrs,PT.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,PT.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(PT.n.md[2:length(yrs)] / PT.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("PTKstoch.RData")
rm(list = ls())



## SIMOSTHENURUS (SO)
load("SO.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*SO.gen.l, 0) - 1)
SO.b.lp.vec2 <- rnorm(n=t, mean=SO.b.lp, sd=0.05*SO.b.lp)
SO.b.lp.vec <- ifelse(SO.b.lp.vec2 < 10, 10, SO.b.lp.vec2)
plot(SO.b.lp.vec, type="l", lty=2)
mean(log(SO.b.lp.vec[2:length(SO.b.lp.vec)] / SO.b.lp.vec[1:(length(SO.b.lp.vec)-1)]))

## set storage matrices & vectors
SO.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
SO.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  SO.popmat <- SO.popmat.orig
  
  SO.n.mat <- matrix(0, nrow=SO.age.max+1,ncol=(t+1))
  SO.n.mat[,1] <- SO.init.vec
  SO.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    SO.s.alpha <- estBetaParams(SO.Sx, SO.s.sd.vec^2)$alpha
    SO.s.beta <- estBetaParams(SO.Sx, SO.s.sd.vec^2)$beta
    SO.s.stoch <- rbeta(length(SO.s.alpha), SO.s.alpha, SO.s.beta)
    
    if (rbinom(1, 1, 0.14/SO.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      SO.s.stoch <- SO.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    SO.fert.stch <- rnorm(length(SO.popmat[,1]), SO.pred.p.mm, SO.m.sd.vec)
    
    SO.totN.i <- sum(SO.n.mat[,i], na.rm=T)
    SO.pred.red.vec[i] <- SO.a.lp/(1+(SO.totN.i/SO.b.lp.vec[i])^SO.c.lp)
    
    diag(SO.popmat[2:(SO.age.max+1),]) <- (SO.s.stoch[-(SO.age.max+1)])*SO.pred.red.vec[i]
    SO.popmat[SO.age.max+1,SO.age.max+1] <- (SO.s.stoch[SO.age.max+1])*SO.pred.red.vec[i]
    SO.popmat[1,] <- ifelse(SO.fert.stch < 0, 0, SO.fert.stch)
    SO.n.mat[,i+1] <- SO.popmat %*% SO.n.mat[,i]
    
  } # end i loop
  
  SO.n.sums.mat[e,] <- (as.vector(colSums(SO.n.mat)))
  SO.pred.red.mn[e] <- mean(SO.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
SO.n.sums.matb <- SO.n.sums.mat[, -(1:round(SO.gen.l*2, 0))]

# total N
SO.n.md <- apply(SO.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
SO.n.up <- apply(SO.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
SO.n.lo <- apply(SO.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(SO.n.sums.matb)[2])
plot(yrs, SO.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(SO.n.lo), max(SO.n.up)))
lines(yrs,SO.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,SO.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(SO.n.md[2:length(yrs)] / SO.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("SOKstoch.RData")
rm(list = ls())



## METASTHENURUS (MN)
load("MN.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*MN.gen.l, 0) - 1)
MN.b.lp.vec2 <- rnorm(n=t, mean=MN.b.lp, sd=0.05*MN.b.lp)
MN.b.lp.vec <- ifelse(MN.b.lp.vec2 < 10, 10, MN.b.lp.vec2)
plot(MN.b.lp.vec, type="l", lty=2)
mean(log(MN.b.lp.vec[2:length(MN.b.lp.vec)] / MN.b.lp.vec[1:(length(MN.b.lp.vec)-1)]))

## set storage matrices & vectors
MN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
MN.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  MN.popmat <- MN.popmat.orig
  
  MN.n.mat <- matrix(0, nrow=MN.age.max+1,ncol=(t+1))
  MN.n.mat[,1] <- MN.init.vec
  MN.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    MN.s.alpha <- estBetaParams(MN.Sx, MN.s.sd.vec^2)$alpha
    MN.s.beta <- estBetaParams(MN.Sx, MN.s.sd.vec^2)$beta
    MN.s.stoch <- rbeta(length(MN.s.alpha), MN.s.alpha, MN.s.beta)
    
    if (rbinom(1, 1, 0.14/MN.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      MN.s.stoch <- MN.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    MN.fert.stch <- rnorm(length(MN.popmat[,1]), MN.pred.p.mm, MN.m.sd.vec)
    
    MN.totN.i <- sum(MN.n.mat[,i], na.rm=T)
    MN.pred.red.vec[i] <- MN.a.lp/(1+(MN.totN.i/MN.b.lp.vec[i])^MN.c.lp)
    
    diag(MN.popmat[2:(MN.age.max+1),]) <- (MN.s.stoch[-(MN.age.max+1)])*MN.pred.red.vec[i]
    MN.popmat[MN.age.max+1,MN.age.max+1] <- (MN.s.stoch[MN.age.max+1])*MN.pred.red.vec[i]
    MN.popmat[1,] <- ifelse(MN.fert.stch < 0, 0, MN.fert.stch)
    MN.n.mat[,i+1] <- MN.popmat %*% MN.n.mat[,i]
    
  } # end i loop
  
  MN.n.sums.mat[e,] <- (as.vector(colSums(MN.n.mat)))
  MN.pred.red.mn[e] <- mean(MN.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
MN.n.sums.matb <- MN.n.sums.mat[, -(1:round(MN.gen.l*2, 0))]

# total N
MN.n.md <- apply(MN.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
MN.n.up <- apply(MN.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
MN.n.lo <- apply(MN.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(MN.n.sums.matb)[2])
plot(yrs, MN.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(MN.n.lo), max(MN.n.up)))
lines(yrs,MN.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,MN.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(MN.n.md[2:length(yrs)] / MN.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("MNKstoch.RData")
rm(list = ls())



## OSPHRANTER (OR)
load("OR.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*OR.gen.l, 0) - 1)
OR.b.lp.vec2 <- rnorm(n=t, mean=OR.b.lp, sd=0.05*OR.b.lp)
OR.b.lp.vec <- ifelse(OR.b.lp.vec2 < 10, 10, OR.b.lp.vec2)
plot(OR.b.lp.vec, type="l", lty=2)
mean(log(OR.b.lp.vec[2:length(OR.b.lp.vec)] / OR.b.lp.vec[1:(length(OR.b.lp.vec)-1)]))

## set storage matrices & vectors
OR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
OR.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  OR.popmat <- OR.popmat.orig
  
  OR.n.mat <- matrix(0, nrow=OR.age.max+1,ncol=(t+1))
  OR.n.mat[,1] <- OR.init.vec
  OR.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    OR.s.alpha <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$alpha
    OR.s.beta <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$beta
    OR.s.stoch <- rbeta(length(OR.s.alpha), OR.s.alpha, OR.s.beta)
    
    if (rbinom(1, 1, 0.14/OR.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      OR.s.stoch <- OR.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    OR.fert.stch <- rnorm(length(OR.popmat[,1]), OR.pred.p.mm, OR.m.sd.vec)
    
    OR.totN.i <- sum(OR.n.mat[,i], na.rm=T)
    OR.pred.red.vec[i] <- OR.a.lp/(1+(OR.totN.i/OR.b.lp.vec[i])^OR.c.lp)
    
    diag(OR.popmat[2:(OR.age.max+1),]) <- (OR.s.stoch[-(OR.age.max+1)])*OR.pred.red.vec[i]
    OR.popmat[OR.age.max+1,OR.age.max+1] <- (OR.s.stoch[OR.age.max+1])*OR.pred.red.vec[i]
    OR.popmat[1,] <- ifelse(OR.fert.stch < 0, 0, OR.fert.stch)
    OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]
    
  } # end i loop
  
  OR.n.sums.mat[e,] <- (as.vector(colSums(OR.n.mat)))
  OR.pred.red.mn[e] <- mean(OR.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
OR.n.sums.matb <- OR.n.sums.mat[, -(1:round(OR.gen.l*2, 0))]

# total N
OR.n.md <- apply(OR.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
OR.n.up <- apply(OR.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
OR.n.lo <- apply(OR.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(OR.n.sums.matb)[2])
plot(yrs, OR.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(OR.n.lo), max(OR.n.up)))
lines(yrs,OR.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,OR.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(OR.n.md[2:length(yrs)] / OR.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("ORKstoch.RData")
rm(list = ls())



## NOTAMACROPUS (NR)
load("NR.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*NR.gen.l, 0) - 1)
NR.b.lp.vec2 <- rnorm(n=t, mean=NR.b.lp, sd=0.05*NR.b.lp)
NR.b.lp.vec <- ifelse(NR.b.lp.vec2 < 10, 10, NR.b.lp.vec2)
plot(NR.b.lp.vec, type="l", lty=2)
mean(log(NR.b.lp.vec[2:length(NR.b.lp.vec)] / NR.b.lp.vec[1:(length(NR.b.lp.vec)-1)]))

## set storage matrices & vectors
NR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
NR.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  NR.popmat <- NR.popmat.orig
  
  NR.n.mat <- matrix(0, nrow=NR.age.max+1,ncol=(t+1))
  NR.n.mat[,1] <- NR.init.vec
  NR.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    NR.s.alpha <- estBetaParams(NR.Sx, NR.s.sd.vec^2)$alpha
    NR.s.beta <- estBetaParams(NR.Sx, NR.s.sd.vec^2)$beta
    NR.s.stoch <- rbeta(length(NR.s.alpha), NR.s.alpha, NR.s.beta)
    
    if (rbinom(1, 1, 0.14/NR.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      NR.s.stoch <- NR.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    NR.fert.stch <- rnorm(length(NR.popmat[,1]), NR.pred.p.mm, NR.m.sd.vec)
    
    NR.totN.i <- sum(NR.n.mat[,i], na.rm=T)
    NR.pred.red.vec[i] <- NR.a.lp/(1+(NR.totN.i/NR.b.lp.vec[i])^NR.c.lp)
    
    diag(NR.popmat[2:(NR.age.max+1),]) <- (NR.s.stoch[-(NR.age.max+1)])*NR.pred.red.vec[i]
    NR.popmat[NR.age.max+1,NR.age.max+1] <- (NR.s.stoch[NR.age.max+1])*NR.pred.red.vec[i]
    NR.popmat[1,] <- ifelse(NR.fert.stch < 0, 0, NR.fert.stch)
    NR.n.mat[,i+1] <- NR.popmat %*% NR.n.mat[,i]
    
  } # end i loop
  
  NR.n.sums.mat[e,] <- (as.vector(colSums(NR.n.mat)))
  NR.pred.red.mn[e] <- mean(NR.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
NR.n.sums.matb <- NR.n.sums.mat[, -(1:round(NR.gen.l*2, 0))]

# total N
NR.n.md <- apply(NR.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
NR.n.up <- apply(NR.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
NR.n.lo <- apply(NR.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(NR.n.sums.matb)[2])
plot(yrs, NR.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(NR.n.lo), max(NR.n.up)))
lines(yrs,NR.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,NR.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(NR.n.md[2:length(yrs)] / NR.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("NRKstoch.RData")
rm(list = ls())



## GENYORNIS (GN)
load("GN.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*GN.gen.l, 0) - 1)
GN.b.lp.vec2 <- rnorm(n=t, mean=GN.b.lp, sd=0.05*GN.b.lp)
GN.b.lp.vec <- ifelse(GN.b.lp.vec2 < 10, 10, GN.b.lp.vec2)
plot(GN.b.lp.vec, type="l", lty=2)
mean(log(GN.b.lp.vec[2:length(GN.b.lp.vec)] / GN.b.lp.vec[1:(length(GN.b.lp.vec)-1)]))

## set storage matrices & vectors
GN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
GN.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  GN.popmat <- GN.popmat.orig
  
  GN.n.mat <- matrix(0, nrow=GN.age.max+1,ncol=(t+1))
  GN.n.mat[,1] <- GN.init.vec
  GN.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    GN.s.alpha <- estBetaParams(GN.Sx, GN.s.sd.vec^2)$alpha
    GN.s.beta <- estBetaParams(GN.Sx, GN.s.sd.vec^2)$beta
    GN.s.stoch <- rbeta(length(GN.s.alpha), GN.s.alpha, GN.s.beta)
    
    if (rbinom(1, 1, 0.14/GN.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      GN.s.stoch <- GN.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    GN.fert.stch <- rnorm(length(GN.popmat[,1]), GN.pred.p.mm, GN.m.sd.vec)
    
    GN.totN.i <- sum(GN.n.mat[,i], na.rm=T)
    GN.pred.red.vec[i] <- GN.a.lp/(1+(GN.totN.i/GN.b.lp.vec[i])^GN.c.lp)
    
    diag(GN.popmat[2:(GN.age.max+1),]) <- (GN.s.stoch[-(GN.age.max+1)])*GN.pred.red.vec[i]
    GN.popmat[GN.age.max+1,GN.age.max+1] <- (GN.s.stoch[GN.age.max+1])*GN.pred.red.vec[i]
    GN.popmat[1,] <- ifelse(GN.fert.stch < 0, 0, GN.fert.stch)
    GN.n.mat[,i+1] <- GN.popmat %*% GN.n.mat[,i]
    
  } # end i loop
  
  GN.n.sums.mat[e,] <- (as.vector(colSums(GN.n.mat)))
  GN.pred.red.mn[e] <- mean(GN.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
GN.n.sums.matb <- GN.n.sums.mat[, -(1:round(GN.gen.l*2, 0))]

# total N
GN.n.md <- apply(GN.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
GN.n.up <- apply(GN.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
GN.n.lo <- apply(GN.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(GN.n.sums.matb)[2])
plot(yrs, GN.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(GN.n.lo), max(GN.n.up)))
lines(yrs,GN.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,GN.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(GN.n.md[2:length(yrs)] / GN.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("GNKstoch.RData")
rm(list = ls())



## DROMAIUS (DN)
load("DN.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*DN.gen.l, 0) - 1)
DN.b.lp.vec2 <- rnorm(n=t, mean=DN.b.lp, sd=0.05*DN.b.lp)
DN.b.lp.vec <- ifelse(DN.b.lp.vec2 < 10, 10, DN.b.lp.vec2)
plot(DN.b.lp.vec, type="l", lty=2)
mean(log(DN.b.lp.vec[2:length(DN.b.lp.vec)] / DN.b.lp.vec[1:(length(DN.b.lp.vec)-1)]))

## set storage matrices & vectors
DN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DN.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  DN.popmat <- DN.popmat.orig
  
  DN.n.mat <- matrix(0, nrow=DN.age.max+1,ncol=(t+1))
  DN.n.mat[,1] <- DN.init.vec
  DN.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    DN.s.alpha <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$alpha
    DN.s.beta <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$beta
    DN.s.stoch <- rbeta(length(DN.s.alpha), DN.s.alpha, DN.s.beta)
    
    if (rbinom(1, 1, 0.14/DN.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      DN.s.stoch <- DN.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    DN.fert.stch <- rnorm(length(DN.popmat[,1]), DN.pred.p.mm, DN.m.sd.vec)
    
    DN.totN.i <- sum(DN.n.mat[,i], na.rm=T)
    DN.pred.red.vec[i] <- DN.a.lp/(1+(DN.totN.i/DN.b.lp.vec[i])^DN.c.lp)
    
    diag(DN.popmat[2:(DN.age.max+1),]) <- (DN.s.stoch[-(DN.age.max+1)])*DN.pred.red.vec[i]
    DN.popmat[DN.age.max+1,DN.age.max+1] <- (DN.s.stoch[DN.age.max+1])*DN.pred.red.vec[i]
    DN.popmat[1,] <- ifelse(DN.fert.stch < 0, 0, DN.fert.stch)
    DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]
    
  } # end i loop
  
  DN.n.sums.mat[e,] <- (as.vector(colSums(DN.n.mat)))
  DN.pred.red.mn[e] <- mean(DN.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
DN.n.sums.matb <- DN.n.sums.mat[, -(1:round(DN.gen.l*2, 0))]

# total N
DN.n.md <- apply(DN.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
DN.n.up <- apply(DN.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DN.n.lo <- apply(DN.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(DN.n.sums.matb)[2])
plot(yrs, DN.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(DN.n.lo), max(DN.n.up)))
lines(yrs,DN.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,DN.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(DN.n.md[2:length(yrs)] / DN.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("DNKstoch.RData")
rm(list = ls())



## ALECTURA (AL)
load("AL.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*AL.gen.l, 0) - 1)
AL.b.lp.vec2 <- rnorm(n=t, mean=AL.b.lp, sd=0.05*AL.b.lp)
AL.b.lp.vec <- ifelse(AL.b.lp.vec2 < 10, 10, AL.b.lp.vec2)
plot(AL.b.lp.vec, type="l", lty=2)
mean(log(AL.b.lp.vec[2:length(AL.b.lp.vec)] / AL.b.lp.vec[1:(length(AL.b.lp.vec)-1)]))

## set storage matrices & vectors
AL.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
AL.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  AL.popmat <- AL.popmat.orig
  
  AL.n.mat <- matrix(0, nrow=AL.age.max+1,ncol=(t+1))
  AL.n.mat[,1] <- AL.init.vec
  AL.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    AL.s.alpha <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$alpha
    AL.s.beta <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$beta
    AL.s.stoch <- rbeta(length(AL.s.alpha), AL.s.alpha, AL.s.beta)
    
    if (rbinom(1, 1, 0.14/AL.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      AL.s.stoch <- AL.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    AL.fert.stch <- rnorm(length(AL.popmat[,1]), AL.pred.p.mm, AL.m.sd.vec)
    
    AL.totN.i <- sum(AL.n.mat[,i], na.rm=T)
    AL.pred.red.vec[i] <- AL.a.lp/(1+(AL.totN.i/AL.b.lp.vec[i])^AL.c.lp)
    
    diag(AL.popmat[2:(AL.age.max+1),]) <- (AL.s.stoch[-(AL.age.max+1)])*AL.pred.red.vec[i]
    AL.popmat[AL.age.max+1,AL.age.max+1] <- (AL.s.stoch[AL.age.max+1])*AL.pred.red.vec[i]
    AL.popmat[1,] <- ifelse(AL.fert.stch < 0, 0, AL.fert.stch)
    AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
    
  } # end i loop
  
  AL.n.sums.mat[e,] <- (as.vector(colSums(AL.n.mat)))
  AL.pred.red.mn[e] <- mean(AL.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
AL.n.sums.matb <- AL.n.sums.mat[, -(1:round(AL.gen.l*2, 0))]

# total N
AL.n.md <- apply(AL.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
AL.n.up <- apply(AL.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
AL.n.lo <- apply(AL.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(AL.n.sums.matb)[2])
plot(yrs, AL.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(AL.n.lo), max(AL.n.up)))
lines(yrs,AL.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,AL.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(AL.n.md[2:length(yrs)] / AL.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("ALKstoch.RData")
rm(list = ls())



## THYLACOLEO (TC)
load("TC.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*TC.gen.l, 0) - 1)
TC.b.lp.vec2 <- rnorm(n=t, mean=TC.b.lp, sd=0.05*TC.b.lp)
TC.b.lp.vec <- ifelse(TC.b.lp.vec2 < 10, 10, TC.b.lp.vec2)
plot(TC.b.lp.vec, type="l", lty=2)
mean(log(TC.b.lp.vec[2:length(TC.b.lp.vec)] / TC.b.lp.vec[1:(length(TC.b.lp.vec)-1)]))

## set storage matrices & vectors
TC.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TC.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  TC.popmat <- TC.popmat.orig
  
  TC.n.mat <- matrix(0, nrow=TC.age.max+1,ncol=(t+1))
  TC.n.mat[,1] <- TC.init.vec
  TC.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    TC.s.alpha <- estBetaParams(TC.Sx, TC.s.sd.vec^2)$alpha
    TC.s.beta <- estBetaParams(TC.Sx, TC.s.sd.vec^2)$beta
    TC.s.stoch <- rbeta(length(TC.s.alpha), TC.s.alpha, TC.s.beta)
    
    if (rbinom(1, 1, 0.14/TC.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      TC.s.stoch <- TC.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    TC.fert.stch <- rnorm(length(TC.popmat[,1]), TC.pred.p.mm, TC.m.sd.vec)
    
    TC.totN.i <- sum(TC.n.mat[,i], na.rm=T)
    TC.pred.red.vec[i] <- TC.a.lp/(1+(TC.totN.i/TC.b.lp.vec[i])^TC.c.lp)
    
    diag(TC.popmat[2:(TC.age.max+1),]) <- (TC.s.stoch[-(TC.age.max+1)])*TC.pred.red.vec[i]
    TC.popmat[TC.age.max+1,TC.age.max+1] <- (TC.s.stoch[TC.age.max+1])*TC.pred.red.vec[i]
    TC.popmat[1,] <- ifelse(TC.fert.stch < 0, 0, TC.fert.stch)
    TC.n.mat[,i+1] <- TC.popmat %*% TC.n.mat[,i]
    
  } # end i loop
  
  TC.n.sums.mat[e,] <- (as.vector(colSums(TC.n.mat)))
  TC.pred.red.mn[e] <- mean(TC.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
TC.n.sums.matb <- TC.n.sums.mat[, -(1:round(TC.gen.l*2, 0))]

# total N
TC.n.md <- apply(TC.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
TC.n.up <- apply(TC.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TC.n.lo <- apply(TC.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(TC.n.sums.matb)[2])
plot(yrs, TC.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(TC.n.lo), max(TC.n.up)))
lines(yrs,TC.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,TC.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(TC.n.md[2:length(yrs)] / TC.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("TCKstoch.RData")
rm(list = ls())



## THYLACINUS (TH)
load("TH.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*TH.gen.l, 0) - 1)
TH.b.lp.vec2 <- rnorm(n=t, mean=TH.b.lp, sd=0.05*TH.b.lp)
TH.b.lp.vec <- ifelse(TH.b.lp.vec2 < 10, 10, TH.b.lp.vec2)
plot(TH.b.lp.vec, type="l", lty=2)
mean(log(TH.b.lp.vec[2:length(TH.b.lp.vec)] / TH.b.lp.vec[1:(length(TH.b.lp.vec)-1)]))

## set storage matrices & vectors
TH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TH.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  TH.popmat <- TH.popmat.orig
  
  TH.n.mat <- matrix(0, nrow=TH.age.max+1,ncol=(t+1))
  TH.n.mat[,1] <- TH.init.vec
  TH.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    TH.s.alpha <- estBetaParams(TH.Sx, TH.s.sd.vec^2)$alpha
    TH.s.beta <- estBetaParams(TH.Sx, TH.s.sd.vec^2)$beta
    TH.s.stoch <- rbeta(length(TH.s.alpha), TH.s.alpha, TH.s.beta)
    
    if (rbinom(1, 1, 0.14/TH.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      TH.s.stoch <- TH.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    TH.fert.stch <- rnorm(length(TH.popmat[,1]), TH.pred.p.mm, TH.m.sd.vec)
    
    TH.totN.i <- sum(TH.n.mat[,i], na.rm=T)
    TH.pred.red.vec[i] <- TH.a.lp/(1+(TH.totN.i/TH.b.lp.vec[i])^TH.c.lp)
    
    diag(TH.popmat[2:(TH.age.max+1),]) <- (TH.s.stoch[-(TH.age.max+1)])*TH.pred.red.vec[i]
    TH.popmat[TH.age.max+1,TH.age.max+1] <- 0
    TH.popmat[1,] <- ifelse(TH.fert.stch < 0, 0, TH.fert.stch)
    TH.n.mat[,i+1] <- TH.popmat %*% TH.n.mat[,i]
    
  } # end i loop
  
  TH.n.sums.mat[e,] <- (as.vector(colSums(TH.n.mat)))
  TH.pred.red.mn[e] <- mean(TH.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
TH.n.sums.matb <- TH.n.sums.mat[, -(1:round(TH.gen.l*2, 0))]

# total N
TH.n.md <- apply(TH.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
TH.n.up <- apply(TH.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TH.n.lo <- apply(TH.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(TH.n.sums.matb)[2])
plot(yrs, TH.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(TH.n.lo), max(TH.n.up)))
lines(yrs,TH.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,TH.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(TH.n.md[2:length(yrs)] / TH.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("THKstoch.RData")
rm(list = ls())



## SARCOPHILUS (SH)
load("SH.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*SH.gen.l, 0) - 1)
SH.b.lp.vec2 <- rnorm(n=t, mean=SH.b.lp, sd=0.05*SH.b.lp)
SH.b.lp.vec <- ifelse(SH.b.lp.vec2 < 10, 10, SH.b.lp.vec2)
plot(SH.b.lp.vec, type="l", lty=2)
mean(log(SH.b.lp.vec[2:length(SH.b.lp.vec)] / SH.b.lp.vec[1:(length(SH.b.lp.vec)-1)]))

## set storage matrices & vectors
SH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
SH.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  SH.popmat <- SH.popmat.orig
  
  SH.n.mat <- matrix(0, nrow=SH.age.max+1,ncol=(t+1))
  SH.n.mat[,1] <- SH.init.vec
  SH.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    SH.s.alpha <- estBetaParams(SH.Sx, SH.s.sd.vec^2)$alpha
    SH.s.beta <- estBetaParams(SH.Sx, SH.s.sd.vec^2)$beta
    SH.s.stoch <- rbeta(length(SH.s.alpha), SH.s.alpha, SH.s.beta)
    
    if (rbinom(1, 1, 0.14/SH.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      SH.s.stoch <- SH.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    SH.fert.stch <- rnorm(length(SH.popmat[,1]), SH.m.vec, SH.m.sd.vec)
    
    SH.totN.i <- sum(SH.n.mat[,i], na.rm=T)
    SH.pred.red.vec[i] <- SH.a.lp/(1+(SH.totN.i/SH.b.lp.vec[i])^SH.c.lp)
    
    diag(SH.popmat[2:(SH.age.max+1),]) <- (SH.s.stoch[-1])*SH.pred.red.vec[i]
    SH.popmat[SH.age.max+1,SH.age.max+1] <- 0
    SH.popmat[1,] <- ifelse(SH.fert.stch < 0, 0, SH.fert.stch)
    SH.n.mat[,i+1] <- SH.popmat %*% SH.n.mat[,i]
    
  } # end i loop
  
  SH.n.sums.mat[e,] <- (as.vector(colSums(SH.n.mat)))
  SH.pred.red.mn[e] <- mean(SH.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
SH.n.sums.matb <- SH.n.sums.mat[, -(1:round(SH.gen.l*2, 0))]

# total N
SH.n.md <- apply(SH.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
SH.n.up <- apply(SH.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
SH.n.lo <- apply(SH.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(SH.n.sums.matb)[2])
plot(yrs, SH.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(SH.n.lo), max(SH.n.up)))
lines(yrs,SH.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,SH.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(SH.n.md[2:length(yrs)] / SH.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("SHKstoch.RData")
rm(list = ls())



## DASYURUS (DM)
load("DM.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*DM.gen.l, 0) - 1)
DM.b.lp.vec2 <- rnorm(n=t, mean=DM.b.lp, sd=0.05*DM.b.lp)
DM.b.lp.vec <- ifelse(DM.b.lp.vec2 < 10, 10, DM.b.lp.vec2)
plot(DM.b.lp.vec, type="l", lty=2)
mean(log(DM.b.lp.vec[2:length(DM.b.lp.vec)] / DM.b.lp.vec[1:(length(DM.b.lp.vec)-1)]))

## set storage matrices & vectors
DM.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DM.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  DM.popmat <- DM.popmat.orig
  
  DM.n.mat <- matrix(0, nrow=DM.age.max+1,ncol=(t+1))
  DM.n.mat[,1] <- DM.init.vec
  DM.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    DM.s.alpha <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$alpha
    DM.s.beta <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$beta
    DM.s.stoch <- rbeta(length(DM.s.alpha), DM.s.alpha, DM.s.beta)
    
    if (rbinom(1, 1, 0.14/DM.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      DM.s.stoch <- DM.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    DM.fert.stch <- rnorm(length(DM.popmat[,1]), DM.m.vec, DM.m.sd.vec)
    
    DM.totN.i <- sum(DM.n.mat[,i], na.rm=T)
    DM.pred.red.vec[i] <- DM.a.lp/(1+(DM.totN.i/DM.b.lp.vec[i])^DM.c.lp)
    
    diag(DM.popmat[2:(DM.age.max+1),]) <- (DM.s.stoch[-1])*DM.pred.red.vec[i]
    DM.popmat[DM.age.max+1,DM.age.max+1] <- 0
    DM.popmat[1,] <- ifelse(DM.fert.stch < 0, 0, DM.fert.stch)
    DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
    
  } # end i loop
  
  DM.n.sums.mat[e,] <- (as.vector(colSums(DM.n.mat)))
  DM.pred.red.mn[e] <- mean(DM.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
DM.n.sums.matb <- DM.n.sums.mat[, -(1:round(DM.gen.l*2, 0))]

# total N
DM.n.md <- apply(DM.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
DM.n.up <- apply(DM.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DM.n.lo <- apply(DM.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(DM.n.sums.matb)[2])
plot(yrs, DM.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(DM.n.lo), max(DM.n.up)))
lines(yrs,DM.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,DM.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(DM.n.md[2:length(yrs)] / DM.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("DMKstoch.RData")
rm(list = ls())



## TACHYGLOSSUS (TA)
load("TA.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*TA.gen.l, 0) - 1)
TA.b.lp.vec2 <- rnorm(n=t, mean=TA.b.lp, sd=0.05*TA.b.lp)
TA.b.lp.vec <- ifelse(TA.b.lp.vec2 < 10, 10, TA.b.lp.vec2)
plot(TA.b.lp.vec, type="l", lty=2)
mean(log(TA.b.lp.vec[2:length(TA.b.lp.vec)] / TA.b.lp.vec[1:(length(TA.b.lp.vec)-1)]))

## set storage matrices & vectors
TA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TA.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  TA.popmat <- TA.popmat.orig
  
  TA.n.mat <- matrix(0, nrow=TA.age.max+1,ncol=(t+1))
  TA.n.mat[,1] <- TA.init.vec
  TA.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    TA.s.alpha <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$alpha
    TA.s.beta <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$beta
    TA.s.stoch <- rbeta(length(TA.s.alpha), TA.s.alpha, TA.s.beta)
    
    if (rbinom(1, 1, 0.14/TA.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      TA.s.stoch <- TA.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    TA.fert.stch <- rnorm(length(TA.popmat[,1]), TA.pred.p.mm, TA.m.sd.vec)
    
    TA.totN.i <- sum(TA.n.mat[,i], na.rm=T)
    TA.pred.red.vec[i] <- TA.a.lp/(1+(TA.totN.i/TA.b.lp.vec[i])^TA.c.lp)
    
    diag(TA.popmat[2:(TA.age.max+1),]) <- (TA.s.stoch[-(TA.age.max+1)])*TA.pred.red.vec[i]
    TA.popmat[TA.age.max+1,TA.age.max+1] <- 0
    TA.popmat[1,] <- ifelse(TA.fert.stch < 0, 0, TA.fert.stch)
    TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
    
  } # end i loop
  
  TA.n.sums.mat[e,] <- (as.vector(colSums(TA.n.mat)))
  TA.pred.red.mn[e] <- mean(TA.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
TA.n.sums.matb <- TA.n.sums.mat[, -(1:round(TA.gen.l*2, 0))]

# total N
TA.n.md <- apply(TA.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
TA.n.up <- apply(TA.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TA.n.lo <- apply(TA.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(TA.n.sums.matb)[2])
plot(yrs, TA.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(TA.n.lo), max(TA.n.up)))
lines(yrs,TA.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,TA.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(TA.n.md[2:length(yrs)] / TA.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("TAKstoch.RData")
rm(list = ls())



## MEGALIBGWILIA (MR)
load("MR.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*MR.gen.l, 0) - 1)
MR.b.lp.vec2 <- rnorm(n=t, mean=MR.b.lp, sd=0.05*MR.b.lp)
MR.b.lp.vec <- ifelse(MR.b.lp.vec2 < 10, 10, MR.b.lp.vec2)
plot(MR.b.lp.vec, type="l", lty=2)
mean(log(MR.b.lp.vec[2:length(MR.b.lp.vec)] / MR.b.lp.vec[1:(length(MR.b.lp.vec)-1)]))

## set storage matrices & vectors
MR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
MR.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  MR.popmat <- MR.popmat.orig
  
  MR.n.mat <- matrix(0, nrow=MR.age.max+1,ncol=(t+1))
  MR.n.mat[,1] <- MR.init.vec
  MR.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    MR.s.alpha <- estBetaParams(MR.Sx, MR.s.sd.vec^2)$alpha
    MR.s.beta <- estBetaParams(MR.Sx, MR.s.sd.vec^2)$beta
    MR.s.stoch <- rbeta(length(MR.s.alpha), MR.s.alpha, MR.s.beta)
    
    if (rbinom(1, 1, 0.14/MR.gen.l) == 1) { # catastrophe
      cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
      cat.beta <- estBetaParams(0.5, 0.05^2)$beta
      MR.s.stoch <- MR.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
    
    # stochastic fertility sampler (Gaussian)
    MR.fert.stch <- rnorm(length(MR.popmat[,1]), MR.pred.p.mm, MR.m.sd.vec)
    
    MR.totN.i <- sum(MR.n.mat[,i], na.rm=T)
    MR.pred.red.vec[i] <- MR.a.lp/(1+(MR.totN.i/MR.b.lp.vec[i])^MR.c.lp)
    
    diag(MR.popmat[2:(MR.age.max+1),]) <- (MR.s.stoch[-(MR.age.max+1)])*MR.pred.red.vec[i]
    MR.popmat[MR.age.max+1,MR.age.max+1] <- 0
    MR.popmat[1,] <- ifelse(MR.fert.stch < 0, 0, MR.fert.stch)
    MR.n.mat[,i+1] <- MR.popmat %*% MR.n.mat[,i]
    
  } # end i loop
  
  MR.n.sums.mat[e,] <- (as.vector(colSums(MR.n.mat)))
  MR.pred.red.mn[e] <- mean(MR.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

# remove first 2 generations as burn-in
MR.n.sums.matb <- MR.n.sums.mat[, -(1:round(MR.gen.l*2, 0))]

# total N
MR.n.md <- apply(MR.n.sums.matb, MARGIN=2, median, na.rm=T) # mean over all iterations
MR.n.up <- apply(MR.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
MR.n.lo <- apply(MR.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

# plot
yrs <- 1:(dim(MR.n.sums.matb)[2])
plot(yrs, MR.n.md, type="l", lty=1, lwd=0.8, ylim=c(min(MR.n.lo), max(MR.n.up)))
lines(yrs,MR.n.up, lty=2, col="red", lwd=0.6)
lines(yrs,MR.n.lo, lty=2, col="red", lwd=0.6)

# mean rate of decline
rmed1 <- na.omit(log(MR.n.md[2:length(yrs)] / MR.n.md[1:(length(yrs)-1)]))
which.inf <- which(is.infinite(rmed1) == T)
if (length(which.inf) == 0) {
  mean(rmed1)
} else {
  mean(rmed1[-which.inf])
}
print("target = -0.01")
save.image("MRKstoch.RData")
rm(list = ls())
