#############################################################################################################################
## density-feedback simulations
## base Leslie matrix models
## assumptions check - effect of increasing variation in juvenile survival relative to adults
## Diprotodon only
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

## body mass estimates
DP.mass <- 2786 # Diprotodon optatum (Wroe et al. 2004 PRSB 271:S34-S36)
PA.mass <- 1000 # Palorchestes azael (Richards et al. 2019 (PLoS One https://doi.org/10.1371/journal.pone.0221824)
ZT.mass <- 500 # Zygomaturus trilobus (Johnson et al. 2006 Australia’s Mammal Extinctions. 278 pp. Cambridge  University Press, Melbourne)
PH.mass <- 200 # Phascolonus gigas (Johnson et al. 2006)
VU.mass <- 25 # Vombatus ursinus (Saran et al. 2011 Pacific Conservation Biology 17:310-319)

PG.mass <- 250 # Procoptodon goliah (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1442-9993.2004.01389.x)
SS.mass <- 150 # Sthenurus stirlingi (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1442-9993.2004.01389.x)
PT.mass <- 130 # Protemnodon anak (Johnson et al. 2006)
SO.mass <- 120 # Simosthenurus occidentalis (Johnson 2006)
MN.mass <- 55 # Metasthenurus newtonae (Johnson et al. 2006)
OR.mass <- 25 # female red kangaroo (Osphranter rufus) Croft, D.B. & Clacy, T.F. (2008) Macropus rufus. The Mammals of Australia (eds S. van Dyck & R. Strahan), pp. 352-354. Reed New Holland, Sydney
NR.mass <- 14 # Notamacropus rufogriseus females; Strahan R. (1983) The Australian Museum Complete Book of Australian Mammals, 1st edn. Angus & Robertson Publishers, Sydney

GN.mass <- 200 # Genyornis newtoni (max) (Johnson et al. 2006)
DN.mass <- 55 # Dromaius novaehollandiae (Sales et al. 2007 Avian and Poultry Biology Reviews 18:1–20)
AL.mass <- 2.2 # Alectura lathami (brush turkey) adult female (Jones, R.W.R.J. Dekker, C.S. Roselaar. The Megapodes: Megapodiidae, Oxford University Press, Oxford, New York, Tokyo (1995))

TC.mass <- 110 # Thylacoleo carnifex (Johnson et al. 2006)
TH.mass <- 20 # Thylacinus (Jones and Stoddart 1998 Journal of Zoology 246:239-246; Lachish et al. 2009)
SH.mass <- 6.1 # Sarcophilus (female; Guiler 1978 in Bradshaw & Brook 2005)
DM.mass <- 2.0 # 2.0 from Belcher et al. AM Mammals (3rd edition; 2008); 1.68/1.7 Dasyurus maculatus spotted-tailed quolls (Körtner et al. 2004-Wildl Res; Green 2008)

TA.mass <- 4.0 # mean(c(3.8,3.4,3.7,(mean(c(3.9,7))))) # short-beaked echidna Tachyglossus aculeatus (Nicol & Andersen 2007)
MR.mass <- 11 # 10-12 kg Megalibgwilia ramsayi (Johnson 2006)

## remove everything
#rm(list = ls())

# libraries
library(plotly)
library(ggpubr)

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# damping ratio (ρ): dimensionless measure of convergence to stable growth; smaller numbers represent slower convergence (Capdevila et al. 2020-TREE)
# "... populations in which reproduction is spread out symmetrically over a wide range of ages should converge to the stable age distribution
# more rapidaly than those with tightly restricted ages at reproduction, or skewed distributionsof age at reproduction" (Caswell 2001, p. 98)
damping.ratio <- function(x) {
  eigvals <- eigen(x)$values
  lmax <- which.max(Re(eigvals))
  lambda1 <- Re(eigvals[lmax])
  lambda2 <- sort(Mod(eigvals), decreasing=T)[2]
  dr <- lambda1/lambda2
  return(dr)}

AICc <- function(...) {
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
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


######################################################
## BASE MODELS
######################################################

# MACROPOD CORRECTIONS
# mass
NR.mass <- 5.1 # Notamacropus rufogriseus
#DL.mass <- 6.48 # Dendrolagus lumholtzi
#DMa.mass <- 9.2 # Dendrolagus matscheiei
#DB.mass <- 9.3 # Dendrolagus bennettianus
OR.mass <- 25 # Osphranter rufus
DOL.mass <- 3.57 # Dorcopsis luctuosa
LC.mass <- 3 # Lagorchestes conspicillatus
LH.mass <- 1.3 # L. hirsutus
LF.mass <- 1.8 # L. fasciatus
MAg.mass <- 11 # Macropus agilis
MAn.mass <- 17.5 # M. antilopus
MDo.mass <- 6.5 # M. dorsalis
#MEu.mass <- 5.5 # M. eugenii
MFu.mass <- 16 # M. fuliginosus
MGi.mass <- 17.8 # M. giganteus
MPa.mass <- 3.55 # M. parma
MParr.mass <- 11 # M. parryi
MRo.mass <- 15.6 # M. robustus
PAs.mass <- 4.3 # Petrogale assimilis
PI.mass <- 4.2 # P. inomata
PP.mass <- 6.3 # P. penicillata
PPe.mass <- 5.2 # P. persephone
PX.mass <- 7 # P. xanthopus
TB.mass <- 3.9 # Thylogale billardienii
TS.mass <- 4.1 # T. stigmatica
TT.mass <- 3.8 # T. thetis
WB.mass <- 13.0 # Wallabia bicolor

# inter-birth interval
NR.IBI <- 286
#DL.IBI <- 435
#DMa.IBI <- 410
#DB.IBI <- 365
OR.IBI <- 241
DOL.IBI <- 191
LC.IBI <- 153
LH.IBI <- 125
LF.IBI <- 365
MAg.IBI <- 220
MAn.IBI <- 270
MDo.IBI <- 211
#MEu.IBI <- 365
MFu.IBI <- 372
MGi.IBI <- 363
MPa.IBI <- 213
MParr.IBI <- 266
MRo.IBI <- mean(c(256,264))
PAs.IBI <- 200
PI.IBI <- 210
PP.IBI <- 205
PPe.IBI <- 209
PX.IBI <- 196
TB.IBI <- 204
TS.IBI <- 185
TT.IBI <- 182
WB.IBI <- 256

MACROPOD.IBI <- c(NR.IBI,OR.IBI,DOL.IBI,LC.IBI,LH.IBI,LF.IBI,MAg.IBI,MAn.IBI,MDo.IBI,MFu.IBI,MGi.IBI,MPa.IBI,MParr.IBI,MRo.IBI,PAs.IBI,PI.IBI,PP.IBI,PPe.IBI,PX.IBI,TB.IBI,TS.IBI,TT.IBI,WB.IBI)
MACROPOD.mass <- c(NR.mass,OR.mass,DOL.mass,LC.mass,LH.mass,LF.mass,MAg.mass,MAn.mass,MDo.mass,MFu.mass,MGi.mass,MPa.mass,MParr.mass,MRo.mass,PAs.mass,PI.mass,PP.mass,PPe.mass,PX.mass,TB.mass,TS.mass,TT.mass,WB.mass)

plot(log10(MACROPOD.mass), MACROPOD.IBI, pch=19)
linreg.ER(log10(MACROPOD.mass), MACROPOD.IBI)
IBImass.fit <- lm(MACROPOD.IBI ~ log10(MACROPOD.mass))
summary(IBImass.fit)
abline(h=mean(MACROPOD.IBI), lty=2)
abline(IBImass.fit, lty=2, col="red")

MACROPOD.F.corr1 <- 365/mean(MACROPOD.IBI)
MACROPOD.F.corr1

MACROPOD.F.corr.a <- coef(IBImass.fit)[1]
MACROPOD.F.corr.b <- coef(IBImass.fit)[2]



# age at first breeding (from Fisher et al. 2001)
NR.alpha <- mean(c(412,420,390))
#DL.alpha <- mean(c(720))
#DMa.alpha <- mean(c(960,730))
OR.alpha <- mean(c(913,613))
DOL.alpha <- 450
LC.alpha <- 363
LH.alpha <- mean(c(315,345))
LF.alpha <- 365
MAg.alpha <- mean(c(405,357))
MAn.alpha <- 802
MDo.alpha <- 420
#MEu.alpha <- 270
MFu.alpha <- mean(c(604,420))
MGi.alpha <- mean(c(660,600))
MPa.alpha <- 424
MParr.alpha <- mean(c(630,586))
MRo.alpha <- mean(c(535,613,909))
PAs.alpha <- 525
PI.alpha <- 540
PP.alpha <- 540
PPe.alpha <- 683
PX.alpha <- 541
TB.alpha <- 392
TS.alpha <- mean(c(341,336))
TT.alpha <- mean(c(510,600))
WB.alpha <- 450

# predicted age at first breeding
NR.alpha.pred <- (exp(-1.34 + (0.214*log(NR.mass*1000))))
#DL.alpha.pred <- (exp(-1.34 + (0.214*log(DL.mass*1000))))
#DMa.alpha.pred <- (exp(-1.34 + (0.214*log(DMa.mass*1000))))
OR.alpha.pred <- (exp(-1.34 + (0.214*log(OR.mass*1000))))
DOL.alpha.pred <- (exp(-1.34 + (0.214*log(DOL.mass*1000))))
LC.alpha.pred <- (exp(-1.34 + (0.214*log(LC.mass*1000))))
LH.alpha.pred <- (exp(-1.34 + (0.214*log(LH.mass*1000))))
LF.alpha.pred <- (exp(-1.34 + (0.214*log(LF.mass*1000))))
MAg.alpha.pred <- (exp(-1.34 + (0.214*log(MAg.mass*1000))))
MAn.alpha.pred <- (exp(-1.34 + (0.214*log(MAn.mass*1000))))
MDo.alpha.pred <- (exp(-1.34 + (0.214*log(MDo.mass*1000))))
#MEu.alpha.pred <- (exp(-1.34 + (0.214*log(MEu.mass*1000))))
MFu.alpha.pred <- (exp(-1.34 + (0.214*log(MFu.mass*1000))))
MGi.alpha.pred <- (exp(-1.34 + (0.214*log(MGi.mass*1000))))
MPa.alpha.pred <- (exp(-1.34 + (0.214*log(MPa.mass*1000))))
MParr.alpha.pred <- (exp(-1.34 + (0.214*log(MParr.mass*1000))))
MRo.alpha.pred <- (exp(-1.34 + (0.214*log(MRo.mass*1000))))
PAs.alpha.pred <- (exp(-1.34 + (0.214*log(PAs.mass*1000))))
PI.alpha.pred <- (exp(-1.34 + (0.214*log(PI.mass*1000))))
PP.alpha.pred <- (exp(-1.34 + (0.214*log(PP.mass*1000))))
PPe.alpha.pred <- (exp(-1.34 + (0.214*log(PPe.mass*1000))))
PX.alpha.pred <- (exp(-1.34 + (0.214*log(PX.mass*1000))))
TB.alpha.pred <- (exp(-1.34 + (0.214*log(TB.mass*1000))))
TS.alpha.pred <- (exp(-1.34 + (0.214*log(TS.mass*1000))))
TT.alpha.pred <- (exp(-1.34 + (0.214*log(TT.mass*1000))))
WB.alpha.pred <- (exp(-1.34 + (0.214*log(WB.mass*1000))))


MACROPOD.alpha.vec <- c(NR.alpha,OR.alpha,DOL.alpha,LC.alpha,LH.alpha,LF.alpha,MAg.alpha,MAn.alpha,MDo.alpha,MFu.alpha,MGi.alpha,MPa.alpha,MParr.alpha,MRo.alpha,PAs.alpha,PI.alpha,PP.alpha,PPe.alpha,PX.alpha,TB.alpha,TS.alpha,TT.alpha,WB.alpha)
MACROPOD.alpha.pred.vec <- c(NR.alpha.pred,OR.alpha.pred,DOL.alpha.pred,LC.alpha.pred,LH.alpha.pred,LF.alpha.pred,MAg.alpha.pred,MAn.alpha.pred,MDo.alpha.pred,MFu.alpha.pred,MGi.alpha.pred,MPa.alpha.pred,MParr.alpha.pred,MRo.alpha.pred,PAs.alpha.pred,PI.alpha.pred,PP.alpha.pred,PPe.alpha.pred,PX.alpha.pred,TB.alpha.pred,TS.alpha.pred,TT.alpha.pred,WB.alpha.pred)
MACROPOD.alpha.vec.yr <- MACROPOD.alpha.vec/365
mean((MACROPOD.alpha.pred.vec - MACROPOD.alpha.vec.yr) / MACROPOD.alpha.pred.vec)
MACROPOD.alpha.corr <- mean(MACROPOD.alpha.vec)/365 / mean(MACROPOD.alpha.pred.vec)
MACROPOD.alpha.corr

plot(MACROPOD.mass, MACROPOD.alpha.vec, pch=19)


# VOMBATIFORM CORRECTIONS
# mass
PC.mass <- 5.1
VU.mass <- 25 
LL.mass <- 26
LK.mass <- 31

# inter-birth interval
PC.IBI <- 383 # koala interbirth interval (Martin, R.W. & Handasyde, K.A. 1995  Koala  In The Australian Museum complete book of Australian Mammals (ed. R. Strahan) pp.196-198. Sydney: Reed Books)
VU.IBI <- 730 # Vombatus ursinus
LL.IBI <- 365 # southern hairy-nosed wombat
LK.IBI <- 548 # northern hairy-nosed wombat

# age at first breeding (from Fisher et al. 2001)
PC.alpha <- 2
VU.alpha <- 2
LL.alpha <- round(mean(c(1095,540))/365, 0)

# predicted age at first breeding
PC.alpha.pred <- ceiling(exp(-1.34 + (0.214*log(PC.mass*1000))))
VU.alpha.pred <- ceiling(exp(-1.34 + (0.214*log(VU.mass*1000))))
LL.alpha.pred <- ceiling(exp(-1.34 + (0.214*log(LL.mass*1000))))

VOMBAT.alpha.corr <- mean(c(PC.alpha,VU.alpha,LL.alpha)) / mean(c(PC.alpha.pred,VU.alpha.pred,LL.alpha.pred))
VOMBAT.alpha.corr

VOMBAT.IBI <- c(PC.IBI,VU.IBI,LL.IBI,LK.IBI)
VOMBAT.mass <- c(PC.mass,VU.mass,LL.mass,LK.mass)

plot((VOMBAT.mass), VOMBAT.IBI, pch=19)
abline(h=mean(VOMBAT.IBI), lty=2)

VOMBAT.F.corr <- 365/mean(VOMBAT.IBI)
VOMBAT.F.corr

# use macropod relationship to project vombatiforms
mass.range.vec <- seq(1,DP.mass,1)
VOMBAT.IBI.pred <- mean(VOMBAT.IBI) + (MACROPOD.F.corr.b * log10(mass.range.vec))
plot(log10(mass.range.vec),VOMBAT.IBI.pred, type="l", ylim=c(min(VOMBAT.IBI),max(VOMBAT.IBI.pred)))
abline(h=mean(VOMBAT.IBI), lty=2)
points(log10(VOMBAT.mass), VOMBAT.IBI, pch=19)
VOMBAT.EXT.mass <- c(DP.mass, PA.mass, ZT.mass, PH.mass)
VOMBAT.IBI.EXT.pred <- mean(VOMBAT.IBI) + (MACROPOD.F.corr.b * log10(VOMBAT.EXT.mass))
points(log10(VOMBAT.EXT.mass), VOMBAT.IBI.EXT.pred, pch=19, col="red")

# as above, but anchor relationship to VU
VU.IBI.pred <- as.numeric(mean(VOMBAT.IBI) + (MACROPOD.F.corr.b * log10(VU.mass)))
VOMBAT.IBI.pred2 <- (mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(mass.range.vec))
plot(log10(mass.range.vec),VOMBAT.IBI.pred2, type="l", ylim=c(min(VOMBAT.IBI),max(VOMBAT.IBI.pred2)))
points(log10(VOMBAT.mass), VOMBAT.IBI, pch=19)
VOMBAT.IBI.EXT.pred2 <- (mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(VOMBAT.EXT.mass))
points(log10(VOMBAT.EXT.mass), VOMBAT.IBI.EXT.pred2, pch=19, col="red")



############################
## DIPROTODON (optatum) (DP)
## sources: Brook & Johnson 2006 (Alcheringa 30:39-48, http://doi.org/10.1080/03115510609506854)

# mass
DP.mass <- 2786 # kg (Wroe et al. 2004 PRSB 271:S34-S36)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
DP.rm.pred <- 10^(0.6914 - (0.2622*log10(DP.mass*1000)))
DP.rm.pred
DP.lm.pred <- exp(DP.rm.pred)
DP.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
DP.D.pred <- (10^(4.196 - (0.74*log10(DP.mass*1000))))/2 # divided by 2 for females only
DP.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
DP.age.max1 <- round(10^(0.89 + (0.13*log10(DP.mass*1000))), 0)
DP.age.max <- round(DP.age.max1 * 26/29, 0) # corrected for over-estimate derived from Vombatus

## age vector
DP.age.vec <- 0:DP.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
DP.F.pred1 <- exp(2.719 - (0.211*log(DP.mass*1000)))/2 # divided by 2 for females
DP.F.pred1
#DP.F.pred <- DP.F.pred1 * VOMBAT.F.corr
#DP.F.pred <- DP.F.pred1 * (365/as.numeric(mean(VOMBAT.IBI) + (MACROPOD.F.corr.b * log10(DP.mass))))
DP.F.pred <- DP.F.pred1 * (365/as.numeric((mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(DP.mass))))
DP.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
DP.alpha1 <- ceiling(exp(-1.34 + (0.214*log(DP.mass*1000))))
DP.alpha1
DP.alpha <- DP.alpha1


## define m function with age
DP.m.vec <- c(rep(0, DP.alpha-1), rep(0.75*DP.F.pred, round(DP.alpha/2,0)), rep(DP.F.pred, (DP.age.max+1-((DP.alpha-1+round(DP.alpha/2,0))))))
DP.m.sd.vec <- 0.05*DP.m.vec
plot(DP.age.vec, DP.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
DP.m.dat <- data.frame(DP.age.vec, DP.m.vec)
param.init <- c(0.3, 6, -5)
DP.fit.logp <- nls(DP.m.vec ~ a / (1+(DP.age.vec/b)^c), 
                data = DP.m.dat,
                algorithm = "port",
                start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DP.fit.logp.summ <- summary(DP.fit.logp)
plot(DP.age.vec, DP.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
DP.age.vec.cont <- seq(0,max(DP.age.vec),1)
DP.pred.p.m <- coef(DP.fit.logp)[1] / (1+(DP.age.vec.cont/coef(DP.fit.logp)[2])^coef(DP.fit.logp)[3])
DP.pred.p.mm <- ifelse(DP.pred.p.m > 1, 1, DP.pred.p.m)
lines(DP.age.vec.cont, DP.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
DP.s.tran <- ln.a.s + b.s*log(DP.mass*1000) + log(1)
DP.s.ad.yr <- exp(-exp(DP.s.tran))
DP.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1*DP.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.0 # rate of mortality decline (also known as bt)
a2 <- 1 - DP.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.9e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.05 # rate of mortality increase
longev <- DP.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
#plot(x,l.x,type="l")
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
DP.Sx <- c(0.99*DP.s.ad.yr, 1 - qx)
plot(x, DP.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
DP.sd.set <- 0.05
DP.sd.juv <- seq(3*DP.sd.set, DP.sd.set, -((3*DP.sd.set - DP.sd.set) / (DP.alpha+1)))
DP.sd.age <- c(DP.sd.juv, rep(DP.sd.set, length(DP.Sx) - length(DP.sd.juv)))
DP.s.sd.vec <- DP.sd.age*DP.Sx # assume 3 times more variation in juveniles compared to adults

## create matrix
DP.popmat <- matrix(data = 0, nrow=DP.age.max+1, ncol=DP.age.max+1)
diag(DP.popmat[2:(DP.age.max+1),]) <- DP.Sx[-(DP.age.max+1)]
DP.popmat[DP.age.max+1,DP.age.max+1] <- DP.Sx[DP.age.max+1]
DP.popmat[1,] <- DP.pred.p.mm
colnames(DP.popmat) <- c(0:DP.age.max)
rownames(DP.popmat) <- c(0:DP.age.max)
DP.popmat.orig <- DP.popmat ## save original matrix

## matrix properties
max.lambda(DP.popmat.orig) ## 1-yr lambda
DP.lm.pred
max.r(DP.popmat.orig) # rate of population change, 1-yr
DP.ssd <- stable.stage.dist(DP.popmat.orig) ## stable stage distribution
plot(DP.age.vec, DP.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(DP.popmat.orig, DP.age.max) # reproductive value
DP.gen.l <- G.val(DP.popmat.orig, DP.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km = 250,000 km^2; equates to approximately 10% larger than State of Victoria (227,444 km^2)
DP.pop.found <- round(area*DP.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
DP.init.vec <- DP.ssd * DP.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*DP.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

DP.tot.F <- sum(DP.popmat.orig[1,])
DP.popmat <- DP.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
DP.n.mat <- matrix(0, nrow=DP.age.max+1,ncol=(t+1))
DP.n.mat[,1] <- DP.init.vec

## set up projection loop
for (i in 1:t) {
  DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]
}

DP.n.pred <- colSums(DP.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(DP.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
DP.K.max <- 1*DP.pop.found
DP.K.vec <- c(1, DP.K.max/2, 0.75*DP.K.max, DP.K.max) 
DP.red.vec <- c(1,0.98,0.96,0.9383)
plot(DP.K.vec, DP.red.vec,pch=19,type="b")
DP.Kred.dat <- data.frame(DP.K.vec, DP.red.vec)

# logistic power function a/(1+(x/b)^c)
DP.param.init <- c(1, 2*DP.K.max, 2)
DP.fit.lp <- nls(DP.red.vec ~ a/(1+(DP.K.vec/b)^c), 
              data = DP.Kred.dat,
              algorithm = "port",
              start = c(a = DP.param.init[1], b = DP.param.init[2], c = DP.param.init[3]),
              trace = TRUE,      
              nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DP.fit.lp.summ <- summary(DP.fit.lp)
plot(DP.K.vec, DP.red.vec, pch=19,xlab="N",ylab="reduction factor")
DP.K.vec.cont <- seq(1,2*DP.pop.found,1)
DP.pred.lp.fx <- coef(DP.fit.lp)[1]/(1+(DP.K.vec.cont/coef(DP.fit.lp)[2])^coef(DP.fit.lp)[3])
lines(DP.K.vec.cont, DP.pred.lp.fx, lty=3,lwd=3,col="red")

DP.a.lp <- coef(DP.fit.lp)[1]
DP.b.lp <- coef(DP.fit.lp)[2]
DP.c.lp <- coef(DP.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
DP.n.mat <- matrix(0, nrow=DP.age.max+1, ncol=(t+1))
DP.n.mat[,1] <- DP.init.vec
DP.popmat <- DP.popmat.orig

## set up projection loop
for (i in 1:t) {
  DP.totN.i <- sum(DP.n.mat[,i])
  DP.pred.red <- as.numeric(DP.a.lp/(1+(DP.totN.i/DP.b.lp)^DP.c.lp))
  diag(DP.popmat[2:(DP.age.max+1),]) <- (DP.Sx[-(DP.age.max+1)])*DP.pred.red
  DP.popmat[DP.age.max+1,DP.age.max+1] <- (DP.Sx[DP.age.max+1])*DP.pred.red
  DP.popmat[1,] <- DP.pred.p.mm
  DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]
}

DP.n.pred <- colSums(DP.n.mat)
plot(yrs, DP.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=DP.pop.found, lty=2, col="red", lwd=2)


## stochastic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

DP.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DP.s.arr <- DP.m.arr <- array(data=NA, dim=c(t+1, DP.age.max+1, iter))
DP.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  DP.popmat <- DP.popmat.orig
  
  DP.n.mat <- DP.s.mat <- DP.m.mat <- matrix(0, nrow=DP.age.max+1,ncol=(t+1))
  DP.n.mat[,1] <- DP.init.vec
  DP.s.mat[,1] <- DP.Sx
  DP.m.mat[,1] <- DP.pred.p.mm
  DP.pred.red.vec <- rep(0, t)
  
  for (i in 1:t) {
    # stochastic survival values
    DP.s.alpha <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$alpha
    DP.s.beta <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$beta
    DP.s.stoch <- rbeta(length(DP.s.alpha), DP.s.alpha, DP.s.beta)
    
    # stochastic fertility sampler (Gaussian)
    DP.fert.stch <- rnorm(length(DP.popmat[,1]), DP.pred.p.mm, DP.m.sd.vec)
    DP.m.arr[i,,e] <- ifelse(DP.fert.stch < 0, 0, DP.fert.stch)
    
    DP.totN.i <- sum(DP.n.mat[,i], na.rm=T)
    DP.pred.red.vec[i] <- DP.a.lp/(1+(DP.totN.i/DP.b.lp)^DP.c.lp)
    
    diag(DP.popmat[2:(DP.age.max+1),]) <- (DP.s.stoch[-(DP.age.max+1)])*DP.pred.red.vec[i]
    DP.popmat[DP.age.max+1,DP.age.max+1] <- (DP.s.stoch[DP.age.max+1])*DP.pred.red.vec[i]
    DP.popmat[1,] <- DP.m.arr[i,,e]
    DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]

    DP.s.arr[i,,e] <- DP.s.stoch * DP.pred.red.vec[i]
    
  } # end i loop
  
  DP.n.sums.mat[e,] <- ((as.vector(colSums(DP.n.mat))/DP.pop.found))
  DP.pred.red.mn[e] <- mean(DP.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

DP.n.md <- apply(DP.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
DP.n.up <- apply(DP.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DP.n.lo <- apply(DP.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

save.image("DPvSsd.RData")


###########################
## Stable no catastrophe ##
###########################

## DIPROTODON (DP)
load("DPvSsd.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*DP.gen.l, 0) - 1)

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
    
    # stochastic fertility sampler (Gaussian)
    DP.fert.stch <- rnorm(length(DP.popmat[,1]), DP.pred.p.mm, DP.m.sd.vec)
    
    DP.totN.i <- sum(DP.n.mat[,i], na.rm=T)
    DP.pred.red.vec[i] <- DP.a.lp/(1+(DP.totN.i/DP.b.lp)^DP.c.lp)
    
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

save.image("DPvSsdstablenocat.RData")
rm(list = ls())




###################################
## Stable 40G (with catastrophe) ##
###################################

## DIPROTODON (DP)
load("DPvSsd.RData")
iter <- 10000
itdiv <- iter/10
t <- (round(40*DP.gen.l, 0) - 1)

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
    DP.pred.red.vec[i] <- DP.a.lp/(1+(DP.totN.i/DP.b.lp)^DP.c.lp)
    
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

save.image("DPvSsdstable40G.RData")
rm(list = ls())




#############################################
## HARVEST TO CAUSE DECLINE by r ~ -0.01   ##
## i.e., median trajectory hits r ~ -0.01  ##
## and median N goes ~ extinct by 40 G     ##
#############################################

  ## DIPROTODON (DP)
  load("DPvSsd.RData")
  harv.lim <- 1.04 # 0.95 gives target ~ -0.001
  iter <- 10000
  itdiv <- iter/10
  prop.red <- 0.0153 * harv.lim
  t <- (round(40*DP.gen.l, 0) - 1)
  
  ## set storage matrices & vectors
  DP.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  DP.pred.red.mn <- rep(0, iter)
  harvest <- round(prop.red*DP.init.vec, 0)
  sum(harvest)
  
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
      DP.pred.red.vec[i] <- DP.a.lp/(1+(DP.totN.i/DP.b.lp)^DP.c.lp)
      
      diag(DP.popmat[2:(DP.age.max+1),]) <- (DP.s.stoch[-(DP.age.max+1)])*DP.pred.red.vec[i]
      DP.popmat[DP.age.max+1,DP.age.max+1] <- (DP.s.stoch[DP.age.max+1])*DP.pred.red.vec[i]
      DP.popmat[1,] <- ifelse(DP.fert.stch < 0, 0, DP.fert.stch)
      DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]
      DP.n.mat[,i+1] <- ifelse((DP.n.mat[,i+1] - harvest) < 0, 0, (DP.n.mat[,i+1] - harvest)) # proportional harvest
      
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
  #save.image("DPvSsddeclSm.RData")
  save.image("DPvSsddecl.RData")
  rm(list = ls())
  
  
  
  
  #####################################################
  ## CARRYING CAPACITY DECLINES TO CAUSE r ~ -0.01   ##
  ## i.e., median trajectory hits r ~ -0.01          ##
  ## and median N goes ~ extinct by 40 G             ##
  #####################################################
  
  ## DIPROTODON (DP)
  load("DPvSsd.RData")
  iter <- 10000
  itdiv <- iter/10
  t <- (round(40*DP.gen.l, 0) - 1)
  divK <- 2.1 # 2.1 for r ~ -0.001
  DP.b.lp.vec1 <- seq(DP.b.lp, DP.b.lp/divK, by = -((DP.b.lp - DP.b.lp/divK)/(t-1)))
  DP.b.lp.vec2 <- rnorm(n=length(DP.b.lp.vec1), mean=DP.b.lp.vec1, sd=0.05*DP.b.lp.vec1[1])
    #DP.b.lp.vec2 <- rnorm(n=t, mean=DP.b.lp, sd=0.05*DP.b.lp)
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
  save.image("DPvSsddeclKsm.RData")
  rm(list = ls())
  
   
  
  
  #################################################################
  ## CARRYING CAPACITY STABLE, STOCHASTICALLY RESAMPLED,         ##
  ## CONSTANT VARIANCE, OR VARIANCE DOUBLES BY END OF PROJECTION ##
  #################################################################
  
  ## DIPROTODON (DP)
  load("DPvSsd.RData")
  iter <- 10000
  itdiv <- iter/10
  t <- (round(40*DP.gen.l, 0) - 1)
  divK <- 2.1 # 1800 # 2.1 for r ~ -0.001
  #DP.b.lp.vec1 <- seq(DP.b.lp, DP.b.lp/divK, by = -((DP.b.lp - DP.b.lp/divK)/(t-1)))
  #DP.b.lp.vec2 <- rnorm(n=length(DP.b.lp.vec1), mean=DP.b.lp.vec1, sd=0.05*DP.b.lp.vec1[1]) # increasing variance
  DP.b.lp.vec2 <- rnorm(n=t, mean=DP.b.lp, sd=0.05*DP.b.lp) # stable variance
  
  SDpc <- 0.05
  SD.incr <- seq(from=SDpc, to=2*SDpc, by=(2*SDpc-SDpc)/(t-1))
  DP.b.lp.vec2 <- rep(NA,t)
  for (f in 1:t) {
    DP.b.lp.vec2[f] <- rnorm(1, mean=DP.b.lp, sd=SD.incr[f]*DP.b.lp)
  }
  
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
  
  #save.image("DPvSsdKstochVarInc.RData")
  save.image("DPvSsdKstoch.RData")
  rm(list = ls())
  
  


  
  ############################################################
  ## NO DENSITY FEEDBACK; CATASTROPHE PROBABLITY SUFFICIENT ##
  ## TO CAUSE APPROXIMATE STABILITY OVER 40 G               ##
  ############################################################
  
  ## DIPROTODON (DP)
  load("DPvSsd.RData")
  iter <- 10000
  itdiv <- iter/10
  t <- (round(40*DP.gen.l, 0) - 1)
  catProbInc <- 18
  cat.prob <- ifelse(((0.14/DP.gen.l) * catProbInc) > 0.99, 0.99, ((0.14/DP.gen.l) * catProbInc))
  cat.prob
  catMag <- 0.5
  
  ## set storage matrices & vectors
  DP.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))

  for (e in 1:iter) {
    DP.popmat <- DP.popmat.orig
    
    DP.n.mat <- matrix(0, nrow=DP.age.max+1,ncol=(t+1))
    DP.n.mat[,1] <- DP.init.vec

    for (i in 1:t) {
      # stochastic survival values
      DP.s.alpha <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$alpha
      DP.s.beta <- estBetaParams(DP.Sx, DP.s.sd.vec^2)$beta
      DP.s.stoch <- rbeta(length(DP.s.alpha), DP.s.alpha, DP.s.beta)
      
      if (rbinom(1, 1, cat.prob) == 1) { # catastrophe
        cat.alpha <- estBetaParams(catMag, 0.05^2)$alpha
        cat.beta <- estBetaParams(catMag, 0.05^2)$beta
        DP.s.stoch <- DP.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertility sampler (Gaussian)
      DP.fert.stch <- rnorm(length(DP.popmat[,1]), DP.pred.p.mm, DP.m.sd.vec)
      
      DP.totN.i <- sum(DP.n.mat[,i], na.rm=T)

      diag(DP.popmat[2:(DP.age.max+1),]) <- (DP.s.stoch[-(DP.age.max+1)])
      DP.popmat[DP.age.max+1,DP.age.max+1] <- (DP.s.stoch[DP.age.max+1])
      DP.popmat[1,] <- ifelse(DP.fert.stch < 0, 0, DP.fert.stch)
      DP.n.mat[,i+1] <- DP.popmat %*% DP.n.mat[,i]
      
    } # end i loop
    
    DP.n.sums.mat[e,] <- (as.vector(colSums(DP.n.mat)))

    if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  # remove first 2 generations as burn-in
  DP.n.sums.matb <- DP.n.sums.mat[, -(1:round(DP.gen.l*2, 0))]
  
  # total N
  DP.n.md <- apply(DP.n.sums.matb, MARGIN=2, mean, na.rm=T) # mean over all iterations
  DP.n.up <- apply(DP.n.sums.matb, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  DP.n.lo <- apply(DP.n.sums.matb, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # plot
  yrs <- 1:(dim(DP.n.sums.matb)[2])
  plot(yrs, DP.n.md, type="l", lty=1, lwd=0.8)

  # mean rate of decline
  rmed1 <- na.omit(log(DP.n.md[2:length(yrs)] / DP.n.md[1:(length(yrs)-1)]))
  which.inf <- which(is.infinite(rmed1) == T)
  if (length(which.inf) == 0) {
    mean(rmed1)
  } else {
    mean(rmed1[-which.inf])
  }

  save.image("DPvSsdstableCat.RData")
  rm(list = ls())
  


    
  ################################################################
  ## ESTIMATE PHENOMENOLOGICAL GOMPERTZ DENSITY-FEEDBACK SIGNAL ##
  ## (choose which scenario to load)
  ################################################################
  load("DPvSsdstablenocat.RData") # stable, no catastrophe
  load("DPvSsdstable40G.RData") # stable, with catastrophe
  load("DPvSsddeclSm.RData") # r = 0.001
  load("DPvSsddecl.RData") # r = -0.01
  load("DPvSsdKstoch.RData") # stochastic K (K constant)
  load("DPvSsdKstochVarInc.RData") # stochastic K (K variance increasing)
  load("DPvSsddeclKsm.RData") # stochastic K, but declining at -0.001
  load("DPvSsdstableCat.RData") # no density feedback; stable from catastrophes
  
  
  ## return Gompertz signal
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
  DP.out
  
  
  
