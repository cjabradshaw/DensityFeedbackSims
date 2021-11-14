#############################################################################################################################
## density-feedback simulations
## base Leslie matrix models (no catastrophe)
## Corey Bradshaw & Salvador Herrando-Pérez
## Flinders University & Museo Nacional de Ciencias Naturales
#############################################################################################################################

# libraries
library(plotly)
library(ggpubr)

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

# MACROPOD CORRECTIONS (see http://doi.org/10.7554/eLife.63870 for justification)
# mass
NR.mass <- 5.1 # Notamacropus rufogriseus
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
OR.IBI <- 241
DOL.IBI <- 191
LC.IBI <- 153
LH.IBI <- 125
LF.IBI <- 365
MAg.IBI <- 220
MAn.IBI <- 270
MDo.IBI <- 211
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
OR.alpha <- mean(c(913,613))
DOL.alpha <- 450
LC.alpha <- 363
LH.alpha <- mean(c(315,345))
LF.alpha <- 365
MAg.alpha <- mean(c(405,357))
MAn.alpha <- 802
MDo.alpha <- 420
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
OR.alpha.pred <- (exp(-1.34 + (0.214*log(OR.mass*1000))))
DOL.alpha.pred <- (exp(-1.34 + (0.214*log(DOL.mass*1000))))
LC.alpha.pred <- (exp(-1.34 + (0.214*log(LC.mass*1000))))
LH.alpha.pred <- (exp(-1.34 + (0.214*log(LH.mass*1000))))
LF.alpha.pred <- (exp(-1.34 + (0.214*log(LF.mass*1000))))
MAg.alpha.pred <- (exp(-1.34 + (0.214*log(MAg.mass*1000))))
MAn.alpha.pred <- (exp(-1.34 + (0.214*log(MAn.mass*1000))))
MDo.alpha.pred <- (exp(-1.34 + (0.214*log(MDo.mass*1000))))
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


# VOMBATIFORM CORRECTIONS  (see http://doi.org/10.7554/eLife.63870 for justification)
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
DP.age.max1 <- round(10^(0.89 + (0.13*log10(DP.mass*1000))), 0)
DP.age.max <- round(DP.age.max1 * 26/29, 0) # corrected for over-estimate derived from Vombatus

## age vector
DP.age.vec <- 0:DP.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
DP.F.pred1 <- exp(2.719 - (0.211*log(DP.mass*1000)))/2 # divided by 2 for females
DP.F.pred1
DP.F.pred <- DP.F.pred1 * (365/as.numeric((mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(DP.mass))))
DP.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
DP.alpha1 <- ceiling(exp(-1.34 + (0.214*log(DP.mass*1000))))
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
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
DP.Sx <- c(0.99*DP.s.ad.yr, 1 - qx)
plot(x, DP.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
DP.s.sd.vec <- 0.05*DP.Sx

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





############################
## PALORCHESTES (azael) (PA)
## sources: Brook & Johnson 2006 (Alcheringa 30:39-48, http://doi.org/10.1080/03115510609506854)
##          Richards et al. 2019 (PLoS One https://doi.org/10.1371/journal.pone.0221824)

# mass
PA.mass <- 1000 # kg (Richards et al. 2019 (PLoS One https://doi.org/10.1371/journal.pone.0221824))

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
PA.rm.pred <- 10^(0.6914 - (0.2622*log10(PA.mass*1000)))
PA.rm.pred
PA.lm.pred <- exp(PA.rm.pred)
PA.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
PA.D.pred <- (10^(4.196 - (0.74*log10(PA.mass*1000))))/2 # divided by 2 for females only
PA.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
PA.age.max1 <- round(10^(0.89 + (0.13*log10(PA.mass*1000))), 0)
PA.age.max1
PA.age.max <- round(PA.age.max1 * 26/29, 0) # corrected for over-estimate derived from Vombatus
PA.age.max

## age vector
PA.age.vec <- 0:PA.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
PA.F.pred1 <- exp(2.719 - (0.211*log(PA.mass*1000)))/2 # divided by 2 for females
PA.F.pred1
PA.F.pred <- PA.F.pred1 * (365/as.numeric((mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(PA.mass))))
PA.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
PA.alpha1 <- ceiling(exp(-1.34 + (0.214*log(PA.mass*1000))))
PA.alpha1
PA.alpha <- PA.alpha1
PA.alpha

## define m function with age
PA.m.vec <- c(rep(0, PA.alpha-1), rep(0.75*PA.F.pred, round(PA.alpha/2,0)), rep(PA.F.pred, (PA.age.max+1-((PA.alpha-1+round(PA.alpha/2,0))))))
PA.m.sd.vec <- 0.05*PA.m.vec
plot(PA.age.vec, PA.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
PA.m.dat <- data.frame(PA.age.vec, PA.m.vec)
param.init <- c(0.3, 7, -5)
PA.fit.logp <- nls(PA.m.vec ~ a / (1+(PA.age.vec/b)^c), 
                   data = PA.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PA.fit.logp.summ <- summary(PA.fit.logp)
plot(PA.age.vec, PA.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
PA.age.vec.cont <- seq(0,max(PA.age.vec),1)
PA.pred.p.m <- coef(PA.fit.logp)[1] / (1+(PA.age.vec.cont/coef(PA.fit.logp)[2])^coef(PA.fit.logp)[3])
PA.pred.p.mm <- ifelse(PA.pred.p.m > 1, 1, PA.pred.p.m)
lines(PA.age.vec.cont, PA.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
PA.s.tran <- ln.a.s + b.s*log(PA.mass*1000) + log(1)
PA.s.ad.yr <- exp(-exp(PA.s.tran))
PA.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (PA.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.5 # rate of mortality decline (also known as bt)
a2 <- 1 - PA.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.7e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.05 # rate of mortality increase
longev <- PA.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
PA.Sx <- c(0.99*PA.s.ad.yr, 1 - qx)
plot(x, PA.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
PA.s.sd.vec <- 0.05*PA.Sx

## create matrix
PA.popmat <- matrix(data = 0, nrow=PA.age.max+1, ncol=PA.age.max+1)
diag(PA.popmat[2:(PA.age.max+1),]) <- PA.Sx[-(PA.age.max+1)]
PA.popmat[PA.age.max+1,PA.age.max+1] <- PA.Sx[PA.age.max+1]
PA.popmat[1,] <- PA.pred.p.mm
colnames(PA.popmat) <- c(0:PA.age.max)
rownames(PA.popmat) <- c(0:PA.age.max)
PA.popmat.orig <- PA.popmat ## save original matrix

## matrix properties
max.lambda(PA.popmat.orig) ## 1-yr lambda
PA.lm.pred
max.r(PA.popmat.orig) # rate of population change, 1-yr
PA.ssd <- stable.stage.dist(PA.popmat.orig) ## stable stage distribution
plot(PA.age.vec, PA.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(PA.popmat.orig, PA.age.max) # reproductive value
PA.gen.l <- G.val(PA.popmat.orig, PA.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
PA.pop.found <- round(area*PA.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
PA.init.vec <- PA.ssd * PA.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*PA.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

PA.tot.F <- sum(PA.popmat.orig[1,])
PA.popmat <- PA.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
PA.n.mat <- matrix(0, nrow=PA.age.max+1,ncol=(t+1))
PA.n.mat[,1] <- PA.init.vec

## set up projection loop
for (i in 1:t) {
  PA.n.mat[,i+1] <- PA.popmat %*% PA.n.mat[,i]
}

PA.n.pred <- colSums(PA.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(PA.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
PA.K.max <- 1*PA.pop.found
PA.K.vec <- c(1, PA.K.max/2, 0.75*PA.K.max, PA.K.max) 
PA.red.vec <- c(1,0.982,0.957,0.9206)
plot(PA.K.vec, PA.red.vec,pch=19,type="b")
PA.Kred.dat <- data.frame(PA.K.vec, PA.red.vec)

# logistic power function a/(1+(x/b)^c)
PA.param.init <- c(1, 2*PA.K.max, 2)
PA.fit.lp <- nls(PA.red.vec ~ a/(1+(PA.K.vec/b)^c), 
                 data = PA.Kred.dat,
                 algorithm = "port",
                 start = c(a = PA.param.init[1], b = PA.param.init[2], c = PA.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PA.fit.lp.summ <- summary(PA.fit.lp)
plot(PA.K.vec, PA.red.vec, pch=19,xlab="N",ylab="reduction factor")
PA.K.vec.cont <- seq(1,2*PA.pop.found,1)
PA.pred.lp.fx <- coef(PA.fit.lp)[1]/(1+(PA.K.vec.cont/coef(PA.fit.lp)[2])^coef(PA.fit.lp)[3])
lines(PA.K.vec.cont, PA.pred.lp.fx, lty=3,lwd=3,col="red")

PA.a.lp <- coef(PA.fit.lp)[1]
PA.b.lp <- coef(PA.fit.lp)[2]
PA.c.lp <- coef(PA.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
PA.n.mat <- matrix(0, nrow=PA.age.max+1, ncol=(t+1))
PA.n.mat[,1] <- PA.init.vec
PA.popmat <- PA.popmat.orig

## set up projection loop
for (i in 1:t) {
  PA.totN.i <- sum(PA.n.mat[,i])
  PA.pred.red <- as.numeric(PA.a.lp/(1+(PA.totN.i/PA.b.lp)^PA.c.lp))
  diag(PA.popmat[2:(PA.age.max+1),]) <- (PA.Sx[-(PA.age.max+1)])*PA.pred.red
  PA.popmat[PA.age.max+1,PA.age.max+1] <- (PA.Sx[PA.age.max+1])*PA.pred.red
  PA.popmat[1,] <- PA.pred.p.mm
  PA.n.mat[,i+1] <- PA.popmat %*% PA.n.mat[,i]
}

PA.n.pred <- colSums(PA.n.mat)
plot(yrs, PA.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=PA.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

PA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PA.s.arr <- PA.m.arr <- array(data=NA, dim=c(t+1, PA.age.max+1, iter))
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
    
    # stochastic fertilty sampler (gaussian)
    PA.fert.stch <- rnorm(length(PA.popmat[,1]), PA.pred.p.mm, PA.m.sd.vec)
    PA.m.arr[i,,e] <- ifelse(PA.fert.stch < 0, 0, PA.fert.stch)
    
    PA.totN.i <- sum(PA.n.mat[,i], na.rm=T)
    PA.pred.red.vec[i] <- PA.a.lp/(1+(PA.totN.i/PA.b.lp)^PA.c.lp)
    
    diag(PA.popmat[2:(PA.age.max+1),]) <- (PA.s.stoch[-(PA.age.max+1)])*PA.pred.red.vec[i]
    PA.popmat[PA.age.max+1,PA.age.max+1] <- (PA.s.stoch[PA.age.max+1])*PA.pred.red.vec[i]
    PA.popmat[1,] <- PA.m.arr[i,,e]
    PA.n.mat[,i+1] <- PA.popmat %*% PA.n.mat[,i]

    PA.s.arr[i,,e] <- PA.s.stoch * PA.pred.red.vec[i]
    
  } # end i loop
  
  PA.n.sums.mat[e,] <- ((as.vector(colSums(PA.n.mat))/PA.pop.found))
  PA.pred.red.mn[e] <- mean(PA.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

PA.n.md <- apply(PA.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
PA.n.up <- apply(PA.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PA.n.lo <- apply(PA.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations




##############################
## ZYGOMATURUS (trilobus) (ZT)
## sources: Brook & Johnson 2006 (Alcheringa 30:39-48, http://doi.org/10.1080/03115510609506854)
##          Richards et al. 2019 (PLoS One https://doi.org/10.1371/journal.pone.0221824)

# mass
ZT.mass <- 500 # kg (Johnson et al. 2006 Australia’s Mammal Extinctions. 278 pp. Cambridge  University Press, Melbourne)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
ZT.rm.pred <- 10^(0.6914 - (0.2622*log10(ZT.mass*1000)))
ZT.rm.pred
ZT.lm.pred <- exp(ZT.rm.pred)
ZT.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
ZT.D.pred <- (10^(4.196 - (0.74*log10(ZT.mass*1000))))/2 # divided by 2 for females only
ZT.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
ZT.age.max1 <- round(10^(0.89 + (0.13*log10(ZT.mass*1000))), 0)
ZT.age.max1
ZT.age.max <- round(ZT.age.max1 * 26/29, 0) # corrected for over-estimate derived from Vombatus
ZT.age.max

## age vector
ZT.age.vec <- 0:ZT.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
ZT.F.pred1 <- exp(2.719 - (0.211*log(ZT.mass*1000)))/2 # divided by 2 for females
ZT.F.pred1
ZT.F.pred <- ZT.F.pred1 * (365/as.numeric((mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(ZT.mass))))
ZT.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
ZT.alpha1 <- ceiling(exp(-1.34 + (0.214*log(ZT.mass*1000))))
ZT.alpha1
ZT.alpha <- ZT.alpha1
ZT.alpha

## define m function with age
ZT.m.vec <- c(rep(0, ZT.alpha-1), rep(0.75*ZT.F.pred, round(ZT.alpha/2,0)), rep(ZT.F.pred, (ZT.age.max+1-((ZT.alpha-1+round(ZT.alpha/2,0))))))
ZT.m.sd.vec <- 0.05*ZT.m.vec
plot(ZT.age.vec, ZT.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
ZT.m.dat <- data.frame(ZT.age.vec, ZT.m.vec)
param.init <- c(0.1, 2, -5)
ZT.fit.logp <- nls(ZT.m.vec ~ a / (1+(ZT.age.vec/b)^c), 
                   data = ZT.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
ZT.fit.logp.summ <- summary(ZT.fit.logp)
plot(ZT.age.vec, ZT.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
ZT.age.vec.cont <- seq(0,max(ZT.age.vec),1)
ZT.pred.p.m <- coef(ZT.fit.logp)[1] / (1+(ZT.age.vec.cont/coef(ZT.fit.logp)[2])^coef(ZT.fit.logp)[3])
ZT.pred.p.mm <- ifelse(ZT.pred.p.m > 1, 1, ZT.pred.p.m)
lines(ZT.age.vec.cont, ZT.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
ZT.s.tran <- ln.a.s + b.s*log(ZT.mass*1000) + log(1)
ZT.s.ad.yr <- exp(-exp(ZT.s.tran))
ZT.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.01*ZT.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.7 # rate of mortality decline (also known as bt)
a2 <- 1 - ZT.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.02 # rate of mortality increase
longev <- ZT.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
ZT.Sx <- c(0.995*ZT.s.ad.yr, 1 - qx)
plot(x, ZT.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
ZT.s.sd.vec <- 0.05*ZT.Sx

## create matrix
ZT.popmat <- matrix(data = 0, nrow=ZT.age.max+1, ncol=ZT.age.max+1)
diag(ZT.popmat[2:(ZT.age.max+1),]) <- ZT.Sx[-(ZT.age.max+1)]
ZT.popmat[ZT.age.max+1,ZT.age.max+1] <- ZT.Sx[ZT.age.max+1]
ZT.popmat[1,] <- ZT.pred.p.mm
colnames(ZT.popmat) <- c(0:ZT.age.max)
rownames(ZT.popmat) <- c(0:ZT.age.max)
ZT.popmat.orig <- ZT.popmat ## save original matrix

## matrix properties
max.lambda(ZT.popmat.orig) ## 1-yr lambda
ZT.lm.pred
max.r(ZT.popmat.orig) # rate of population change, 1-yr
ZT.ssd <- stable.stage.dist(ZT.popmat.orig) ## stable stage distribution
plot(ZT.age.vec, ZT.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(ZT.popmat.orig, ZT.age.max) # reproductive value
ZT.gen.l <- G.val(ZT.popmat.orig, ZT.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
ZT.pop.found <- round(area*ZT.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
ZT.init.vec <- ZT.ssd * ZT.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*ZT.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

ZT.tot.F <- sum(ZT.popmat.orig[1,])
ZT.popmat <- ZT.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
ZT.n.mat <- matrix(0, nrow=ZT.age.max+1,ncol=(t+1))
ZT.n.mat[,1] <- ZT.init.vec

## set up projection loop
for (i in 1:t) {
  ZT.n.mat[,i+1] <- ZT.popmat %*% ZT.n.mat[,i]
}

ZT.n.pred <- colSums(ZT.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(ZT.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
ZT.K.max <- 1*ZT.pop.found
ZT.K.vec <- c(1, ZT.K.max/2, 0.75*ZT.K.max, ZT.K.max) 
ZT.red.vec <- c(1,0.975,0.945,0.9018)
plot(ZT.K.vec, ZT.red.vec,pch=19,type="b")
ZT.Kred.dat <- data.frame(ZT.K.vec, ZT.red.vec)

# logistic power function a/(1+(x/b)^c)
ZT.param.init <- c(1, 2*ZT.K.max, 2)
ZT.fit.lp <- nls(ZT.red.vec ~ a/(1+(ZT.K.vec/b)^c), 
                 data = ZT.Kred.dat,
                 algorithm = "port",
                 start = c(a = ZT.param.init[1], b = ZT.param.init[2], c = ZT.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
ZT.fit.lp.summ <- summary(ZT.fit.lp)
plot(ZT.K.vec, ZT.red.vec, pch=19,xlab="N",ylab="reduction factor")
ZT.K.vec.cont <- seq(1,2*ZT.pop.found,1)
ZT.pred.lp.fx <- coef(ZT.fit.lp)[1]/(1+(ZT.K.vec.cont/coef(ZT.fit.lp)[2])^coef(ZT.fit.lp)[3])
lines(ZT.K.vec.cont, ZT.pred.lp.fx, lty=3,lwd=3,col="red")

ZT.a.lp <- coef(ZT.fit.lp)[1]
ZT.b.lp <- coef(ZT.fit.lp)[2]
ZT.c.lp <- coef(ZT.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
ZT.n.mat <- matrix(0, nrow=ZT.age.max+1, ncol=(t+1))
ZT.n.mat[,1] <- ZT.init.vec
ZT.popmat <- ZT.popmat.orig

## set up projection loop
for (i in 1:t) {
  ZT.totN.i <- sum(ZT.n.mat[,i])
  ZT.pred.red <- as.numeric(ZT.a.lp/(1+(ZT.totN.i/ZT.b.lp)^ZT.c.lp))
  diag(ZT.popmat[2:(ZT.age.max+1),]) <- (ZT.Sx[-(ZT.age.max+1)])*ZT.pred.red
  ZT.popmat[ZT.age.max+1,ZT.age.max+1] <- (ZT.Sx[ZT.age.max+1])*ZT.pred.red
  ZT.popmat[1,] <- ZT.pred.p.mm
  ZT.n.mat[,i+1] <- ZT.popmat %*% ZT.n.mat[,i]
}

ZT.n.pred <- colSums(ZT.n.mat)
plot(yrs, ZT.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=ZT.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

ZT.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
ZT.s.arr <- ZT.m.arr <- array(data=NA, dim=c(t+1, ZT.age.max+1, iter))
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
    
    # stochastic fertilty sampler (gaussian)
    ZT.fert.stch <- rnorm(length(ZT.popmat[,1]), ZT.pred.p.mm, ZT.m.sd.vec)
    ZT.m.arr[i,,e] <- ifelse(ZT.fert.stch < 0, 0, ZT.fert.stch)
    
    ZT.totN.i <- sum(ZT.n.mat[,i], na.rm=T)
    ZT.pred.red.vec[i] <- ZT.a.lp/(1+(ZT.totN.i/ZT.b.lp)^ZT.c.lp)
    
    diag(ZT.popmat[2:(ZT.age.max+1),]) <- (ZT.s.stoch[-(ZT.age.max+1)])*ZT.pred.red.vec[i]
    ZT.popmat[ZT.age.max+1,ZT.age.max+1] <- (ZT.s.stoch[ZT.age.max+1])*ZT.pred.red.vec[i]
    ZT.popmat[1,] <- ZT.m.arr[i,,e]
    ZT.n.mat[,i+1] <- ZT.popmat %*% ZT.n.mat[,i]

    ZT.s.arr[i,,e] <- ZT.s.stoch * ZT.pred.red.vec[i]
    
  } # end i loop
  
  ZT.n.sums.mat[e,] <- ((as.vector(colSums(ZT.n.mat))/ZT.pop.found))
  ZT.pred.red.mn[e] <- mean(ZT.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

ZT.n.md <- apply(ZT.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
ZT.n.up <- apply(ZT.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
ZT.n.lo <- apply(ZT.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations





##############################
## PHASCOLONUS (gigas) (PH)

# mass
PH.mass <- 200 # Phascolonus gigas (Johnson et al. 2006)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
PH.rm.pred <- 10^(0.6914 - (0.2622*log10(PH.mass*1000)))
PH.rm.pred
PH.lm.pred <- exp(PH.rm.pred)
PH.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
PH.D.pred <- (10^(4.196 - (0.74*log10(PH.mass*1000))))/2 # divided by 2 for females only
PH.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
PH.age.max1 <- round(10^(0.89 + (0.13*log10(PH.mass*1000))), 0)
PH.age.max1
PH.age.max <- round(PH.age.max1 * 26/29, 0) # corrected for over-estimate derived from Vombatus
PH.age.max

## age vector
PH.age.vec <- 0:PH.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
PH.F.pred1 <- exp(2.719 - (0.211*log(PH.mass*1000)))/2 # divided by 2 for females
PH.F.pred1
PH.F.pred <- PH.F.pred1 * (365/as.numeric((mean(VOMBAT.IBI) + (VU.IBI - VU.IBI.pred)) + (MACROPOD.F.corr.b * log10(PH.mass))))
PH.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
PH.alpha1 <- ceiling(exp(-1.34 + (0.214*log(PH.mass*1000))))
PH.alpha1
PH.alpha <- PH.alpha1
PH.alpha

## define m function with age
PH.m.vec <- c(rep(0, PH.alpha-1), rep(0.75*PH.F.pred, round(PH.alpha/2,0)), rep(PH.F.pred, (PH.age.max+1-((PH.alpha-1+round(PH.alpha/2,0))))))
PH.m.sd.vec <- 0.05*PH.m.vec
plot(PH.age.vec, PH.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
PH.m.dat <- data.frame(PH.age.vec, PH.m.vec)
param.init <- c(0.5, 2, -4)
PH.fit.logp <- nls(PH.m.vec ~ a / (1+(PH.age.vec/b)^c), 
                   data = PH.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PH.fit.logp.summ <- summary(PH.fit.logp)
plot(PH.age.vec, PH.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
PH.age.vec.cont <- seq(0,max(PH.age.vec),1)
PH.pred.p.m <- coef(PH.fit.logp)[1] / (1+(PH.age.vec.cont/coef(PH.fit.logp)[2])^coef(PH.fit.logp)[3])
PH.pred.p.mm <- ifelse(PH.pred.p.m > 1, 1, PH.pred.p.m)
lines(PH.age.vec.cont, PH.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
PH.s.tran <- ln.a.s + b.s*log(PH.mass*1000) + log(1)
PH.s.ad.yr <- exp(-exp(PH.s.tran))
PH.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.01*PH.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 3.1 # rate of mortality decline (also known as bt)
a2 <- 1 - PH.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.02 # rate of mortality increase
longev <- PH.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
PH.Sx <- c(0.99*PH.s.ad.yr, 1 - qx)
plot(x, PH.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
PH.s.sd.vec <- 0.05*PH.Sx

## create matrix
PH.popmat <- matrix(data = 0, nrow=PH.age.max+1, ncol=PH.age.max+1)
diag(PH.popmat[2:(PH.age.max+1),]) <- PH.Sx[-(PH.age.max+1)]
PH.popmat[PH.age.max+1,PH.age.max+1] <- PH.Sx[PH.age.max+1]
PH.popmat[1,] <- PH.pred.p.mm
colnames(PH.popmat) <- c(0:PH.age.max)
rownames(PH.popmat) <- c(0:PH.age.max)
PH.popmat.orig <- PH.popmat ## save original matrix

## matrix properties
max.lambda(PH.popmat.orig) ## 1-yr lambda
PH.lm.pred
max.r(PH.popmat.orig) # rate of population change, 1-yr
PH.ssd <- stable.stage.dist(PH.popmat.orig) ## stable stage distribution
plot(PH.age.vec, PH.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(PH.popmat.orig, PH.age.max) # reproductive value
PH.gen.l <- G.val(PH.popmat.orig, PH.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
PH.pop.found <- round(area*PH.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
PH.init.vec <- PH.ssd * PH.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*PH.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

PH.tot.F <- sum(PH.popmat.orig[1,])
PH.popmat <- PH.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
PH.n.mat <- matrix(0, nrow=PH.age.max+1,ncol=(t+1))
PH.n.mat[,1] <- PH.init.vec

## set up projection loop
for (i in 1:t) {
  PH.n.mat[,i+1] <- PH.popmat %*% PH.n.mat[,i]
}

PH.n.pred <- colSums(PH.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(PH.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
PH.K.max <- 1*PH.pop.found
PH.K.vec <- c(1, PH.K.max/2, 0.75*PH.K.max, PH.K.max) 
PH.red.vec <- c(1,0.98,0.945,0.874)
plot(PH.K.vec, PH.red.vec,pch=19,type="b")
PH.Kred.dat <- data.frame(PH.K.vec, PH.red.vec)

# logistic power function a/(1+(x/b)^c)
PH.param.init <- c(1, 2*PH.K.max, 2)
PH.fit.lp <- nls(PH.red.vec ~ a/(1+(PH.K.vec/b)^c), 
                 data = PH.Kred.dat,
                 algorithm = "port",
                 start = c(a = PH.param.init[1], b = PH.param.init[2], c = PH.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PH.fit.lp.summ <- summary(PH.fit.lp)
plot(PH.K.vec, PH.red.vec, pch=19,xlab="N",ylab="reduction factor")
PH.K.vec.cont <- seq(1,2*PH.pop.found,1)
PH.pred.lp.fx <- coef(PH.fit.lp)[1]/(1+(PH.K.vec.cont/coef(PH.fit.lp)[2])^coef(PH.fit.lp)[3])
lines(PH.K.vec.cont, PH.pred.lp.fx, lty=3,lwd=3,col="red")

PH.a.lp <- coef(PH.fit.lp)[1]
PH.b.lp <- coef(PH.fit.lp)[2]
PH.c.lp <- coef(PH.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
PH.n.mat <- matrix(0, nrow=PH.age.max+1, ncol=(t+1))
PH.n.mat[,1] <- PH.init.vec
PH.popmat <- PH.popmat.orig

## set up projection loop
for (i in 1:t) {
  PH.totN.i <- sum(PH.n.mat[,i])
  PH.pred.red <- as.numeric(PH.a.lp/(1+(PH.totN.i/PH.b.lp)^PH.c.lp))
  diag(PH.popmat[2:(PH.age.max+1),]) <- (PH.Sx[-(PH.age.max+1)])*PH.pred.red
  PH.popmat[PH.age.max+1,PH.age.max+1] <- (PH.Sx[PH.age.max+1])*PH.pred.red
  PH.popmat[1,] <- PH.pred.p.mm
  PH.n.mat[,i+1] <- PH.popmat %*% PH.n.mat[,i]
}

PH.n.pred <- colSums(PH.n.mat)
plot(yrs, PH.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=PH.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

PH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PH.s.arr <- PH.m.arr <- array(data=NA, dim=c(t+1, PH.age.max+1, iter))
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
    
    # stochastic fertilty sampler (gaussian)
    PH.fert.stch <- rnorm(length(PH.popmat[,1]), PH.pred.p.mm, PH.m.sd.vec)
    PH.m.arr[i,,e] <- ifelse(PH.fert.stch < 0, 0, PH.fert.stch)
    
    PH.totN.i <- sum(PH.n.mat[,i], na.rm=T)
    PH.pred.red.vec[i] <- PH.a.lp/(1+(PH.totN.i/PH.b.lp)^PH.c.lp)
    
    diag(PH.popmat[2:(PH.age.max+1),]) <- (PH.s.stoch[-(PH.age.max+1)])*PH.pred.red.vec[i]
    PH.popmat[PH.age.max+1,PH.age.max+1] <- (PH.s.stoch[PH.age.max+1])*PH.pred.red.vec[i]
    PH.popmat[1,] <- PH.m.arr[i,,e]
    PH.n.mat[,i+1] <- PH.popmat %*% PH.n.mat[,i]

    PH.s.arr[i,,e] <- PH.s.stoch * PH.pred.red.vec[i]
    
  } # end i loop
  
  PH.n.sums.mat[e,] <- ((as.vector(colSums(PH.n.mat))/PH.pop.found))
  PH.pred.red.mn[e] <- mean(PH.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

PH.n.md <- apply(PH.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
PH.n.up <- apply(PH.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PH.n.lo <- apply(PH.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## VOMBATUS (ursinus) (VU)
## sources: Roger et al. 2011 Popul Ecol 53:215-227

# mass
VU.mass <- 25 # Vombatus ursinus (Saran et al. 2011 Pacific Conservation Biology 17:310-319)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
VU.rm.pred <- 10^(0.6914 - (0.2622*log10(VU.mass*1000)))
VU.rm.pred
VU.lm.pred <- exp(VU.rm.pred)
VU.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
VU.D.pred <- (10^(4.196 - (0.74*log10(VU.mass*1000))))/2 # divided by 2 for females only
VU.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
VU.age.max <- round(10^(0.89 + (0.13*log10(VU.mass*1000))), 0)
VU.age.max
VU.age.max <- 26 # reset McIlroy 2008

## age vector
VU.age.vec <- 0:VU.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
VU.F.pred <- exp(2.719 - (0.211*log(VU.mass*1000)))/2 # divided by 2 for females
VU.F.pred
VU.F.pred <- 1/2/2 # inter-birth interval of 2 years / 2 for daughters only (McIlroy 1995)

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
VU.alpha <- ceiling(exp(-1.34 + (0.214*log(VU.mass*1000))))
VU.alpha
VU.alpha <- 2  # reset according to Roger et al. (2011)

## define m function with age
## pouch young per year
VU.pypy <- 0.5

## proportion of females breeding each year
VU.pfbr <- 0.84

## proportion young that are female
VU.sr <- 0.5 

## m vector
#VU.m.vec <- c(0, 0.3*VU.pypy*VU.pfbr*VU.sr, 0.9*VU.pypy*VU.pfbr*VU.sr, rep(VU.pypy*VU.pfbr*VU.sr,24))
VU.m.vec <- c(0, 0.3*VU.F.pred, 0.9*VU.F.pred, rep(VU.F.pred, 24))
VU.m.sd.vec <- 0.05*VU.m.vec
plot(VU.age.vec, VU.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
VU.m.dat <- data.frame(VU.age.vec, VU.m.vec)
param.init <- c(0.2, 2, -3)
VU.fit.logp <- nls(VU.m.vec ~ a / (1+(VU.age.vec/b)^c), 
                   data = VU.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
VU.fit.logp.summ <- summary(VU.fit.logp)
plot(VU.age.vec, VU.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
VU.age.vec.cont <- seq(0,max(VU.age.vec),1)
VU.pred.p.m <- coef(VU.fit.logp)[1] / (1+(VU.age.vec.cont/coef(VU.fit.logp)[2])^coef(VU.fit.logp)[3])
VU.pred.p.mm <- ifelse(VU.pred.p.m > 1, 1, VU.pred.p.m)
VU.pred.p.mm[2] <- 0
lines(VU.age.vec.cont, VU.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
VU.s.tran <- ln.a.s + b.s*log(VU.mass*1000) + log(1)
VU.s.ad.yr <- exp(-exp(VU.s.tran))
VU.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.05*VU.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 3.1 # rate of mortality decline (also known as bt)
a2 <- 1 - VU.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.02 # rate of mortality increase
longev <- VU.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
VU.Sx <- c(0.99*VU.s.ad.yr, 1 - qx)
plot(x, VU.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
VU.s.sd.vec <- 0.05*VU.Sx

## create matrix
VU.popmat <- matrix(data = 0, nrow=VU.age.max+1, ncol=VU.age.max+1)
diag(VU.popmat[2:(VU.age.max+1),]) <- VU.Sx[-(VU.age.max+1)]
VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
VU.popmat[1,] <- VU.pred.p.mm
colnames(VU.popmat) <- c(0:VU.age.max)
rownames(VU.popmat) <- c(0:VU.age.max)
VU.popmat.orig <- VU.popmat ## save original matrix

## matrix properties
max.lambda(VU.popmat.orig) ## 1-yr lambda
VU.lm.pred
max.r(VU.popmat.orig) # rate of population change, 1-yr
VU.ssd <- stable.stage.dist(VU.popmat.orig) ## stable stage distribution
plot(VU.age.vec, VU.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(VU.popmat.orig, VU.age.max) # reproductive value
VU.gen.l <- G.val(VU.popmat.orig, VU.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
VU.pop.found <- round(area*VU.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
VU.init.vec <- VU.ssd * VU.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*VU.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

VU.tot.F <- sum(VU.popmat.orig[1,])
VU.popmat <- VU.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
VU.n.mat <- matrix(0, nrow=VU.age.max+1,ncol=(t+1))
VU.n.mat[,1] <- VU.init.vec

## set up projection loop
for (i in 1:t) {
  VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]
}

VU.n.pred <- colSums(VU.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(VU.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
VU.K.max <- 1*VU.pop.found
VU.K.vec <- c(1, VU.K.max/2, 0.75*VU.K.max, VU.K.max) 
VU.red.vec <- c(1,0.985,0.951,0.873)
plot(VU.K.vec, VU.red.vec,pch=19,type="b")
VU.Kred.dat <- data.frame(VU.K.vec, VU.red.vec)

# logistic power function a/(1+(x/b)^c)
VU.param.init <- c(1, 2*VU.K.max, 2)
VU.fit.lp <- nls(VU.red.vec ~ a/(1+(VU.K.vec/b)^c), 
                 data = VU.Kred.dat,
                 algorithm = "port",
                 start = c(a = VU.param.init[1], b = VU.param.init[2], c = VU.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
VU.fit.lp.summ <- summary(VU.fit.lp)
plot(VU.K.vec, VU.red.vec, pch=19,xlab="N",ylab="reduction factor")
VU.K.vec.cont <- seq(1,2*VU.pop.found,1)
VU.pred.lp.fx <- coef(VU.fit.lp)[1]/(1+(VU.K.vec.cont/coef(VU.fit.lp)[2])^coef(VU.fit.lp)[3])
lines(VU.K.vec.cont, VU.pred.lp.fx, lty=3,lwd=3,col="red")

VU.a.lp <- coef(VU.fit.lp)[1]
VU.b.lp <- coef(VU.fit.lp)[2]
VU.c.lp <- coef(VU.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
VU.n.mat <- matrix(0, nrow=VU.age.max+1, ncol=(t+1))
VU.n.mat[,1] <- VU.init.vec
VU.popmat <- VU.popmat.orig

## set up projection loop
for (i in 1:t) {
  VU.totN.i <- sum(VU.n.mat[,i])
  VU.pred.red <- as.numeric(VU.a.lp/(1+(VU.totN.i/VU.b.lp)^VU.c.lp))
  diag(VU.popmat[2:(VU.age.max+1),]) <- (VU.Sx[-(VU.age.max+1)])*VU.pred.red
  VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
  VU.popmat[1,] <- VU.pred.p.mm
  VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]
}

VU.n.pred <- colSums(VU.n.mat)
plot(yrs, VU.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=VU.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

VU.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
VU.s.arr <- VU.m.arr <- array(data=NA, dim=c(t+1, VU.age.max+1, iter))
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
    
    # stochastic fertilty sampler (gaussian)
    VU.fert.stch <- rnorm(length(VU.popmat[,1]), VU.pred.p.mm, VU.m.sd.vec)
    VU.m.arr[i,,e] <- ifelse(VU.fert.stch < 0, 0, VU.fert.stch)
    
    VU.totN.i <- sum(VU.n.mat[,i], na.rm=T)
    VU.pred.red.vec[i] <- VU.a.lp/(1+(VU.totN.i/VU.b.lp)^VU.c.lp)
    
    diag(VU.popmat[2:(VU.age.max+1),]) <- (VU.s.stoch[-(VU.age.max+1)])*VU.pred.red.vec[i]
    VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
    VU.popmat[1,] <- VU.m.arr[i,,e]
    VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]

    VU.s.arr[i,,e] <- VU.s.stoch * VU.pred.red.vec[i]
    
  } # end i loop
  
  VU.n.sums.mat[e,] <- ((as.vector(colSums(VU.n.mat))/VU.pop.found))
  VU.pred.red.mn[e] <- mean(VU.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

VU.n.md <- apply(VU.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
VU.n.up <- apply(VU.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
VU.n.lo <- apply(VU.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## PROCOPTODON (goliah) (PG)

## mass
PG.mass <- 250 # Procoptodon goliah (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1442-9993.2004.01389.x)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
PG.rm.pred <- 10^(0.6914 - (0.2622*log10(PG.mass*1000)))
PG.rm.pred
PG.lm.pred <- exp(PG.rm.pred)
PG.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
PG.D.pred <- (10^(4.196 - (0.74*log10(PG.mass*1000))))/2 # divided by 2 for females only
PG.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
PG.age.max1 <- round(10^(0.89 + (0.13*log10(PG.mass*1000))), 0)
PG.age.max1
PG.age.max <- ceiling(13/30*PG.age.max1) # (amount overestimated for OR from Jonzén et al. 2010-J Anim Ecol)
PG.age.max

## age vector
PG.age.vec <- 0:PG.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
PG.F.pred1 <- exp(2.719 - (0.211*log(PG.mass*1000)))/2 # divided by 2 for females
PG.F.pred1

# correction 1 for macropodids
PG.F.pred2 <- PG.F.pred1 * MACROPOD.F.corr1
PG.F.pred2

# correction 2 for macropodids (allometrically predicted)
PG.F.pred <- as.numeric(PG.F.pred1 * (365/(MACROPOD.F.corr.a + MACROPOD.F.corr.b*log10(PG.mass))))
PG.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
PG.alpha1 <- ceiling(exp(-1.34 + (0.214*log(PG.mass*1000))))
PG.alpha1
PG.alpha <- round(PG.alpha1 * MACROPOD.alpha.corr, 0)
PG.alpha

## define m function with age
PG.m.vec <- c(rep(0, PG.alpha-1), rep(0.5*PG.F.pred, round(PG.alpha/2,0)), rep(PG.F.pred, (PG.age.max+1-((PG.alpha-1+round(PG.alpha/2,0))))))
PG.m.sd.vec <- 0.05*PG.m.vec
plot(PG.age.vec, PG.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
PG.m.dat <- data.frame(PG.age.vec, PG.m.vec)
param.init <- c(0.6, 4, -5)
PG.fit.logp <- nls(PG.m.vec ~ a / (1+(PG.age.vec/b)^c), 
                   data = PG.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PG.fit.logp.summ <- summary(PG.fit.logp)
plot(PG.age.vec, PG.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
PG.age.vec.cont <- seq(0,max(PG.age.vec),1)
PG.pred.p.m <- coef(PG.fit.logp)[1] / (1+(PG.age.vec.cont/coef(PG.fit.logp)[2])^coef(PG.fit.logp)[3])
PG.pred.p.mm <- ifelse(PG.pred.p.m > 1, 1, PG.pred.p.m)
lines(PG.age.vec.cont, PG.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
PG.s.tran <- ln.a.s + b.s*log(PG.mass*1000) + log(1)
PG.s.ad.yr <- exp(-exp(PG.s.tran))
PG.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.93*PG.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.3 # rate of mortality decline (also known as bt)
a2 <- 1 - PG.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.4e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.17 # rate of mortality increase
longev <- PG.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
PG.Sx <- c(0.89*PG.s.ad.yr, 1 - qx)
plot(x, PG.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
PG.s.sd.vec <- 0.05*PG.Sx

## create matrix
PG.popmat <- matrix(data = 0, nrow=PG.age.max+1, ncol=PG.age.max+1)
diag(PG.popmat[2:(PG.age.max+1),]) <- PG.Sx[-(PG.age.max+1)]
PG.popmat[PG.age.max+1,PG.age.max+1] <- PG.Sx[PG.age.max+1]
PG.popmat[1,] <- PG.pred.p.mm
colnames(PG.popmat) <- c(0:PG.age.max)
rownames(PG.popmat) <- c(0:PG.age.max)
PG.popmat.orig <- PG.popmat ## save original matrix

## matrix properties
max.lambda(PG.popmat.orig) ## 1-yr lambda
PG.lm.pred
max.r(PG.popmat.orig) # rate of population change, 1-yr
PG.ssd <- stable.stage.dist(PG.popmat.orig) ## stable stage distribution
plot(PG.age.vec, PG.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(PG.popmat.orig, PG.age.max) # reproductive value
PG.gen.l <- G.val(PG.popmat.orig, PG.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
PG.pop.found <- round(area*PG.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
PG.init.vec <- PG.ssd * PG.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*PG.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

PG.tot.F <- sum(PG.popmat.orig[1,])
PG.popmat <- PG.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
PG.n.mat <- matrix(0, nrow=PG.age.max+1,ncol=(t+1))
PG.n.mat[,1] <- PG.init.vec

## set up projection loop
for (i in 1:t) {
  PG.n.mat[,i+1] <- PG.popmat %*% PG.n.mat[,i]
}

PG.n.pred <- colSums(PG.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(PG.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
PG.K.max <- 1*PG.pop.found
PG.K.vec <- c(1, PG.K.max/2, 0.75*PG.K.max, PG.K.max) 
PG.red.vec <- c(1,0.95,0.885,0.805)
plot(PG.K.vec, PG.red.vec,pch=19,type="b")
PG.Kred.dat <- data.frame(PG.K.vec, PG.red.vec)

# logistic power function a/(1+(x/b)^c)
PG.param.init <- c(1, 2*PG.K.max, 2)
PG.fit.lp <- nls(PG.red.vec ~ a/(1+(PG.K.vec/b)^c), 
                 data = PG.Kred.dat,
                 algorithm = "port",
                 start = c(a = PG.param.init[1], b = PG.param.init[2], c = PG.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PG.fit.lp.summ <- summary(PG.fit.lp)
plot(PG.K.vec, PG.red.vec, pch=19,xlab="N",ylab="reduction factor")
PG.K.vec.cont <- seq(1,2*PG.pop.found,1)
PG.pred.lp.fx <- coef(PG.fit.lp)[1]/(1+(PG.K.vec.cont/coef(PG.fit.lp)[2])^coef(PG.fit.lp)[3])
lines(PG.K.vec.cont, PG.pred.lp.fx, lty=3,lwd=3,col="red")

PG.a.lp <- coef(PG.fit.lp)[1]
PG.b.lp <- coef(PG.fit.lp)[2]
PG.c.lp <- coef(PG.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
PG.n.mat <- matrix(0, nrow=PG.age.max+1, ncol=(t+1))
PG.n.mat[,1] <- PG.init.vec
PG.popmat <- PG.popmat.orig

## set up projection loop
for (i in 1:t) {
  PG.totN.i <- sum(PG.n.mat[,i])
  PG.pred.red <- as.numeric(PG.a.lp/(1+(PG.totN.i/PG.b.lp)^PG.c.lp))
  diag(PG.popmat[2:(PG.age.max+1),]) <- (PG.Sx[-(PG.age.max+1)])*PG.pred.red
  PG.popmat[PG.age.max+1,PG.age.max+1] <- (PG.Sx[PG.age.max+1])*PG.pred.red
  PG.popmat[1,] <- PG.pred.p.mm
  PG.n.mat[,i+1] <- PG.popmat %*% PG.n.mat[,i]
}

PG.n.pred <- colSums(PG.n.mat)
plot(yrs, PG.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=PG.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

PG.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PG.s.arr <- PG.m.arr <- array(data=NA, dim=c(t+1, PG.age.max+1, iter))
PG.pred.red.mn <- rep(0,iter)

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
    
    # stochastic fertilty sampler (gaussian)
    PG.fert.stch <- rnorm(length(PG.popmat[,1]), PG.pred.p.mm, PG.m.sd.vec)
    PG.m.arr[i,,e] <- ifelse(PG.fert.stch < 0, 0, PG.fert.stch)
    
    PG.totN.i <- sum(PG.n.mat[,i], na.rm=T)
    PG.pred.red.vec[i] <- PG.a.lp/(1+(PG.totN.i/PG.b.lp)^PG.c.lp)
    
    diag(PG.popmat[2:(PG.age.max+1),]) <- (PG.s.stoch[-(PG.age.max+1)])*PG.pred.red.vec[i]
    PG.popmat[PG.age.max+1,PG.age.max+1] <- (PG.s.stoch[PG.age.max+1])*PG.pred.red.vec[i]
    PG.popmat[1,] <- PG.m.arr[i,,e]
    PG.n.mat[,i+1] <- PG.popmat %*% PG.n.mat[,i]
    
    PG.s.arr[i,,e] <- PG.s.stoch * PG.pred.red.vec[i]
    
  } # end i loop
  
  PG.n.sums.mat[e,] <- ((as.vector(colSums(PG.n.mat))/PG.pop.found))
  PG.pred.red.mn[e] <- mean(PG.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

PG.n.md <- apply(PG.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
PG.n.up <- apply(PG.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PG.n.lo <- apply(PG.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## STHENURUS (stirlingi) (SS)

## mass
SS.mass <- 150 # Sthenurus stirlingi (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1442-9993.2004.01389.x)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
SS.rm.pred <- 10^(0.6914 - (0.2622*log10(SS.mass*1000)))
SS.rm.pred
SS.lm.pred <- exp(SS.rm.pred)
SS.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
SS.D.pred <- (10^(4.196 - (0.74*log10(SS.mass*1000))))/2 # divided by 2 for females only
SS.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
SS.age.max1 <- round(10^(0.89 + (0.13*log10(SS.mass*1000))), 0)
SS.age.max1
SS.age.max <- ceiling(13/30*SS.age.max1) # (amount overestimated for OR from Jonzén et al. 2010-J Anim Ecol)
SS.age.max

## age vector
SS.age.vec <- 0:SS.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
SS.F.pred1 <- exp(2.719 - (0.211*log(SS.mass*1000)))/2 # divided by 2 for females
SS.F.pred1

# correction 1 for macropodids
SS.F.pred2 <- SS.F.pred1 * MACROPOD.F.corr1
SS.F.pred2

# correction 2 for macropodids (allometrically predicted)
SS.F.pred <- as.numeric(SS.F.pred1 * (365/(MACROPOD.F.corr.a + MACROPOD.F.corr.b*log10(SS.mass))))
SS.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
SS.alpha1 <- ceiling(exp(-1.34 + (0.214*log(SS.mass*1000))))
SS.alpha1
SS.alpha <- round(SS.alpha1 * MACROPOD.alpha.corr, 0)
SS.alpha

## define m function with age
SS.m.vec <- c(rep(0, SS.alpha-1), rep(0.5*SS.F.pred, round(SS.alpha/2,0)), rep(SS.F.pred, (SS.age.max+1-((SS.alpha-1+round(SS.alpha/2,0))))))
SS.m.sd.vec <- 0.05*SS.m.vec
plot(SS.age.vec, SS.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
SS.m.dat <- data.frame(SS.age.vec, SS.m.vec)
param.init <- c(0.6, 4, -5)
SS.fit.logp <- nls(SS.m.vec ~ a / (1+(SS.age.vec/b)^c), 
                   data = SS.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
SS.fit.logp.summ <- summary(SS.fit.logp)
plot(SS.age.vec, SS.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
SS.age.vec.cont <- seq(0,max(SS.age.vec),1)
SS.pred.p.m <- coef(SS.fit.logp)[1] / (1+(SS.age.vec.cont/coef(SS.fit.logp)[2])^coef(SS.fit.logp)[3])
SS.pred.p.mm <- ifelse(SS.pred.p.m > 1, 1, SS.pred.p.m)
lines(SS.age.vec.cont, SS.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
SS.s.tran <- ln.a.s + b.s*log(SS.mass*1000) + log(1)
SS.s.ad.yr <- exp(-exp(SS.s.tran))
SS.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.95*SS.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.7 # rate of mortality decline (also known as bt)
a2 <- 1 - SS.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.3e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.07 # rate of mortality increase
longev <- SS.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
SS.Sx <- c(0.925*SS.s.ad.yr, 1 - qx)
plot(x, SS.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
SS.s.sd.vec <- 0.05*SS.Sx

## create matrix
SS.popmat <- matrix(data = 0, nrow=SS.age.max+1, ncol=SS.age.max+1)
diag(SS.popmat[2:(SS.age.max+1),]) <- SS.Sx[-(SS.age.max+1)]
SS.popmat[SS.age.max+1,SS.age.max+1] <- SS.Sx[SS.age.max+1]
SS.popmat[1,] <- SS.pred.p.mm
colnames(SS.popmat) <- c(0:SS.age.max)
rownames(SS.popmat) <- c(0:SS.age.max)
SS.popmat.orig <- SS.popmat ## save original matrix

## matrix properties
max.lambda(SS.popmat.orig) ## 1-yr lambda
SS.lm.pred
max.r(SS.popmat.orig) # rate of population change, 1-yr
SS.ssd <- stable.stage.dist(SS.popmat.orig) ## stable stage distribution
plot(SS.age.vec, SS.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(SS.popmat.orig, SS.age.max) # reproductive value
SS.gen.l <- G.val(SS.popmat.orig, SS.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
SS.pop.found <- round(area*SS.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
SS.init.vec <- SS.ssd * SS.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*SS.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

SS.tot.F <- sum(SS.popmat.orig[1,])
SS.popmat <- SS.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
SS.n.mat <- matrix(0, nrow=SS.age.max+1,ncol=(t+1))
SS.n.mat[,1] <- SS.init.vec

## set up projection loop
for (i in 1:t) {
  SS.n.mat[,i+1] <- SS.popmat %*% SS.n.mat[,i]
}

SS.n.pred <- colSums(SS.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(SS.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
SS.K.max <- 1*SS.pop.found
SS.K.vec <- c(1, SS.K.max/2, 0.75*SS.K.max, SS.K.max) 
SS.red.vec <- c(1,0.935,0.855,0.78)
plot(SS.K.vec, SS.red.vec,pch=19,type="b")
SS.Kred.dat <- data.frame(SS.K.vec, SS.red.vec)

# logistic power function a/(1+(x/b)^c)
SS.param.init <- c(1, 2*SS.K.max, 2)
SS.fit.lp <- nls(SS.red.vec ~ a/(1+(SS.K.vec/b)^c), 
                 data = SS.Kred.dat,
                 algorithm = "port",
                 start = c(a = SS.param.init[1], b = SS.param.init[2], c = SS.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
SS.fit.lp.summ <- summary(SS.fit.lp)
plot(SS.K.vec, SS.red.vec, pch=19,xlab="N",ylab="reduction factor")
SS.K.vec.cont <- seq(1,2*SS.pop.found,1)
SS.pred.lp.fx <- coef(SS.fit.lp)[1]/(1+(SS.K.vec.cont/coef(SS.fit.lp)[2])^coef(SS.fit.lp)[3])
lines(SS.K.vec.cont, SS.pred.lp.fx, lty=3,lwd=3,col="red")

SS.a.lp <- coef(SS.fit.lp)[1]
SS.b.lp <- coef(SS.fit.lp)[2]
SS.c.lp <- coef(SS.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
SS.n.mat <- matrix(0, nrow=SS.age.max+1, ncol=(t+1))
SS.n.mat[,1] <- SS.init.vec
SS.popmat <- SS.popmat.orig

## set up projection loop
for (i in 1:t) {
  SS.totN.i <- sum(SS.n.mat[,i])
  SS.pred.red <- as.numeric(SS.a.lp/(1+(SS.totN.i/SS.b.lp)^SS.c.lp))
  diag(SS.popmat[2:(SS.age.max+1),]) <- (SS.Sx[-(SS.age.max+1)])*SS.pred.red
  SS.popmat[SS.age.max+1,SS.age.max+1] <- (SS.Sx[SS.age.max+1])*SS.pred.red
  SS.popmat[1,] <- SS.pred.p.mm
  SS.n.mat[,i+1] <- SS.popmat %*% SS.n.mat[,i]
}

SS.n.pred <- colSums(SS.n.mat)
plot(yrs, SS.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=SS.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

SS.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
SS.s.arr <- SS.m.arr <- array(data=NA, dim=c(t+1, SS.age.max+1, iter))
SS.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  SS.popmat <- SS.popmat.orig
  
  SS.n.mat <- matrix(0, nrow=SS.age.max+1,ncol=(t+1))
  SS.n.mat[,1] <- SS.init.vec
  SS.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    SS.s.alpha <- estBetaParams(SS.Sx, SS.s.sd.vec^2)$alpha
    SS.s.beta <- estBetaParams(SS.Sx, SS.s.sd.vec^2)$beta
    SS.s.stoch <- rbeta(length(SS.s.alpha), SS.s.alpha, SS.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    SS.fert.stch <- rnorm(length(SS.popmat[,1]), SS.pred.p.mm, SS.m.sd.vec)
    SS.m.arr[i,,e] <- ifelse(SS.fert.stch < 0, 0, SS.fert.stch)
    
    SS.totN.i <- sum(SS.n.mat[,i], na.rm=T)
    SS.pred.red.vec[i] <- SS.a.lp/(1+(SS.totN.i/SS.b.lp)^SS.c.lp)
    
    diag(SS.popmat[2:(SS.age.max+1),]) <- (SS.s.stoch[-(SS.age.max+1)])*SS.pred.red.vec[i]
    SS.popmat[SS.age.max+1,SS.age.max+1] <- (SS.s.stoch[SS.age.max+1])*SS.pred.red.vec[i]
    SS.popmat[1,] <- SS.m.arr[i,,e]
    SS.n.mat[,i+1] <- SS.popmat %*% SS.n.mat[,i]
    
    SS.s.arr[i,,e] <- SS.s.stoch * SS.pred.red.vec[i]
    
  } # end i loop
  
  SS.n.sums.mat[e,] <- ((as.vector(colSums(SS.n.mat))/SS.pop.found))
  SS.pred.red.mn[e] <- mean(SS.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

SS.n.md <- apply(SS.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
SS.n.up <- apply(SS.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
SS.n.lo <- apply(SS.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## PROTEMNODON (anak) (PT)

# mass
PT.mass <- 130 # Protemnodon anak (Johnson et al. 2006)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
PT.rm.pred <- 10^(0.6914 - (0.2622*log10(PT.mass*1000)))
PT.rm.pred
PT.lm.pred <- exp(PT.rm.pred)
PT.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
PT.D.pred <- (10^(4.196 - (0.74*log10(PT.mass*1000))))/2 # divided by 2 for females only
PT.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
PT.age.max1 <- round(10^(0.89 + (0.13*log10(PT.mass*1000))), 0)
PT.age.max1
PT.age.max <- ceiling(13/30*PT.age.max1) # (amount overestimated for OR from Jonzén et al. 2010-J Anim Ecol)
PT.age.max

## age vector
PT.age.vec <- 0:PT.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
PT.F.pred1 <- exp(2.719 - (0.211*log(PT.mass*1000)))/2 # divided by 2 for females
PT.F.pred1

# correction 1 for macropodids
PT.F.pred2 <- PT.F.pred1 * MACROPOD.F.corr1
PT.F.pred2

# correction 2 for macropodids (allometrically predicted)
PT.F.pred <- as.numeric(PT.F.pred1 * (365/(MACROPOD.F.corr.a + MACROPOD.F.corr.b*log10(PT.mass))))
PT.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
PT.alpha1 <- ceiling(exp(-1.34 + (0.214*log(PT.mass*1000))))
PT.alpha1
PT.alpha <- round(PT.alpha1 * MACROPOD.alpha.corr, 0)
PT.alpha

## define m function with age
PT.m.vec <- c(rep(0, PT.alpha-1), rep(0.5*PT.F.pred, round(PT.alpha/2,0)), rep(PT.F.pred, (PT.age.max+1-((PT.alpha-1+round(PT.alpha/2,0))))))
PT.m.sd.vec <- 0.05*PT.m.vec
plot(PT.age.vec, PT.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
PT.m.dat <- data.frame(PT.age.vec, PT.m.vec)
param.init <- c(0.6, 4, -5)
PT.fit.logp <- nls(PT.m.vec ~ a / (1+(PT.age.vec/b)^c), 
                   data = PT.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PT.fit.logp.summ <- summary(PT.fit.logp)
plot(PT.age.vec, PT.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
PT.age.vec.cont <- seq(0,max(PT.age.vec),1)
PT.pred.p.m <- coef(PT.fit.logp)[1] / (1+(PT.age.vec.cont/coef(PT.fit.logp)[2])^coef(PT.fit.logp)[3])
PT.pred.p.mm <- ifelse(PT.pred.p.m > 1, 1, PT.pred.p.m)
lines(PT.age.vec.cont, PT.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
PT.s.tran <- ln.a.s + b.s*log(PT.mass*1000) + log(1)
PT.s.ad.yr <- exp(-exp(PT.s.tran))
PT.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.965*PT.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.9 # rate of mortality decline (also known as bt)
a2 <- 1 - PT.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.05 # rate of mortality increase
longev <- PT.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
PT.Sx <- c(0.935*PT.s.ad.yr, 1 - qx)
plot(x, PT.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
PT.s.sd.vec <- 0.05*PT.Sx

## create matrix
PT.popmat <- matrix(data = 0, nrow=PT.age.max+1, ncol=PT.age.max+1)
diag(PT.popmat[2:(PT.age.max+1),]) <- PT.Sx[-(PT.age.max+1)]
PT.popmat[PT.age.max+1,PT.age.max+1] <- PT.Sx[PT.age.max+1]
PT.popmat[1,] <- PT.pred.p.mm
colnames(PT.popmat) <- c(0:PT.age.max)
rownames(PT.popmat) <- c(0:PT.age.max)
PT.popmat.orig <- PT.popmat ## save original matrix

## matrix properties
max.lambda(PT.popmat.orig) ## 1-yr lambda
PT.lm.pred
max.r(PT.popmat.orig) # rate of population change, 1-yr
PT.ssd <- stable.stage.dist(PT.popmat.orig) ## stable stage distribution
plot(PT.age.vec, PT.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(PT.popmat.orig, PT.age.max) # reproductive value
PT.gen.l <- G.val(PT.popmat.orig, PT.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
PT.pop.found <- round(area*PT.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
PT.init.vec <- PT.ssd * PT.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*PT.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

PT.tot.F <- sum(PT.popmat.orig[1,])
PT.popmat <- PT.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
PT.n.mat <- matrix(0, nrow=PT.age.max+1,ncol=(t+1))
PT.n.mat[,1] <- PT.init.vec

## set up projection loop
for (i in 1:t) {
  PT.n.mat[,i+1] <- PT.popmat %*% PT.n.mat[,i]
}

PT.n.pred <- colSums(PT.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(PT.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
PT.K.max <- 1*PT.pop.found
PT.K.vec <- c(1, PT.K.max/2, 0.75*PT.K.max, PT.K.max) 
PT.red.vec <- c(1,0.935,0.855,0.77)
plot(PT.K.vec, PT.red.vec,pch=19,type="b")
PT.Kred.dat <- data.frame(PT.K.vec, PT.red.vec)

# logistic power function a/(1+(x/b)^c)
PT.param.init <- c(1, 2*PT.K.max, 2)
PT.fit.lp <- nls(PT.red.vec ~ a/(1+(PT.K.vec/b)^c), 
                 data = PT.Kred.dat,
                 algorithm = "port",
                 start = c(a = PT.param.init[1], b = PT.param.init[2], c = PT.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PT.fit.lp.summ <- summary(PT.fit.lp)
plot(PT.K.vec, PT.red.vec, pch=19,xlab="N",ylab="reduction factor")
PT.K.vec.cont <- seq(1,2*PT.pop.found,1)
PT.pred.lp.fx <- coef(PT.fit.lp)[1]/(1+(PT.K.vec.cont/coef(PT.fit.lp)[2])^coef(PT.fit.lp)[3])
lines(PT.K.vec.cont, PT.pred.lp.fx, lty=3,lwd=3,col="red")

PT.a.lp <- coef(PT.fit.lp)[1]
PT.b.lp <- coef(PT.fit.lp)[2]
PT.c.lp <- coef(PT.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
PT.n.mat <- matrix(0, nrow=PT.age.max+1, ncol=(t+1))
PT.n.mat[,1] <- PT.init.vec
PT.popmat <- PT.popmat.orig

## set up projection loop
for (i in 1:t) {
  PT.totN.i <- sum(PT.n.mat[,i])
  PT.pred.red <- as.numeric(PT.a.lp/(1+(PT.totN.i/PT.b.lp)^PT.c.lp))
  diag(PT.popmat[2:(PT.age.max+1),]) <- (PT.Sx[-(PT.age.max+1)])*PT.pred.red
  PT.popmat[PT.age.max+1,PT.age.max+1] <- (PT.Sx[PT.age.max+1])*PT.pred.red
  PT.popmat[1,] <- PT.pred.p.mm
  PT.n.mat[,i+1] <- PT.popmat %*% PT.n.mat[,i]
}

PT.n.pred <- colSums(PT.n.mat)
plot(yrs, PT.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=PT.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

PT.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PT.s.arr <- PT.m.arr <- array(data=NA, dim=c(t+1, PT.age.max+1, iter))
PT.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  PT.popmat <- PT.popmat.orig
  
  PT.n.mat <- matrix(0, nrow=PT.age.max+1,ncol=(t+1))
  PT.n.mat[,1] <- PT.init.vec
  PT.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    PT.s.alpha <- estBetaParams(PT.Sx, PT.s.sd.vec^2)$alpha
    PT.s.beta <- estBetaParams(PT.Sx, PT.s.sd.vec^2)$beta
    PT.s.stoch <- rbeta(length(PT.s.alpha), PT.s.alpha, PT.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    PT.fert.stch <- rnorm(length(PT.popmat[,1]), PT.pred.p.mm, PT.m.sd.vec)
    PT.m.arr[i,,e] <- ifelse(PT.fert.stch < 0, 0, PT.fert.stch)
    
    PT.totN.i <- sum(PT.n.mat[,i], na.rm=T)
    PT.pred.red.vec[i] <- PT.a.lp/(1+(PT.totN.i/PT.b.lp)^PT.c.lp)
    
    diag(PT.popmat[2:(PT.age.max+1),]) <- (PT.s.stoch[-(PT.age.max+1)])*PT.pred.red.vec[i]
    PT.popmat[PT.age.max+1,PT.age.max+1] <- (PT.s.stoch[PT.age.max+1])*PT.pred.red.vec[i]
    PT.popmat[1,] <- PT.m.arr[i,,e]
    PT.n.mat[,i+1] <- PT.popmat %*% PT.n.mat[,i]
    
    PT.s.arr[i,,e] <- PT.s.stoch * PT.pred.red.vec[i]
    
  } # end i loop
  
  PT.n.sums.mat[e,] <- ((as.vector(colSums(PT.n.mat))/PT.pop.found))
  PT.pred.red.mn[e] <- mean(PT.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

PT.n.md <- apply(PT.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
PT.n.up <- apply(PT.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PT.n.lo <- apply(PT.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## SIMOSTHENURUS (occidentalis) (SO)

# mass
SO.mass <- 120 # Simosthenurus occidentalis (Johnson 2006, in Webb 2008)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
SO.rm.pred <- 10^(0.6914 - (0.2622*log10(SO.mass*1000)))
SO.rm.pred
SO.lm.pred <- exp(SO.rm.pred)
SO.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
SO.D.pred <- (10^(4.196 - (0.74*log10(SO.mass*1000))))/2 # divided by 2 for females only
SO.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
SO.age.max1 <- round(10^(0.89 + (0.13*log10(SO.mass*1000))), 0)
SO.age.max1
SO.age.max <- ceiling(13/30*SO.age.max1) # (amount overestimated for OR from Jonzén et al. 2010-J Anim Ecol)
SO.age.max

## age vector
SO.age.vec <- 0:SO.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
SO.F.pred1<- exp(2.719 - (0.211*log(SO.mass*1000)))/2 # divided by 2 for females
SO.F.pred1

# correction 1 for macropodids
SO.F.pred2 <- SO.F.pred1 * MACROPOD.F.corr1
SO.F.pred2

# correction 2 for macropodids (allometrically predicted)
SO.F.pred <- as.numeric(SO.F.pred1 * (365/(MACROPOD.F.corr.a + MACROPOD.F.corr.b*log10(SO.mass))))
SO.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
SO.alpha1 <- ceiling(exp(-1.34 + (0.214*log(SO.mass*1000))))
SO.alpha1
SO.alpha <- round(SO.alpha1 * MACROPOD.alpha.corr, 0)
SO.alpha

## define m function with age
SO.m.vec <- c(rep(0, SO.alpha-1), rep(0.5*SO.F.pred, round(SO.alpha/2,0)), rep(SO.F.pred, (SO.age.max+1-((SO.alpha-1+round(SO.alpha/2,0))))))
SO.m.sd.vec <- 0.05*SO.m.vec
plot(SO.age.vec, SO.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
SO.m.dat <- data.frame(SO.age.vec, SO.m.vec)
param.init <- c(0.6, 4, -5)
SO.fit.logp <- nls(SO.m.vec ~ a / (1+(SO.age.vec/b)^c), 
                   data = SO.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
SO.fit.logp.summ <- summary(SO.fit.logp)
plot(SO.age.vec, SO.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
SO.age.vec.cont <- seq(0,max(SO.age.vec),1)
SO.pred.p.m <- coef(SO.fit.logp)[1] / (1+(SO.age.vec.cont/coef(SO.fit.logp)[2])^coef(SO.fit.logp)[3])
SO.pred.p.mm <- ifelse(SO.pred.p.m > 1, 1, SO.pred.p.m)
lines(SO.age.vec.cont, SO.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
SO.s.tran <- ln.a.s + b.s*log(SO.mass*1000) + log(1)
SO.s.ad.yr <- exp(-exp(SO.s.tran))
SO.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.95*SO.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.5 # rate of mortality decline (also known as bt)
a2 <- 1 - SO.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.1 # rate of mortality increase
longev <- SO.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
SO.Sx <- c(0.925*SO.s.ad.yr, 1 - qx)
plot(x, SO.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
SO.s.sd.vec <- 0.05*SO.Sx

## create matrix
SO.popmat <- matrix(data = 0, nrow=SO.age.max+1, ncol=SO.age.max+1)
diag(SO.popmat[2:(SO.age.max+1),]) <- SO.Sx[-(SO.age.max+1)]
SO.popmat[SO.age.max+1,SO.age.max+1] <- SO.Sx[SO.age.max+1]
SO.popmat[1,] <- SO.pred.p.mm
colnames(SO.popmat) <- c(0:SO.age.max)
rownames(SO.popmat) <- c(0:SO.age.max)
SO.popmat.orig <- SO.popmat ## save original matrix

## matrix properties
max.lambda(SO.popmat.orig) ## 1-yr lambda
SO.lm.pred
max.r(SO.popmat.orig) # rate of population change, 1-yr
SO.ssd <- stable.stage.dist(SO.popmat.orig) ## stable stage distribution
plot(SO.age.vec, SO.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(SO.popmat.orig, SO.age.max) # reproductive value
SO.gen.l <- G.val(SO.popmat.orig, SO.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
SO.pop.found <- round(area*SO.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
SO.init.vec <- SO.ssd * SO.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*SO.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

SO.tot.F <- sum(SO.popmat.orig[1,])
SO.popmat <- SO.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
SO.n.mat <- matrix(0, nrow=SO.age.max+1,ncol=(t+1))
SO.n.mat[,1] <- SO.init.vec

## set up projection loop
for (i in 1:t) {
  SO.n.mat[,i+1] <- SO.popmat %*% SO.n.mat[,i]
}

SO.n.pred <- colSums(SO.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(SO.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
SO.K.max <- 1*SO.pop.found
SO.K.vec <- c(1, SO.K.max/2, 0.75*SO.K.max, SO.K.max) 
SO.red.vec <- c(1,0.935,0.854,0.769)
plot(SO.K.vec, SO.red.vec,pch=19,type="b")
SO.Kred.dat <- data.frame(SO.K.vec, SO.red.vec)

# logistic power function a/(1+(x/b)^c)
SO.param.init <- c(1, 2*SO.K.max, 2)
SO.fit.lp <- nls(SO.red.vec ~ a/(1+(SO.K.vec/b)^c), 
                 data = SO.Kred.dat,
                 algorithm = "port",
                 start = c(a = SO.param.init[1], b = SO.param.init[2], c = SO.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
SO.fit.lp.summ <- summary(SO.fit.lp)
plot(SO.K.vec, SO.red.vec, pch=19,xlab="N",ylab="reduction factor")
SO.K.vec.cont <- seq(1,2*SO.pop.found,1)
SO.pred.lp.fx <- coef(SO.fit.lp)[1]/(1+(SO.K.vec.cont/coef(SO.fit.lp)[2])^coef(SO.fit.lp)[3])
lines(SO.K.vec.cont, SO.pred.lp.fx, lty=3,lwd=3,col="red")

SO.a.lp <- coef(SO.fit.lp)[1]
SO.b.lp <- coef(SO.fit.lp)[2]
SO.c.lp <- coef(SO.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
SO.n.mat <- matrix(0, nrow=SO.age.max+1, ncol=(t+1))
SO.n.mat[,1] <- SO.init.vec
SO.popmat <- SO.popmat.orig

## set up projection loop
for (i in 1:t) {
  SO.totN.i <- sum(SO.n.mat[,i])
  SO.pred.red <- as.numeric(SO.a.lp/(1+(SO.totN.i/SO.b.lp)^SO.c.lp))
  diag(SO.popmat[2:(SO.age.max+1),]) <- (SO.Sx[-(SO.age.max+1)])*SO.pred.red
  SO.popmat[SO.age.max+1,SO.age.max+1] <- (SO.Sx[SO.age.max+1])*SO.pred.red
  SO.popmat[1,] <- SO.pred.p.mm
  SO.n.mat[,i+1] <- SO.popmat %*% SO.n.mat[,i]
}

SO.n.pred <- colSums(SO.n.mat)
plot(yrs, SO.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=SO.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

SO.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
SO.s.arr <- SO.m.arr <- array(data=NA, dim=c(t+1, SO.age.max+1, iter))
SO.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  SO.popmat <- SO.popmat.orig
  
  SO.n.mat <- matrix(0, nrow=SO.age.max+1,ncol=(t+1))
  SO.n.mat[,1] <- SO.init.vec
  SO.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    SO.s.alpha <- estBetaParams(SO.Sx, SO.s.sd.vec^2)$alpha
    SO.s.beta <- estBetaParams(SO.Sx, SO.s.sd.vec^2)$beta
    SO.s.stoch <- rbeta(length(SO.s.alpha), SO.s.alpha, SO.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    SO.fert.stch <- rnorm(length(SO.popmat[,1]), SO.pred.p.mm, SO.m.sd.vec)
    SO.m.arr[i,,e] <- ifelse(SO.fert.stch < 0, 0, SO.fert.stch)
    
    SO.totN.i <- sum(SO.n.mat[,i], na.rm=T)
    SO.pred.red.vec[i] <- SO.a.lp/(1+(SO.totN.i/SO.b.lp)^SO.c.lp)
    
    diag(SO.popmat[2:(SO.age.max+1),]) <- (SO.s.stoch[-(SO.age.max+1)])*SO.pred.red.vec[i]
    SO.popmat[SO.age.max+1,SO.age.max+1] <- (SO.s.stoch[SO.age.max+1])*SO.pred.red.vec[i]
    SO.popmat[1,] <- SO.m.arr[i,,e]
    SO.n.mat[,i+1] <- SO.popmat %*% SO.n.mat[,i]

    SO.s.arr[i,,e] <- SO.s.stoch * SO.pred.red.vec[i]
    
  } # end i loop
  
  SO.n.sums.mat[e,] <- ((as.vector(colSums(SO.n.mat))/SO.pop.found))
  SO.pred.red.mn[e] <- mean(SO.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

SO.n.md <- apply(SO.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
SO.n.up <- apply(SO.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
SO.n.lo <- apply(SO.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations



##############################
## METASTHENURUS (newtonae) (MN)

# mass
MN.mass <- 55 # Metasthenurus newtonae (Johnson et al. 2006)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
MN.rm.pred <- 10^(0.6914 - (0.2622*log10(MN.mass*1000)))
MN.rm.pred
MN.lm.pred <- exp(MN.rm.pred)
MN.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
MN.D.pred <- (10^(4.196 - (0.74*log10(MN.mass*1000))))/2 # divided by 2 for females only
MN.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
MN.age.max1 <- round(10^(0.89 + (0.13*log10(MN.mass*1000))), 0)
MN.age.max1
MN.age.max <- ceiling(13/30*MN.age.max1) # (amount overestimated for OR from Jonzén et al. 2010-J Anim Ecol)
MN.age.max

## age vector
MN.age.vec <- 0:MN.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
MN.F.pred1 <- exp(2.719 - (0.211*log(MN.mass*1000)))/2 # divided by 2 for females
MN.F.pred1

# correction 1 for macropodids
MN.F.pred2 <- MN.F.pred1 * MACROPOD.F.corr1
MN.F.pred2

# correction 2 for macropodids (allometrically predicted)
MN.F.pred <- as.numeric(MN.F.pred1 * (365/(MACROPOD.F.corr.a + MACROPOD.F.corr.b*log10(MN.mass))))
MN.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
MN.alpha1 <- ceiling(exp(-1.34 + (0.214*log(MN.mass*1000))))
MN.alpha1
MN.alpha <- round(MN.alpha1 * MACROPOD.alpha.corr, 0)
MN.alpha

## define m function with age
MN.m.vec <- c(rep(0, MN.alpha-1), rep(0.5*MN.F.pred, round(MN.alpha/2,0)), rep(MN.F.pred, (MN.age.max+1-((MN.alpha-1+round(MN.alpha/2,0))))))
MN.m.sd.vec <- 0.05*MN.m.vec
plot(MN.age.vec, MN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
MN.m.dat <- data.frame(MN.age.vec, MN.m.vec)
param.init <- c(1, 2, -30)
MN.fit.logp <- nls(MN.m.vec ~ a / (1+(MN.age.vec/b)^c), 
                   data = MN.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
MN.fit.logp.summ <- summary(MN.fit.logp)
plot(MN.age.vec, MN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
MN.age.vec.cont <- seq(0,max(MN.age.vec),1)
MN.pred.p.m <- coef(MN.fit.logp)[1] / (1+(MN.age.vec.cont/coef(MN.fit.logp)[2])^coef(MN.fit.logp)[3])
MN.pred.p.mm <- ifelse(MN.pred.p.m > 1, 1, MN.pred.p.m)
lines(MN.age.vec.cont, MN.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
MN.s.tran <- ln.a.s + b.s*log(MN.mass*1000) + log(1)
MN.s.ad.yr <- exp(-exp(MN.s.tran))
MN.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.79*MN.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.5 # rate of mortality decline (also known as bt)
a2 <- 1 - MN.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.8e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.2 # rate of mortality increase
longev <- MN.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
MN.Sx <- c(0.76*MN.s.ad.yr, 1 - qx)
plot(x, MN.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
MN.s.sd.vec <- 0.05*MN.Sx

## create matrix
MN.popmat <- matrix(data = 0, nrow=MN.age.max+1, ncol=MN.age.max+1)
diag(MN.popmat[2:(MN.age.max+1),]) <- MN.Sx[-(MN.age.max+1)]
MN.popmat[MN.age.max+1,MN.age.max+1] <- MN.Sx[MN.age.max+1]
MN.popmat[1,] <- MN.pred.p.mm
colnames(MN.popmat) <- c(0:MN.age.max)
rownames(MN.popmat) <- c(0:MN.age.max)
MN.popmat.orig <- MN.popmat ## save original matrix

## matrix properties
max.lambda(MN.popmat.orig) ## 1-yr lambda
MN.lm.pred
max.r(MN.popmat.orig) # rate of population change, 1-yr
MN.ssd <- stable.stage.dist(MN.popmat.orig) ## stable stage distribution
plot(MN.age.vec, MN.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(MN.popmat.orig, MN.age.max) # reproductive value
MN.gen.l <- G.val(MN.popmat.orig, MN.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
MN.pop.found <- round(area*MN.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
MN.init.vec <- MN.ssd * MN.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*MN.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

MN.tot.F <- sum(MN.popmat.orig[1,])
MN.popmat <- MN.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
MN.n.mat <- matrix(0, nrow=MN.age.max+1,ncol=(t+1))
MN.n.mat[,1] <- MN.init.vec

## set up projection loop
for (i in 1:t) {
  MN.n.mat[,i+1] <- MN.popmat %*% MN.n.mat[,i]
}

MN.n.pred <- colSums(MN.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(MN.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
MN.K.max <- 1*MN.pop.found
MN.K.vec <- c(1, MN.K.max/2, 0.75*MN.K.max, MN.K.max) 
MN.red.vec <- c(1,0.915,0.81,0.700)
plot(MN.K.vec, MN.red.vec,pch=19,type="b")
MN.Kred.dat <- data.frame(MN.K.vec, MN.red.vec)

# logistic power function a/(1+(x/b)^c)
MN.param.init <- c(1, 2*MN.K.max, 2)
MN.fit.lp <- nls(MN.red.vec ~ a/(1+(MN.K.vec/b)^c), 
                 data = MN.Kred.dat,
                 algorithm = "port",
                 start = c(a = MN.param.init[1], b = MN.param.init[2], c = MN.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
MN.fit.lp.summ <- summary(MN.fit.lp)
plot(MN.K.vec, MN.red.vec, pch=19,xlab="N",ylab="reduction factor")
MN.K.vec.cont <- seq(1,2*MN.pop.found,1)
MN.pred.lp.fx <- coef(MN.fit.lp)[1]/(1+(MN.K.vec.cont/coef(MN.fit.lp)[2])^coef(MN.fit.lp)[3])
lines(MN.K.vec.cont, MN.pred.lp.fx, lty=3,lwd=3,col="red")

MN.a.lp <- coef(MN.fit.lp)[1]
MN.b.lp <- coef(MN.fit.lp)[2]
MN.c.lp <- coef(MN.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
MN.n.mat <- matrix(0, nrow=MN.age.max+1, ncol=(t+1))
MN.n.mat[,1] <- MN.init.vec
MN.popmat <- MN.popmat.orig

## set up projection loop
for (i in 1:t) {
  MN.totN.i <- sum(MN.n.mat[,i])
  MN.pred.red <- as.numeric(MN.a.lp/(1+(MN.totN.i/MN.b.lp)^MN.c.lp))
  diag(MN.popmat[2:(MN.age.max+1),]) <- (MN.Sx[-(MN.age.max+1)])*MN.pred.red
  MN.popmat[MN.age.max+1,MN.age.max+1] <- (MN.Sx[MN.age.max+1])*MN.pred.red
  MN.popmat[1,] <- MN.pred.p.mm
  MN.n.mat[,i+1] <- MN.popmat %*% MN.n.mat[,i]
}

MN.n.pred <- colSums(MN.n.mat)
plot(yrs, MN.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=MN.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

MN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
MN.s.arr <- MN.m.arr <- array(data=NA, dim=c(t+1, MN.age.max+1, iter))
MN.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  MN.popmat <- MN.popmat.orig
  
  MN.n.mat <- matrix(0, nrow=MN.age.max+1,ncol=(t+1))
  MN.n.mat[,1] <- MN.init.vec
  MN.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    MN.s.alpha <- estBetaParams(MN.Sx, MN.s.sd.vec^2)$alpha
    MN.s.beta <- estBetaParams(MN.Sx, MN.s.sd.vec^2)$beta
    MN.s.stoch <- rbeta(length(MN.s.alpha), MN.s.alpha, MN.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    MN.fert.stch <- rnorm(length(MN.popmat[,1]), MN.pred.p.mm, MN.m.sd.vec)
    MN.m.arr[i,,e] <- ifelse(MN.fert.stch < 0, 0, MN.fert.stch)
    
    MN.totN.i <- sum(MN.n.mat[,i], na.rm=T)
    MN.pred.red.vec[i] <- MN.a.lp/(1+(MN.totN.i/MN.b.lp)^MN.c.lp)
    
    diag(MN.popmat[2:(MN.age.max+1),]) <- (MN.s.stoch[-(MN.age.max+1)])*MN.pred.red.vec[i]
    MN.popmat[MN.age.max+1,MN.age.max+1] <- (MN.s.stoch[MN.age.max+1])*MN.pred.red.vec[i]
    MN.popmat[1,] <- MN.m.arr[i,,e]
    MN.n.mat[,i+1] <- MN.popmat %*% MN.n.mat[,i]

    MN.s.arr[i,,e] <- MN.s.stoch * MN.pred.red.vec[i]
    
  } # end i loop
  
  MN.n.sums.mat[e,] <- ((as.vector(colSums(MN.n.mat))/MN.pop.found))
  MN.pred.red.mn[e] <- mean(MN.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

MN.n.md <- apply(MN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
MN.n.up <- apply(MN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
MN.n.lo <- apply(MN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations



##############################
## OSPHRANTER (rufus) (OR)

# mass
OR.mass <- 25 # (Croft & Clacy 2008 Macropus rufus. The Mammals of Australia (eds S. van Dyck & R. Strahan), pp. 352-354. Reed New Holland, Sydney)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
OR.rm.pred <- 10^(0.6914 - (0.2622*log10(OR.mass*1000)))
OR.rm.pred
OR.lm.pred <- exp(OR.rm.pred)
OR.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
OR.D.pred <- (10^(4.196 - (0.74*log10(OR.mass*1000))))/2 # divided by 2 for females only
OR.D.pred # animals/km2
## 17, 0.6, 16/km2 (Fisher database); mean(17,0.6,16) = 17 (8.5/km2 for females only)

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
OR.age.max <- round(10^(0.89 + (0.13*log10(OR.mass*1000))), 0)
OR.age.max
OR.age.max <- 13 # (Jonzén et al. 2010-J Anim Ecol)

## age vector
OR.age.vec <- 0:OR.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
OR.F.pred <- exp(2.719 - (0.211*log(OR.mass*1000)))/2 # divided by 2 for females
OR.F.pred
OR.F.pred <- 1.5/2

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
OR.alpha <- ceiling(exp(-1.34 + (0.214*log(OR.mass*1000))))
OR.alpha
OR.alpha <- 2 # (Jonzén et al. 2010-J Anim Ecol)

## define m function with age
OR.m.vec <- c(rep(0, OR.alpha-1), rep(0.75*OR.F.pred, round(OR.alpha/2,0)), rep(OR.F.pred, (OR.age.max+1-((OR.alpha-1+round(OR.alpha/2,0))))))
OR.m.sd.vec <- 0.05*OR.m.vec
plot(OR.age.vec, OR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# Weibull y = a - b*exp(-c*x^d)
OR.m.dat <- data.frame(OR.age.vec, OR.m.vec)
param.init <- c(0.75, 0.75, 0.69, 5.1)
OR.fit.logp <- nls(OR.m.vec ~ a - b*exp(-c*OR.age.vec^d), 
                   data = OR.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3], d = param.init[4]),
                   trace = TRUE,      
                   nls.control(maxiter = 100000, tol = 1e-05, minFactor = 1/1024))
#OR.fit.logp.summ <- summary(OR.fit.logp)
plot(OR.age.vec, OR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
OR.age.vec.cont <- seq(0,max(OR.age.vec),1)
OR.pred.p.m <- coef(OR.fit.logp)[1] - coef(OR.fit.logp)[2]*exp(-coef(OR.fit.logp)[3]*OR.age.vec.cont^coef(OR.fit.logp)[4])
OR.pred.p.mm <- ifelse(OR.pred.p.m > 1, 1, OR.pred.p.m)
lines(OR.age.vec.cont, OR.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
OR.s.tran <- ln.a.s + b.s*log(OR.mass*1000) + log(1)
OR.s.ad.yr <- exp(-exp(OR.s.tran))
OR.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.99*OR.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.0 # rate of mortality decline (also known as bt)
a2 <- 1 - 1.01*OR.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.2e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.1 # rate of mortality increase
longev <- OR.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
OR.Sx <- c(0.975*OR.s.ad.yr, 1 - qx)
plot(x, OR.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
OR.s.sd.vec <- 0.05*OR.Sx

## create matrix
OR.popmat <- matrix(data = 0, nrow=OR.age.max+1, ncol=OR.age.max+1)
diag(OR.popmat[2:(OR.age.max+1),]) <- OR.Sx[-(OR.age.max+1)]
OR.popmat[OR.age.max+1,OR.age.max+1] <- OR.Sx[OR.age.max+1]
OR.popmat[1,] <- OR.pred.p.mm
colnames(OR.popmat) <- c(0:OR.age.max)
rownames(OR.popmat) <- c(0:OR.age.max)
OR.popmat.orig <- OR.popmat ## save original matrix

## matrix properties
max.lambda(OR.popmat.orig) ## 1-yr lambda
OR.lm.pred
max.r(OR.popmat.orig) # rate of population change, 1-yr
OR.ssd <- stable.stage.dist(OR.popmat.orig) ## stable stage distribution
plot(OR.age.vec, OR.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(OR.popmat.orig, OR.age.max) # reproductive value
OR.gen.l <- G.val(OR.popmat.orig, OR.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
OR.pop.found <- round(area*OR.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
OR.init.vec <- OR.ssd * OR.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*OR.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

OR.tot.F <- sum(OR.popmat.orig[1,])
OR.popmat <- OR.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
OR.n.mat <- matrix(0, nrow=OR.age.max+1,ncol=(t+1))
OR.n.mat[,1] <- OR.init.vec

## set up projection loop
for (i in 1:t) {
  OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]
}

OR.n.pred <- colSums(OR.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(OR.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
OR.K.max <- 1*OR.pop.found
OR.K.vec <- c(1, OR.K.max/2, 0.75*OR.K.max, OR.K.max) 
OR.red.vec <- c(1,0.88,0.77,0.63)
plot(OR.K.vec, OR.red.vec,pch=19,type="b")
OR.Kred.dat <- data.frame(OR.K.vec, OR.red.vec)

# logistic power function a/(1+(x/b)^c)
OR.param.init <- c(1, 2*OR.K.max, 2)
OR.fit.lp <- nls(OR.red.vec ~ a/(1+(OR.K.vec/b)^c), 
                 data = OR.Kred.dat,
                 algorithm = "port",
                 start = c(a = OR.param.init[1], b = OR.param.init[2], c = OR.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
OR.fit.lp.summ <- summary(OR.fit.lp)
plot(OR.K.vec, OR.red.vec, pch=19,xlab="N",ylab="reduction factor")
OR.K.vec.cont <- seq(1,2*OR.pop.found,1)
OR.pred.lp.fx <- coef(OR.fit.lp)[1]/(1+(OR.K.vec.cont/coef(OR.fit.lp)[2])^coef(OR.fit.lp)[3])
lines(OR.K.vec.cont, OR.pred.lp.fx, lty=3,lwd=3,col="red")

OR.a.lp <- coef(OR.fit.lp)[1]
OR.b.lp <- coef(OR.fit.lp)[2]
OR.c.lp <- coef(OR.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
OR.n.mat <- matrix(0, nrow=OR.age.max+1, ncol=(t+1))
OR.n.mat[,1] <- OR.init.vec
OR.popmat <- OR.popmat.orig

## set up projection loop
for (i in 1:t) {
  OR.totN.i <- sum(OR.n.mat[,i])
  OR.pred.red <- as.numeric(OR.a.lp/(1+(OR.totN.i/OR.b.lp)^OR.c.lp))
  diag(OR.popmat[2:(OR.age.max+1),]) <- (OR.Sx[-(OR.age.max+1)])*OR.pred.red
  OR.popmat[OR.age.max+1,OR.age.max+1] <- (OR.Sx[OR.age.max+1])*OR.pred.red
  OR.popmat[1,] <- OR.pred.p.mm
  OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]
}

OR.n.pred <- colSums(OR.n.mat)
plot(yrs, OR.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=OR.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

OR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
OR.s.arr <- OR.m.arr <- array(data=NA, dim=c(t+1, OR.age.max+1, iter))
OR.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  OR.popmat <- OR.popmat.orig
  
  OR.n.mat <- matrix(0, nrow=OR.age.max+1,ncol=(t+1))
  OR.n.mat[,1] <- OR.init.vec
  OR.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    OR.s.alpha <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$alpha
    OR.s.beta <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$beta
    OR.s.stoch <- rbeta(length(OR.s.alpha), OR.s.alpha, OR.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    OR.fert.stch <- rnorm(length(OR.popmat[,1]), OR.pred.p.mm, OR.m.sd.vec)
    OR.m.arr[i,,e] <- ifelse(OR.fert.stch < 0, 0, OR.fert.stch)
    
    OR.totN.i <- sum(OR.n.mat[,i], na.rm=T)
    OR.pred.red.vec[i] <- OR.a.lp/(1+(OR.totN.i/OR.b.lp)^OR.c.lp)
    
    diag(OR.popmat[2:(OR.age.max+1),]) <- (OR.s.stoch[-(OR.age.max+1)])*OR.pred.red.vec[i]
    OR.popmat[OR.age.max+1,OR.age.max+1] <- (OR.s.stoch[OR.age.max+1])*OR.pred.red.vec[i]
    OR.popmat[1,] <- OR.m.arr[i,,e]
    OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]

    OR.s.arr[i,,e] <- OR.s.stoch * OR.pred.red.vec[i]
    
  } # end i loop
  
  OR.n.sums.mat[e,] <- ((as.vector(colSums(OR.n.mat))/OR.pop.found))
  OR.pred.red.mn[e] <- mean(OR.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

OR.n.md <- apply(OR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
OR.n.up <- apply(OR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
OR.n.lo <- apply(OR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations



##############################
## NOTAMACROPUS (rufogriseus) (NR) 
## red-necked (Bennett’s) wallaby (Macropus rufogriseus rugogriseus)
## sources:
## Frith, H.J. and J.H. Calaby. 1969. Kangaroos. F.W. Cheshire, Melbourne, Australia.
## Morcombe, M. 1972. Michael Morcombe's Australian marsupials and other native mammals. Charles Scribner's Sons, New York.
## Strahan, R. 1983. The Australian Museum complete book of Australian mammals. Angus and Robertson, Sydney, Australia.
## Thomas, O. 1888. Catalogue of the marsupalia and monotremata in the collection of the British Museum (Natural History). Taylor and Francis, London, England.
## Troughton, E. 1962. Furred animals of Australia, 7th edition. Halstead Press, Sydney, Australia.
## Walker, E.P. 1964. Mammals of the world. John Hopkins University Press, Baltimore, MD.

# mass
NR.mass <- 14 # females; Strahan R. (1983) The Australian Museum Complete Book of Australian Mammals, 1st edn. Angus & Robertson Publishers, Sydney

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
NR.rm.pred <- 10^(0.6914 - (0.2622*log10(NR.mass*1000)))
NR.rm.pred
NR.lm.pred <- exp(NR.rm.pred)
NR.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
NR.D.pred <- (10^(4.196 - (0.74*log10(NR.mass*1000))))/2 # divided by 2 for females only
NR.D.pred # animals/km2
## 48/km2, 250/km2, 15 km2, 45/km2 from Fisher et al. 2001 database (for females, averages mean(24,125,7.5,22.5) = 24/km2)

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
NR.age.max <- round(mean(c(15, 18.5, 15.2)), 0) # Grzimek, B., ed. 1990. Grzimek's Animal Life Encyclopedia. Mammals I - IV. ed. Series. Grzimek, B. Vol. I-IV. New York: McGraw-Hill Publishing Company.
NR.age.max

## age vector
NR.age.vec <- 0:NR.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
NR.F.pred <- exp(2.719 - (0.211*log(NR.mass*1000)))/2 # divided by 2 for females
NR.F.pred
## inter-birth interval = 286 days
## 365/286 = 1.3
NR.F.pred <- 1.3 * 1/2 * (1.028) # average number of offspring https://animaldiversity.org/accounts/Macropus_rufogriseus/ + 286 day interbirth interval + 2.8% twinning (Catt, D. C. (1977). The breeding biology of Bennett’s wallaby(Macropus rufogriseus jruticus)in South Canterbury, New Zealand. New Zealand Journal of Zoology, 4(4), 401–411. doi:10.1080/03014223.1977.9517965)
NR.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
NR.alpha <- ceiling(exp(-1.34 + (0.214*log(NR.mass*1000))))
NR.alpha
NR.alpha <- 1 # Catt 1977 (14-17 months)

## define m function with age
## 0.667 females produce young in second year
## breeding from 14 months (1.167 years)
NR.m.vec <- c(0, 0.8*NR.F.pred, 0.9*NR.F.pred, rep(NR.F.pred, (NR.age.max+1-3)))
NR.m.sd.vec <- 0.05*NR.m.vec
NR.age.vec <- c(seq(0,16,1))
plot(NR.age.vec, NR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
NR.m.dat <- data.frame(NR.age.vec, NR.m.vec)
param.init <- c(0.7, 0.9, -20)
NR.fit.logp <- nls(NR.m.vec ~ a / (1+(NR.age.vec/b)^c), 
                   data = NR.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
NR.fit.logp.summ <- summary(NR.fit.logp)
plot(NR.age.vec, NR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
NR.age.vec.cont <- seq(0,max(NR.age.vec),1)
NR.pred.p.m <- coef(NR.fit.logp)[1] / (1+(NR.age.vec.cont/coef(NR.fit.logp)[2])^coef(NR.fit.logp)[3])
NR.pred.p.mm <- ifelse(NR.pred.p.m > 1, 1, NR.pred.p.m)
lines(NR.age.vec.cont, NR.pred.p.mm,lty=2,lwd=3,col="red")


## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
NR.s.tran <- ln.a.s + b.s*log(NR.mass*1000) + log(1)
NR.s.ad.yr <- 1.05 * exp(-exp(NR.s.tran))
NR.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.0*NR.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 6.5 # rate of mortality decline (also known as bt)
a2 <- 1 - NR.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.06e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.02 # rate of mortality increase
longev <- NR.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
NR.Sx <- c(0.997*NR.s.ad.yr, 1 - qx)
plot(x, NR.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
NR.s.sd.vec <- 0.05*NR.Sx

## create matrix
NR.popmat <- matrix(data = 0, nrow=NR.age.max+1, ncol=NR.age.max+1)
diag(NR.popmat[2:(NR.age.max+1),]) <- NR.Sx[-(NR.age.max+1)]
NR.popmat[NR.age.max+1,NR.age.max+1] <- NR.Sx[NR.age.max+1]
NR.popmat[1,] <- NR.pred.p.mm
colnames(NR.popmat) <- c(0:NR.age.max)
rownames(NR.popmat) <- c(0:NR.age.max)
NR.popmat.orig <- NR.popmat ## save original matrix

## matrix properties
max.lambda(NR.popmat.orig) ## 1-yr lambda
NR.lm.pred
max.r(NR.popmat.orig) # rate of population change, 1-yr
NR.ssd <- stable.stage.dist(NR.popmat.orig) ## stable stage distribution
plot(NR.age.vec, NR.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(NR.popmat.orig, NR.age.max) # reproductive value
NR.gen.l <- G.val(NR.popmat.orig, NR.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
NR.pop.found <- round(area*NR.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
NR.init.vec <- NR.ssd * NR.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*NR.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

NR.tot.F <- sum(NR.popmat.orig[1,])
NR.popmat <- NR.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
NR.n.mat <- matrix(0, nrow=NR.age.max+1,ncol=(t+1))
NR.n.mat[,1] <- NR.init.vec

## set up projection loop
for (i in 1:t) {
  NR.n.mat[,i+1] <- NR.popmat %*% NR.n.mat[,i]
}

NR.n.pred <- colSums(NR.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(NR.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
NR.K.max <- 1*NR.pop.found
NR.K.vec <- c(1, NR.K.max/2, 0.75*NR.K.max, NR.K.max) 
NR.red.vec <- c(1,0.90,0.78,0.625)
plot(NR.K.vec, NR.red.vec,pch=19,type="b")
NR.Kred.dat <- data.frame(NR.K.vec, NR.red.vec)

# logistic power function a/(1+(x/b)^c)
NR.param.init <- c(1, 2*NR.K.max, 2)
NR.fit.lp <- nls(NR.red.vec ~ a/(1+(NR.K.vec/b)^c), 
                 data = NR.Kred.dat,
                 algorithm = "port",
                 start = c(a = NR.param.init[1], b = NR.param.init[2], c = NR.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
NR.fit.lp.summ <- summary(NR.fit.lp)
plot(NR.K.vec, NR.red.vec, pch=19,xlab="N",ylab="reduction factor")
NR.K.vec.cont <- seq(1,2*NR.pop.found,1)
NR.pred.lp.fx <- coef(NR.fit.lp)[1]/(1+(NR.K.vec.cont/coef(NR.fit.lp)[2])^coef(NR.fit.lp)[3])
lines(NR.K.vec.cont, NR.pred.lp.fx, lty=3,lwd=3,col="red")

NR.a.lp <- coef(NR.fit.lp)[1]
NR.b.lp <- coef(NR.fit.lp)[2]
NR.c.lp <- coef(NR.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
NR.n.mat <- matrix(0, nrow=NR.age.max+1, ncol=(t+1))
NR.n.mat[,1] <- NR.init.vec
NR.popmat <- NR.popmat.orig

## set up projection loop
for (i in 1:t) {
  NR.totN.i <- sum(NR.n.mat[,i])
  NR.pred.red <- as.numeric(NR.a.lp/(1+(NR.totN.i/NR.b.lp)^NR.c.lp))
  diag(NR.popmat[2:(NR.age.max+1),]) <- (NR.Sx[-(NR.age.max+1)])*NR.pred.red
  NR.popmat[NR.age.max+1,NR.age.max+1] <- (NR.Sx[NR.age.max+1])*NR.pred.red
  NR.popmat[1,] <- NR.pred.p.mm
  NR.n.mat[,i+1] <- NR.popmat %*% NR.n.mat[,i]
}

NR.n.pred <- colSums(NR.n.mat)
plot(yrs, NR.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=NR.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

NR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
NR.s.arr <- NR.m.arr <- array(data=NA, dim=c(t+1, NR.age.max+1, iter))
NR.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  NR.popmat <- NR.popmat.orig
  
  NR.n.mat <- matrix(0, nrow=NR.age.max+1,ncol=(t+1))
  NR.n.mat[,1] <- NR.init.vec
  NR.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    NR.s.alpha <- estBetaParams(NR.Sx, NR.s.sd.vec^2)$alpha
    NR.s.beta <- estBetaParams(NR.Sx, NR.s.sd.vec^2)$beta
    NR.s.stoch <- rbeta(length(NR.s.alpha), NR.s.alpha, NR.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    NR.fert.stch <- rnorm(length(NR.popmat[,1]), NR.pred.p.mm, NR.m.sd.vec)
    NR.m.arr[i,,e] <- ifelse(NR.fert.stch < 0, 0, NR.fert.stch)
    
    NR.totN.i <- sum(NR.n.mat[,i], na.rm=T)
    NR.pred.red.vec[i] <- NR.a.lp/(1+(NR.totN.i/NR.b.lp)^NR.c.lp)
    
    diag(NR.popmat[2:(NR.age.max+1),]) <- (NR.s.stoch[-(NR.age.max+1)])*NR.pred.red.vec[i]
    NR.popmat[NR.age.max+1,NR.age.max+1] <- (NR.s.stoch[NR.age.max+1])*NR.pred.red.vec[i]
    NR.popmat[1,] <- NR.m.arr[i,,e]
    NR.n.mat[,i+1] <- NR.popmat %*% NR.n.mat[,i]
    
    NR.s.arr[i,,e] <- NR.s.stoch * NR.pred.red.vec[i]
    
  } # end i loop
  
  NR.n.sums.mat[e,] <- ((as.vector(colSums(NR.n.mat))/NR.pop.found))
  NR.pred.red.mn[e] <- mean(NR.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

NR.n.md <- apply(NR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
NR.n.up <- apply(NR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
NR.n.lo <- apply(NR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations



##############################
## GENYORNIS (newtonii) (GN)

# mass
GN.mass <- 200 # Genyornis newtoni (max) (Johnson et al. 2006)

## large, flightless birds (medium) log10(D) = 3.65 – 0.82×log10(M in g) (https://doi.org/10.1111/ecog.04917)
GN.D.pred <- (10^(3.65 - 0.82*(log10(GN.mass*1000))))/2
GN.D.pred

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
GN.age.max <- round(10^(0.89 + (0.13*log10(GN.mass*1000))), 0)
GN.age.max

## age vector
GN.age.vec <- 0:GN.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.35 - 0.17lnM (all birds)
GN.F.pred1 <- exp(2.35 - (0.17*log(GN.mass*1000)))/2 # divided by 2 for females
GN.F.pred1
GN.eggs <- GN.F.pred1/1.003673 * 11.8 / 2
GN.hatch <- (0.406 * 0.419)
GN.F.pred <- GN.eggs * GN.hatch


## age at primiparity
## alpha (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
GN.alpha <- ceiling(0.214*(GN.mass*1000)^0.303) # for birds | 4-5 years for ostriches (https://doi.org/10.1007/s11250-009-9428-2)
GN.alpha

## define m function with age
GN.m.vec <- c(rep(0, GN.alpha-1), rep(0.5*GN.F.pred, round(GN.alpha/2,0)), rep(GN.F.pred, (GN.age.max+1-((GN.alpha-1+round(GN.alpha/2,0))))))
GN.m.sd.vec <- 0.05*GN.m.vec
plot(GN.age.vec, GN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
GN.m.dat <- data.frame(GN.age.vec, GN.m.vec)
param.init <- c(0.5, 4, -4)
GN.fit.logp <- nls(GN.m.vec ~ a / (1+(GN.age.vec/b)^c), 
                   data = GN.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
GN.fit.logp.summ <- summary(GN.fit.logp)
plot(GN.age.vec, GN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
GN.age.vec.cont <- seq(0,max(GN.age.vec),1)
GN.pred.p.m <- coef(GN.fit.logp)[1] / (1+(GN.age.vec.cont/coef(GN.fit.logp)[2])^coef(GN.fit.logp)[3])
GN.pred.p.mm <- ifelse(GN.pred.p.m > 1, 1, GN.pred.p.m)
lines(GN.age.vec.cont, GN.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -1.78; b.s <- -0.21 # for birds
GN.s.tran <- ln.a.s + b.s*log(GN.mass*1000) + log(1)
GN.s.ad.yr <- exp(-exp(GN.s.tran))
GN.s.ad.yr

## calculate lmax from Dillingham et al. 2016 Ecol Appl 26:322-333
# lmax_DIM occurs when lmax.fun returns 0
lmax.fun <- function(lmax, alpha, s, ar.aT) {
  abs(lmax - exp(ar.aT*(alpha + s/(lmax-s) )^(-1) ))
}
# get.lmax calculates lmax_DIM for scalar alpha and scalar/vector s
get.lmax <- function(alpha, s, ar.aT) {
  # Fixed alpha, allows a single value or vector for s
  lmax <- rep(NA, length(s))
  for (i in 1:length(s)) {
    lmax[i] <- optimize(lmax.fun, c(1, 5), tol = 0.0001, alpha = alpha, s=s[i], ar.aT=ar.aT)$minimum
  }
  return(list(lmax=lmax, alpha=alpha, s=s, ar.aT=ar.aT))
}

GN.lm.pred <- (get.lmax(alpha=GN.alpha, s=GN.s.ad.yr, ar.aT=1.107))$lmax # for birds
GN.lm.pred
GN.rm.pred <- log(GN.lm.pred)

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.63*GN.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 0.3 # rate of mortality decline (also known as bt)
a2 <- 1 - GN.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.5e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.18 # rate of mortality increase
longev <- GN.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
GN.Sx <- c(0.63*GN.s.ad.yr, 1 - qx)
plot(x, GN.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
GN.s.sd.vec <- 0.05*GN.Sx

## create matrix
GN.popmat <- matrix(data = 0, nrow=GN.age.max+1, ncol=GN.age.max+1)
diag(GN.popmat[2:(GN.age.max+1),]) <- GN.Sx[-(GN.age.max+1)]
GN.popmat[GN.age.max+1,GN.age.max+1] <- GN.Sx[GN.age.max+1]
GN.popmat[1,] <- GN.pred.p.mm
colnames(GN.popmat) <- c(0:GN.age.max)
rownames(GN.popmat) <- c(0:GN.age.max)
GN.popmat.orig <- GN.popmat ## save original matrix

## matrix properties
max.lambda(GN.popmat.orig) ## 1-yr lambda
GN.lm.pred
max.r(GN.popmat.orig) # rate of population change, 1-yr
GN.ssd <- stable.stage.dist(GN.popmat.orig) ## stable stage distribution
plot(GN.age.vec, GN.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(GN.popmat.orig, GN.age.max) # reproductive value
GN.gen.l <- G.val(GN.popmat.orig, GN.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
GN.pop.found <- round(area*GN.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
GN.init.vec <- GN.ssd * GN.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*GN.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

GN.tot.F <- sum(GN.popmat.orig[1,])
GN.popmat <- GN.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
GN.n.mat <- matrix(0, nrow=GN.age.max+1,ncol=(t+1))
GN.n.mat[,1] <- GN.init.vec

## set up projection loop
for (i in 1:t) {
  GN.n.mat[,i+1] <- GN.popmat %*% GN.n.mat[,i]
}

GN.n.pred <- colSums(GN.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(GN.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
GN.K.max <- 1*GN.pop.found
GN.K.vec <- c(1, GN.K.max/2, 0.75*GN.K.max, GN.K.max) 
GN.red.vec <- c(1,0.995,0.98,0.9581)
plot(GN.K.vec, GN.red.vec,pch=19,type="b")
GN.Kred.dat <- data.frame(GN.K.vec, GN.red.vec)

# logistic power function a/(1+(x/b)^c)
GN.param.init <- c(1, 2*GN.K.max, 2)
GN.fit.lp <- nls(GN.red.vec ~ a/(1+(GN.K.vec/b)^c), 
                 data = GN.Kred.dat,
                 algorithm = "port",
                 start = c(a = GN.param.init[1], b = GN.param.init[2], c = GN.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
GN.fit.lp.summ <- summary(GN.fit.lp)
plot(GN.K.vec, GN.red.vec, pch=19,xlab="N",ylab="reduction factor")
GN.K.vec.cont <- seq(1,2*GN.pop.found,1)
GN.pred.lp.fx <- coef(GN.fit.lp)[1]/(1+(GN.K.vec.cont/coef(GN.fit.lp)[2])^coef(GN.fit.lp)[3])
lines(GN.K.vec.cont, GN.pred.lp.fx, lty=3,lwd=3,col="red")

GN.a.lp <- coef(GN.fit.lp)[1]
GN.b.lp <- coef(GN.fit.lp)[2]
GN.c.lp <- coef(GN.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
GN.n.mat <- matrix(0, nrow=GN.age.max+1, ncol=(t+1))
GN.n.mat[,1] <- GN.init.vec
GN.popmat <- GN.popmat.orig

## set up projection loop
for (i in 1:t) {
  GN.totN.i <- sum(GN.n.mat[,i])
  GN.pred.red <- as.numeric(GN.a.lp/(1+(GN.totN.i/GN.b.lp)^GN.c.lp))
  diag(GN.popmat[2:(GN.age.max+1),]) <- (GN.Sx[-(GN.age.max+1)])*GN.pred.red
  GN.popmat[GN.age.max+1,GN.age.max+1] <- (GN.Sx[GN.age.max+1])*GN.pred.red
  GN.popmat[1,] <- GN.pred.p.mm
  GN.n.mat[,i+1] <- GN.popmat %*% GN.n.mat[,i]
}

GN.n.pred <- colSums(GN.n.mat)
plot(yrs, GN.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=GN.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

GN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
GN.s.arr <- GN.m.arr <- array(data=NA, dim=c(t+1, GN.age.max+1, iter))
GN.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  GN.popmat <- GN.popmat.orig
  
  GN.n.mat <- matrix(0, nrow=GN.age.max+1,ncol=(t+1))
  GN.n.mat[,1] <- GN.init.vec
  GN.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    GN.s.alpha <- estBetaParams(GN.Sx, GN.s.sd.vec^2)$alpha
    GN.s.beta <- estBetaParams(GN.Sx, GN.s.sd.vec^2)$beta
    GN.s.stoch <- rbeta(length(GN.s.alpha), GN.s.alpha, GN.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    GN.fert.stch <- rnorm(length(GN.popmat[,1]), GN.pred.p.mm, GN.m.sd.vec)
    GN.m.arr[i,,e] <- ifelse(GN.fert.stch < 0, 0, GN.fert.stch)
    
    GN.totN.i <- sum(GN.n.mat[,i], na.rm=T)
    GN.pred.red.vec[i] <- GN.a.lp/(1+(GN.totN.i/GN.b.lp)^GN.c.lp)
    
    diag(GN.popmat[2:(GN.age.max+1),]) <- (GN.s.stoch[-(GN.age.max+1)])*GN.pred.red.vec[i]
    GN.popmat[GN.age.max+1,GN.age.max+1] <- (GN.s.stoch[GN.age.max+1])*GN.pred.red.vec[i]
    GN.popmat[1,] <- GN.m.arr[i,,e]
    GN.n.mat[,i+1] <- GN.popmat %*% GN.n.mat[,i]
    
    GN.s.arr[i,,e] <- GN.s.stoch * GN.pred.red.vec[i]
    
  } # end i loop
  
  GN.n.sums.mat[e,] <- ((as.vector(colSums(GN.n.mat))/GN.pop.found))
  GN.pred.red.mn[e] <- mean(GN.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

GN.n.md <- apply(GN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
GN.n.up <- apply(GN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
GN.n.lo <- apply(GN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations



##############################
## DROMAIUS (novaehollandiae) (DN)

# mass
DN.mass <- 55 # Dromaius novaehollandiae (Sales et al. 2007 Avian and Poultry Biology Reviews 18:1–20)

## large, flightless birds (medium) log10(D) = 3.65 – 0.82×log10(M in g) (https://doi.org/10.1111/ecog.04917)
DN.D.pred <- (10^(3.65 - 0.82*(log10(DN.mass*1000))))/2
DN.D.pred

## emu densities: 0.5 to 2.0 km-2 (Pople et al. 2000 Wildlife Research 27:269 doi:10.1071/wr99030)
## rhea K: 0.6 km-2 (https://ri.conicet.gov.ar/bitstream/handle/11336/7941/CONICET_Digital_Nro.10498.pdf?sequence=1&isAllowed=y)

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
DN.age.max <- round(10^(0.89 + (0.13*log10(DN.mass*1000))), 0)
DN.age.max
DN.age.max <- 17 # adjusted downward

## age vector
DN.age.vec <- 0:DN.age.max

## fertility
# 8-13 eggs/female ostrich (https://doi.org/10.1007/s11250-009-9428-2)
# 11.8 eggs/female ostrich (https://doi.org/10.1007/s11250-009-9331-x)
# ostrich nests with at least 1 egg hatch: 0.406 (https://doi.org/10.1007/s11250-009-9331-x)
# ostrich hatch rate: 0.419 ± 0.12 (https://doi.org/10.1007/s11250-009-9331-x)
# 11.8*0.406*0.419/2 = 1.004
# 6.7 eggs/clutch, 3.4 clutches/year # Sales 2007
# assuming ostrich hatching success: 6.7*3.4*.406*.419/2
DN.eggs <- 6.7 * 3.4 / 2 # /2 for females
DN.hatch <- 0.406 * 0.419
DN.F.pred <- DN.eggs * DN.hatch

## age at primiparity
## alpha (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
DN.alpha <- ceiling(0.214*(DN.mass*1000)^0.303) # for birds | 4-5 years for ostriches (https://doi.org/10.1007/s11250-009-9428-2)
DN.alpha
DN.alpha <- 3 # (http://www.veterinaryworld.org/Vol.2/November/Behavior%20of%20Emu%20bird%20(Dromaius%20novaehollandiae).pdf)

## define m function with age
DN.m.vec <- c(rep(0, DN.alpha-1), rep(0.5*DN.F.pred, round(DN.alpha/2,0)), rep(DN.F.pred, (DN.age.max+1-((DN.alpha-1+round(DN.alpha/2,0))))))
DN.m.sd.vec <- 0.05*DN.m.vec
plot(DN.age.vec, DN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
DN.m.dat <- data.frame(DN.age.vec, DN.m.vec)
param.init <- c(0.5, 4, -4)
DN.fit.logp <- nls(DN.m.vec ~ a / (1+(DN.age.vec/b)^c), 
                   data = DN.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DN.fit.logp.summ <- summary(DN.fit.logp)
plot(DN.age.vec, DN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
DN.age.vec.cont <- seq(0,max(DN.age.vec),1)
DN.pred.p.mm <- coef(DN.fit.logp)[1] / (1+(DN.age.vec.cont/coef(DN.fit.logp)[2])^coef(DN.fit.logp)[3])
#DN.pred.p.mm <- ifelse(DN.pred.p.m > 1, 1, DN.pred.p.m)
lines(DN.age.vec.cont, DN.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -1.78; b.s <- -0.21 # for birds
DN.s.tran <- ln.a.s + b.s*log(DN.mass*1000) + log(1)
DN.s.ad.yr <- exp(-exp(DN.s.tran))
DN.s.ad.yr

# Moa survival: 0: 0.25; 1: 0.74; 2: 0.75; 3: 0.83; 4: 0.93; 5+: 0.95 (Holdaway & Jacombe 2000 Science)

## calculate lmax from Dillingham et al. 2016 Ecol Appl 26:322-333
# lmax_DIM occurs when lmax.fun returns 0
lmax.fun <- function(lmax, alpha, s, ar.aT) {
  abs(lmax - exp(ar.aT*(alpha + s/(lmax-s) )^(-1) ))
}
# get.lmax calculates lmax_DIM for scalar alpha and scalar/vector s
get.lmax <- function(alpha, s, ar.aT) {
  # Fixed alpha, allows a single value or vector for s
  lmax <- rep(NA, length(s))
  for (i in 1:length(s)) {
    lmax[i] <- optimize(lmax.fun, c(1, 5), tol = 0.0001, alpha = alpha, s=s[i], ar.aT=ar.aT)$minimum
  }
  return(list(lmax=lmax, alpha=alpha, s=s, ar.aT=ar.aT))
}

DN.lm.pred <- (get.lmax(alpha=DN.alpha, s=DN.s.ad.yr, ar.aT=1.107))$lmax # for birds
DN.lm.pred
DN.rm.pred <- log(DN.lm.pred)

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.55*DN.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 0.105 # rate of mortality decline (also known as bt)
a2 <- 1 - DN.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 2.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.43 # rate of mortality increase
longev <- DN.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
DN.Sx <- c(0.55*DN.s.ad.yr, 1 - qx)
plot(x, DN.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
DN.s.sd.vec <- 0.05*DN.Sx

## create matrix
DN.popmat <- matrix(data = 0, nrow=DN.age.max+1, ncol=DN.age.max+1)
diag(DN.popmat[2:(DN.age.max+1),]) <- DN.Sx[-(DN.age.max+1)]
DN.popmat[DN.age.max+1,DN.age.max+1] <- DN.Sx[DN.age.max+1]
DN.popmat[1,] <- DN.pred.p.mm
colnames(DN.popmat) <- c(0:DN.age.max)
rownames(DN.popmat) <- c(0:DN.age.max)
DN.popmat.orig <- DN.popmat ## save original matrix

## matrix properties
max.lambda(DN.popmat.orig) ## 1-yr lambda
DN.lm.pred
max.r(DN.popmat.orig) # rate of population change, 1-yr
DN.ssd <- stable.stage.dist(DN.popmat.orig) ## stable stage distribution
plot(DN.age.vec, DN.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(DN.popmat.orig, DN.age.max) # reproductive value
DN.gen.l <- G.val(DN.popmat.orig, DN.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
DN.pop.found <- round(area*DN.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
DN.init.vec <- DN.ssd * DN.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*DN.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

DN.tot.F <- sum(DN.popmat.orig[1,])
DN.popmat <- DN.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
DN.n.mat <- matrix(0, nrow=DN.age.max+1,ncol=(t+1))
DN.n.mat[,1] <- DN.init.vec

## set up projection loop
for (i in 1:t) {
  DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]
}

DN.n.pred <- colSums(DN.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(DN.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
DN.K.max <- 1*DN.pop.found
DN.K.vec <- c(1, DN.K.max/2, 0.75*DN.K.max, DN.K.max) 
DN.red.vec <- c(1,0.97,0.932,0.885)
plot(DN.K.vec, DN.red.vec,pch=19,type="b")
DN.Kred.dat <- data.frame(DN.K.vec, DN.red.vec)

# logistic power function a/(1+(x/b)^c)
DN.param.init <- c(1, 2*DN.K.max, 2)
DN.fit.lp <- nls(DN.red.vec ~ a/(1+(DN.K.vec/b)^c), 
                 data = DN.Kred.dat,
                 algorithm = "port",
                 start = c(a = DN.param.init[1], b = DN.param.init[2], c = DN.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DN.fit.lp.summ <- summary(DN.fit.lp)
plot(DN.K.vec, DN.red.vec, pch=19,xlab="N",ylab="reduction factor")
DN.K.vec.cont <- seq(1,2*DN.pop.found,1)
DN.pred.lp.fx <- coef(DN.fit.lp)[1]/(1+(DN.K.vec.cont/coef(DN.fit.lp)[2])^coef(DN.fit.lp)[3])
lines(DN.K.vec.cont, DN.pred.lp.fx, lty=3,lwd=3,col="red")

DN.a.lp <- coef(DN.fit.lp)[1]
DN.b.lp <- coef(DN.fit.lp)[2]
DN.c.lp <- coef(DN.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
DN.n.mat <- matrix(0, nrow=DN.age.max+1, ncol=(t+1))
DN.n.mat[,1] <- DN.init.vec
DN.popmat <- DN.popmat.orig

## set up projection loop
for (i in 1:t) {
  DN.totN.i <- sum(DN.n.mat[,i])
  DN.pred.red <- as.numeric(DN.a.lp/(1+(DN.totN.i/DN.b.lp)^DN.c.lp))
  diag(DN.popmat[2:(DN.age.max+1),]) <- (DN.Sx[-(DN.age.max+1)])*DN.pred.red
  DN.popmat[DN.age.max+1,DN.age.max+1] <- (DN.Sx[DN.age.max+1])*DN.pred.red
  DN.popmat[1,] <- DN.pred.p.mm
  DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]
}

DN.n.pred <- colSums(DN.n.mat)
plot(yrs, DN.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=DN.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

DN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DN.s.arr <- DN.m.arr <- array(data=NA, dim=c(t+1, DN.age.max+1, iter))
DN.pred.red.mn <- rep(0, iter)

for (e in 1:iter) {
  DN.popmat <- DN.popmat.orig
  
  DN.n.mat <- matrix(0, nrow=DN.age.max+1,ncol=(t+1))
  DN.n.mat[,1] <- DN.init.vec
  DN.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    DN.s.alpha <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$alpha
    DN.s.beta <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$beta
    DN.s.stoch <- rbeta(length(DN.s.alpha), DN.s.alpha, DN.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    DN.fert.stch <- rnorm(length(DN.popmat[,1]), DN.pred.p.mm, DN.m.sd.vec)
    DN.m.arr[i,,e] <- ifelse(DN.fert.stch < 0, 0, DN.fert.stch)
    
    DN.totN.i <- sum(DN.n.mat[,i], na.rm=T)
    DN.pred.red.vec[i] <- DN.a.lp/(1+(DN.totN.i/DN.b.lp)^DN.c.lp)
    
    diag(DN.popmat[2:(DN.age.max+1),]) <- (DN.s.stoch[-(DN.age.max+1)])*DN.pred.red.vec[i]
    DN.popmat[DN.age.max+1,DN.age.max+1] <- (DN.s.stoch[DN.age.max+1])*DN.pred.red.vec[i]
    DN.popmat[1,] <- DN.m.arr[i,,e]
    DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]

    DN.s.arr[i,,e] <- DN.s.stoch * DN.pred.red.vec[i]
    
  } # end i loop
  
  DN.n.sums.mat[e,] <- ((as.vector(colSums(DN.n.mat))/DN.pop.found))
  DN.pred.red.mn[e] <- mean(DN.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

DN.n.md <- apply(DN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
DN.n.up <- apply(DN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DN.n.lo <- apply(DN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## ALECTURA (lathami) (AL)

# mass
AL.mass <- 2.2 # Alectura lathami (brush turkey) adult female (Jones, R.W.R.J. Dekker, C.S. Roselaar. The Megapodes: Megapodiidae, Oxford University Press, Oxford, New York, Tokyo (1995))
## for comparison, female malleefowl range from 1.68-2.235 kg (Leipoa ocellata) Priddel & Wheeler 2003 Wildl Res 30:451-464

## Omnivorous birds (Juanes 1986 Am Nat 128: 921-929)
AL.D.pred <- (10^(1.63 - 0.23*(log10(AL.mass*1000))))/2
AL.D.pred

# non-passerines (tmax=7.02*M^0.174) (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
AL.age.max <- round(7.02*(AL.mass*1000)^0.174,0)

## age vector
AL.age.vec <- 0:AL.age.max

## fertility
## 0.866 eggs hatch successfully (Jones 1988 Aust Ornithol)
## 0.90-0.97 mortality post-hatching to 1 year (Jones 1988)
## 0.88-1.0 mortality first 3 weeks post-hatching (Göth & Vogel 2002 Wildl Res 29:503-511)
## maximum clutch size (33-34) reported in Malleefowl (Frith 1962a Tne Mallee-fowl. Angus & Robertson, Sydney;
## Brickhill 1986 Breeding success of Malleefowl Leipoa ocellata in central New South Wales. Emu 87,42-4
## Booth 1987 Home range and hatching success of Malleefowl Leipoa ocellata Aust. Wild Res. 14,95-104)
## 29-30 eggs in two clutches (Weathers et al. 1990-Emu)
## mean 16.6 eggs/pair (Bode & Brennan 2011): N(mean=16.6, SD=4.5)
## Malleefowl mean clutch size = 15.6; 0.464 of eggs left nest (Brickhill 1987 Emu 87:42-45)
AL.eggs <- 16.6 / 2 # (for females only)
AL.hatch <- 0.866
AL.F.pred <- AL.eggs*AL.hatch
AL.F.pred

## age at primiparity
## alpha (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
AL.alpha <- ceiling(0.214*(AL.mass*1000)^0.303) # for birds
AL.alpha
AL.alpha <- 2 # Bode & Brennan 2011 (Frith 1959)

## define m function with age
AL.m.vec <- c(rep(0, AL.alpha-1), rep(0.5*AL.F.pred, round(AL.alpha/2,0)), 0.75*AL.F.pred, rep(AL.F.pred, (AL.age.max+1-((AL.alpha+round(AL.alpha/2,0))))))
AL.m.sd.vec <- 0.05*AL.m.vec
plot(AL.age.vec, AL.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# Bleasdale yield density function y = x(a + bx^θ)^(-1/θ)
AL.m.dat <- data.frame(AL.age.vec, AL.m.vec)
param.init <- c(8.8077e-03, 4.571e-04, 3.89)
AL.fit.byd <- nls(AL.m.vec ~ AL.age.vec * ((a + b*(AL.age.vec)^theta)^(-1/theta)), 
                   data = AL.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], theta = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
#AL.fit.byd.summ <- summary(AL.fit.byd)
plot(AL.age.vec, AL.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
AL.age.vec.cont <- seq(0,max(AL.age.vec),1)
AL.pred.p.mm <- AL.age.vec.cont * ((coef(AL.fit.byd)[1] + coef(AL.fit.byd)[2]*(AL.age.vec.cont)^coef(AL.fit.byd)[3])^(-1/coef(AL.fit.byd)[3]))
lines(AL.age.vec.cont, AL.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -1.78; b.s <- -0.21 # for birds
AL.s.tran <- ln.a.s + b.s*log(AL.mass*1000) + log(1)
AL.s.ad.yr <- exp(-exp(AL.s.tran))
AL.s.ad.yr

## calculate lmax from Dillingham et al. 2016 Ecol Appl 26:322-333
# lmax_DIM occurs when lmax.fun returns 0
lmax.fun <- function(lmax, alpha, s, ar.aT) {
  abs(lmax - exp(ar.aT*(alpha + s/(lmax-s) )^(-1) ))
}
# get.lmax calculates lmax_DIM for scalar alpha and scalar/vector s
get.lmax <- function(alpha, s, ar.aT) {
  # Fixed alpha, allows a single value or vector for s
  lmax <- rep(NA, length(s))
  for (i in 1:length(s)) {
    lmax[i] <- optimize(lmax.fun, c(1, 5), tol = 0.0001, alpha = alpha, s=s[i], ar.aT=ar.aT)$minimum
  }
  return(list(lmax=lmax, alpha=alpha, s=s, ar.aT=ar.aT))
}

AL.lm.pred <- (get.lmax(alpha=AL.alpha, s=AL.s.ad.yr, ar.aT=1.107))$lmax # for birds
AL.lm.pred
AL.rm.pred <- log(AL.lm.pred)

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.56*AL.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 0.3 # rate of mortality decline (also known as bt)
a2 <- 1 - AL.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 2.5e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.31 # rate of mortality increase
longev <- AL.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
AL.Sx <- c(0.1*AL.s.ad.yr, 1 - qx)
plot(x, AL.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
AL.s.sd.vec <- 0.05*AL.Sx

## create matrix
AL.popmat <- matrix(data = 0, nrow=AL.age.max+1, ncol=AL.age.max+1)
diag(AL.popmat[2:(AL.age.max+1),]) <- AL.Sx[-(AL.age.max+1)]
AL.popmat[AL.age.max+1,AL.age.max+1] <- AL.Sx[AL.age.max+1]
AL.popmat[1,] <- AL.pred.p.mm
colnames(AL.popmat) <- c(0:AL.age.max)
rownames(AL.popmat) <- c(0:AL.age.max)
AL.popmat.orig <- AL.popmat ## save original matrix

## matrix properties
max.lambda(AL.popmat.orig) ## 1-yr lambda
AL.lm.pred
max.r(AL.popmat.orig) # rate of population change, 1-yr
AL.ssd <- stable.stage.dist(AL.popmat.orig) ## stable stage distribution
plot(AL.age.vec, AL.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(AL.popmat.orig, AL.age.max) # reproductive value
AL.gen.l <- G.val(AL.popmat.orig, AL.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
AL.pop.found <- round(area*AL.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
AL.init.vec <- AL.ssd * AL.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*AL.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

AL.tot.F <- sum(AL.popmat.orig[1,])
AL.popmat <- AL.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
AL.n.mat <- matrix(0, nrow=AL.age.max+1,ncol=(t+1))
AL.n.mat[,1] <- AL.init.vec

## set up projection loop
for (i in 1:t) {
  AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
}

AL.n.pred <- colSums(AL.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(AL.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
AL.K.max <- 1*AL.pop.found
AL.K.vec <- c(1, AL.K.max/2, 0.75*AL.K.max, AL.K.max) 
AL.red.vec <- c(1,0.96,0.88,0.805)
plot(AL.K.vec, AL.red.vec,pch=19,type="b")
AL.Kred.dat <- data.frame(AL.K.vec, AL.red.vec)

# logistic power function a/(1+(x/b)^c)
AL.param.init <- c(1, 2*AL.K.max, 2)
AL.fit.lp <- nls(AL.red.vec ~ a/(1+(AL.K.vec/b)^c), 
                 data = AL.Kred.dat,
                 algorithm = "port",
                 start = c(a = AL.param.init[1], b = AL.param.init[2], c = AL.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
AL.fit.lp.summ <- summary(AL.fit.lp)
plot(AL.K.vec, AL.red.vec, pch=19,xlab="N",ylab="reduction factor")
AL.K.vec.cont <- seq(1,2*AL.pop.found,1)
AL.pred.lp.fx <- coef(AL.fit.lp)[1]/(1+(AL.K.vec.cont/coef(AL.fit.lp)[2])^coef(AL.fit.lp)[3])
lines(AL.K.vec.cont, AL.pred.lp.fx, lty=3,lwd=3,col="red")

AL.a.lp <- coef(AL.fit.lp)[1]
AL.b.lp <- coef(AL.fit.lp)[2]
AL.c.lp <- coef(AL.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
AL.n.mat <- matrix(0, nrow=AL.age.max+1, ncol=(t+1))
AL.n.mat[,1] <- AL.init.vec
AL.popmat <- AL.popmat.orig

## set up projection loop
for (i in 1:t) {
  AL.totN.i <- sum(AL.n.mat[,i])
  AL.pred.red <- as.numeric(AL.a.lp/(1+(AL.totN.i/AL.b.lp)^AL.c.lp))
  diag(AL.popmat[2:(AL.age.max+1),]) <- (AL.Sx[-(AL.age.max+1)])*AL.pred.red
  AL.popmat[AL.age.max+1,AL.age.max+1] <- (AL.Sx[AL.age.max+1])*AL.pred.red
  AL.popmat[1,] <- AL.pred.p.mm
  AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
}

AL.n.pred <- colSums(AL.n.mat)
plot(yrs, AL.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=AL.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

AL.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
AL.s.arr <- AL.m.arr <- array(data=NA, dim=c(t+1, AL.age.max+1, iter))
AL.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  AL.popmat <- AL.popmat.orig
  
  AL.n.mat <- matrix(0, nrow=AL.age.max+1,ncol=(t+1))
  AL.n.mat[,1] <- AL.init.vec
  AL.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    AL.s.alpha <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$alpha
    AL.s.beta <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$beta
    AL.s.stoch <- rbeta(length(AL.s.alpha), AL.s.alpha, AL.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    AL.fert.stch <- rnorm(length(AL.popmat[,1]), AL.pred.p.mm, AL.m.sd.vec)
    AL.m.arr[i,,e] <- ifelse(AL.fert.stch < 0, 0, AL.fert.stch)
    
    AL.totN.i <- sum(AL.n.mat[,i], na.rm=T)
    AL.pred.red.vec[i] <- AL.a.lp/(1+(AL.totN.i/AL.b.lp)^AL.c.lp)
    
    diag(AL.popmat[2:(AL.age.max+1),]) <- (AL.s.stoch[-(AL.age.max+1)])*AL.pred.red.vec[i]
    AL.popmat[AL.age.max+1,AL.age.max+1] <- (AL.s.stoch[AL.age.max+1])*AL.pred.red.vec[i]
    AL.popmat[1,] <- AL.m.arr[i,,e]
    AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
    
    AL.s.arr[i,,e] <- AL.s.stoch * AL.pred.red.vec[i]
    
  } # end i loop
  
  AL.n.sums.mat[e,] <- ((as.vector(colSums(AL.n.mat))/AL.pop.found))
  AL.pred.red.mn[e] <- mean(AL.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

AL.n.md <- apply(AL.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
AL.n.up <- apply(AL.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
AL.n.lo <- apply(AL.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## THYLACOLEO (TC)
## sources: Prowse et al. 2014 Ecology 95:693–702 http://dx.doi.org/10.1890/13-0746.1
##          Lachish et al. 2009 Journal of Animal Ecology 78:427–436

# mass
TC.mass <- 110 # Thylacoleo carnifex (Johnson et al. 2006)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
TC.rm.pred <- 10^(0.6914 - (0.2622*log10(TC.mass*1000)))
TC.rm.pred
TC.lm.pred <- exp(TC.rm.pred)
TC.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
TC.D.pred <- (10^(1.77 + (-1.02*log10(TC.mass))))/2 # divided by 2 for females only
TC.D.pred # animals/km2

# (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227); they said log10, but only makes sense at ln
lmu <- c(-2.4098957, -0.32328582, 2.5504842)
lDu <- c(4.3807473,2.2986877,-0.70201576)
lml <- c(-2.4042623,0.26574665,2.3397434)
lDl <- c(0.80104977, -1.3227539, -2.954604)
plot(lmu, lDu, pch=19, xlim=c(min(lml),max(lmu)), ylim=c(min(lDl),max(lDu)))
lmlDu.fit <- lm(lDu~lmu)
abline(lmlDu.fit, lty=2, col="red")
points(lml, lDl, pch=19)
lmlDl.fit <- lm(lDl~lml)
abline(lmlDl.fit, lty=2, col="red")
TC.D.pred.l <- exp(as.numeric(coef(lmlDl.fit)[1] + coef(lmlDl.fit)[2]*log(TC.mass)))/2
TC.D.pred.u <- exp(as.numeric(coef(lmlDu.fit)[1] + coef(lmlDu.fit)[2]*log(TC.mass)))/2
TC.D.pred.u

TC.D.pred <- TC.D.pred.u

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
TC.age.max <- round(10^(0.89 + (0.13*log10(TC.mass*1000))), 0)
TC.age.max
TC.age.max <- 17 # adjust for female lion max longevity

## age vector
TC.age.vec <- 0:TC.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
TC.F.pred <- exp(2.719 - (0.211*log(TC.mass*1000)))/2 # divided by 2 for females
TC.F.pred
TC.F.pred <- 0.5 # most large vombatiforms only have 1 offspring, which is reasonably close to the predicted value above

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
TC.alpha <- ceiling(exp(-1.34 + (0.214*log(TC.mass*1000))))
TC.alpha

## define m function with age
TC.m.vec <- c(rep(0, TC.alpha-1), rep(0.5*TC.F.pred, round(TC.alpha/2,0)), rep(TC.F.pred, (TC.age.max+1-((TC.alpha-1+round(TC.alpha/2,0))))))
TC.m.sd.vec <- 0.05*TC.m.vec
plot(TC.age.vec, TC.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
TC.m.dat <- data.frame(TC.age.vec, TC.m.vec)
param.init <- c(0.5, 4, -4)
TC.fit.logp <- nls(TC.m.vec ~ a / (1+(TC.age.vec/b)^c), 
                   data = TC.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TC.fit.logp.summ <- summary(TC.fit.logp)
plot(TC.age.vec, TC.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
TC.age.vec.cont <- seq(0,max(TC.age.vec),1)
TC.pred.p.m <- coef(TC.fit.logp)[1] / (1+(TC.age.vec.cont/coef(TC.fit.logp)[2])^coef(TC.fit.logp)[3])
TC.pred.p.mm <- ifelse(TC.pred.p.m > 1, 1, TC.pred.p.m)
lines(TC.age.vec.cont, TC.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
TC.s.tran <- ln.a.s + b.s*log(TC.mass*1000) + log(1)
TC.s.ad.yr <- exp(-exp(TC.s.tran))
TC.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 0.97 - (TC.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.4 # rate of mortality decline (also known as bt)
a2 <- 0.97 - TC.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.25e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.03 # rate of mortality increase
longev <- TC.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
TC.Sx <- c(1.02*TC.s.ad.yr, 1 - qx)
plot(x, TC.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
TC.s.sd.vec <- 0.05*TC.Sx

## create matrix
TC.popmat <- matrix(data = 0, nrow=TC.age.max+1, ncol=TC.age.max+1)
diag(TC.popmat[2:(TC.age.max+1),]) <- TC.Sx[-(TC.age.max+1)]
TC.popmat[TC.age.max+1,TC.age.max+1] <- TC.Sx[TC.age.max+1]
TC.popmat[1,] <- TC.pred.p.mm
colnames(TC.popmat) <- c(0:TC.age.max)
rownames(TC.popmat) <- c(0:TC.age.max)
TC.popmat.orig <- TC.popmat ## save original matrix

## matrix properties
max.lambda(TC.popmat.orig) ## 1-yr lambda
TC.lm.pred
max.r(TC.popmat.orig) # rate of population change, 1-yr
TC.ssd <- stable.stage.dist(TC.popmat.orig) ## stable stage distribution
plot(TC.age.vec, TC.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(TC.popmat.orig, TC.age.max) # reproductive value
TC.gen.l <- G.val(TC.popmat.orig, TC.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
TC.pop.found <- round(area*TC.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
TC.init.vec <- TC.ssd * TC.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*TC.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

TC.tot.F <- sum(TC.popmat.orig[1,])
TC.popmat <- TC.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
TC.n.mat <- matrix(0, nrow=TC.age.max+1,ncol=(t+1))
TC.n.mat[,1] <- TC.init.vec

## set up projection loop
for (i in 1:t) {
  TC.n.mat[,i+1] <- TC.popmat %*% TC.n.mat[,i]
}

TC.n.pred <- colSums(TC.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(TC.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
TC.K.max <- 1*TC.pop.found
TC.K.vec <- c(1, TC.K.max/2, 0.75*TC.K.max, TC.K.max) 
TC.red.vec <- c(1,0.95,0.85,0.804)
plot(TC.K.vec, TC.red.vec,pch=19,type="b")
TC.Kred.dat <- data.frame(TC.K.vec, TC.red.vec)

# logistic power function a/(1+(x/b)^c)
TC.param.init <- c(1, 2*TC.K.max, 2)
TC.fit.lp <- nls(TC.red.vec ~ a/(1+(TC.K.vec/b)^c), 
                 data = TC.Kred.dat,
                 algorithm = "port",
                 start = c(a = TC.param.init[1], b = TC.param.init[2], c = TC.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TC.fit.lp.summ <- summary(TC.fit.lp)
plot(TC.K.vec, TC.red.vec, pch=19,xlab="N",ylab="reduction factor")
TC.K.vec.cont <- seq(1,2*TC.pop.found,1)
TC.pred.lp.fx <- coef(TC.fit.lp)[1]/(1+(TC.K.vec.cont/coef(TC.fit.lp)[2])^coef(TC.fit.lp)[3])
lines(TC.K.vec.cont, TC.pred.lp.fx, lty=3,lwd=3,col="red")

TC.a.lp <- coef(TC.fit.lp)[1]
TC.b.lp <- coef(TC.fit.lp)[2]
TC.c.lp <- coef(TC.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
TC.n.mat <- matrix(0, nrow=TC.age.max+1, ncol=(t+1))
TC.n.mat[,1] <- TC.init.vec
TC.popmat <- TC.popmat.orig

## set up projection loop
for (i in 1:t) {
  TC.totN.i <- sum(TC.n.mat[,i])
  TC.pred.red <- as.numeric(TC.a.lp/(1+(TC.totN.i/TC.b.lp)^TC.c.lp))
  diag(TC.popmat[2:(TC.age.max+1),]) <- (TC.Sx[-(TC.age.max+1)])*TC.pred.red
  TC.popmat[TC.age.max+1,TC.age.max+1] <- (TC.Sx[TC.age.max+1])*TC.pred.red
  TC.popmat[1,] <- TC.pred.p.mm
  TC.n.mat[,i+1] <- TC.popmat %*% TC.n.mat[,i]
}

TC.n.pred <- colSums(TC.n.mat)
plot(yrs, TC.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=TC.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

TC.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TC.s.arr <- TC.m.arr <- array(data=NA, dim=c(t+1, TC.age.max+1, iter))
TC.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  TC.popmat <- TC.popmat.orig
  
  TC.n.mat <- matrix(0, nrow=TC.age.max+1,ncol=(t+1))
  TC.n.mat[,1] <- TC.init.vec
  TC.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    TC.s.alpha <- estBetaParams(TC.Sx, TC.s.sd.vec^2)$alpha
    TC.s.beta <- estBetaParams(TC.Sx, TC.s.sd.vec^2)$beta
    TC.s.stoch <- rbeta(length(TC.s.alpha), TC.s.alpha, TC.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    TC.fert.stch <- rnorm(length(TC.popmat[,1]), TC.pred.p.mm, TC.m.sd.vec)
    TC.m.arr[i,,e] <- ifelse(TC.fert.stch < 0, 0, TC.fert.stch)
    
    TC.totN.i <- sum(TC.n.mat[,i], na.rm=T)
    TC.pred.red.vec[i] <- TC.a.lp/(1+(TC.totN.i/TC.b.lp)^TC.c.lp)
    
    diag(TC.popmat[2:(TC.age.max+1),]) <- (TC.s.stoch[-(TC.age.max+1)])*TC.pred.red.vec[i]
    TC.popmat[TC.age.max+1,TC.age.max+1] <- (TC.s.stoch[TC.age.max+1])*TC.pred.red.vec[i]
    TC.popmat[1,] <- TC.m.arr[i,,e]
    TC.n.mat[,i+1] <- TC.popmat %*% TC.n.mat[,i]

    TC.s.arr[i,,e] <- TC.s.stoch * TC.pred.red.vec[i]
    
  } # end i loop
  
  TC.n.sums.mat[e,] <- ((as.vector(colSums(TC.n.mat))/TC.pop.found))
  TC.pred.red.mn[e] <- mean(TC.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

TC.n.md <- apply(TC.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
TC.n.up <- apply(TC.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TC.n.lo <- apply(TC.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## THYLACINUS (TH)
## sources: Prowse et al. 2014 Ecology 95:693–702 http://dx.doi.org/10.1890/13-0746.1
##          Lachish et al. 2009 Journal of Animal Ecology 78:427–436

# mass
TH.mass <- 20 # kg (Jones and Stoddart 1998 Journal of Zoology 246:239-246; Lachish et al. 2009)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
TH.rm.pred <- 10^(0.6914 - (0.2622*log10(TH.mass*1000)))
TH.rm.pred
TH.lm.pred <- exp(TH.rm.pred)
TH.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log(y) = β0 + β1 log(x) (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227)
## for carnivores: intercept: 1.77; slope = −1.02 (upper)
## for carnivores: intercept: -2.37; slope = −1.13 (upper)
TH.D.pred <- (10^(1.77 + (-1.02*log10(TH.mass))))/2 # divided by 2 for females only
TH.D.pred # animals/km2

# (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227); they said log10, but only makes sense at ln
lmu <- c(-2.4098957, -0.32328582, 2.5504842)
lDu <- c(4.3807473,2.2986877,-0.70201576)
lml <- c(-2.4042623,0.26574665,2.3397434)
lDl <- c(0.80104977, -1.3227539, -2.954604)
plot(lmu, lDu, pch=19, xlim=c(min(lml),max(lmu)), ylim=c(min(lDl),max(lDu)))
lmlDu.fit <- lm(lDu~lmu)
abline(lmlDu.fit, lty=2, col="red")
points(lml, lDl, pch=19)
lmlDl.fit <- lm(lDl~lml)
abline(lmlDl.fit, lty=2, col="red")
TH.D.pred.l <- exp(as.numeric(coef(lmlDl.fit)[1] + coef(lmlDl.fit)[2]*log(TH.mass)))/2
TH.D.pred.u <- exp(as.numeric(coef(lmlDu.fit)[1] + coef(lmlDu.fit)[2]*log(TH.mass)))/2
TH.D.pred.u

TH.D.pred <- TH.D.pred.u

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
TH.age.max <- round(10^(0.89 + (0.13*log10(TH.mass*1000))), 0)
TH.age.max
TH.age.max <- 10 # http://www.esapubs.org/archive/ecol/E095/057/appendix-A.php (downward for Dasyurid shorter lifespan)
## ~ matches lifespan of 9 (Collins 1973 Monotremes and marsupials: a reference for zoological institutions  Washington: Smithsonian Institute Press)

## age vector
TH.age.vec <- 0:TH.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
TH.F.pred <- exp(2.719 - (0.211*log(TH.mass*1000)))/2 # divided by 2 for females
TH.F.pred
TH.F.pred <- 3.42*0.91/2 # http://www.esapubs.org/archive/ecol/E095/057/appendix-A.php
## matches litter size = 3 (Rounsevell & Mooney 1995 Thylacine In The Australian Museum complete book of Australian Mammals (ed. R. Strahan) pp. 164-165. Sydney: Reed Books)

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
TH.alpha <- ceiling(exp(-1.34 + (0.214*log(TH.mass*1000))))
TH.alpha
TH.alpha <- 1

## define m function with age
TH.m.vec <- c(0, 0.075*3.42/2, rep(TH.F.pred,9))
TH.m.sd.vec <- 0.05*TH.m.vec
plot(TH.age.vec, TH.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
TH.m.dat <- data.frame(TH.age.vec, TH.m.vec)
param.init <- c(1.55, 1, -30)
TH.fit.logp <- nls(TH.m.vec ~ a / (1+(TH.age.vec/b)^c), 
                   data = TH.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TH.fit.logp.summ <- summary(TH.fit.logp)
plot(TH.age.vec, TH.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
TH.age.vec.cont <- seq(0,max(TH.age.vec),1)
TH.pred.p.m <- coef(TH.fit.logp)[1] / (1+(TH.age.vec.cont/coef(TH.fit.logp)[2])^coef(TH.fit.logp)[3])
#TH.pred.p.mm <- ifelse(TH.pred.p.m > 1, 1, TH.pred.p.m)
TH.pred.p.mm <- TH.pred.p.m
lines(TH.age.vec.cont, TH.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
TH.s.tran <- ln.a.s + b.s*log(TH.mass*1000) + log(1)
TH.s.ad.yr <- exp(-exp(TH.s.tran))
TH.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.80*TH.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.1 # rate of mortality decline (also known as bt)
a2 <- 1 - TH.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.5e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.7 # rate of mortality increase
longev <- TH.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
TH.Sx <- c(0.80*TH.s.ad.yr, 1 - qx)
plot(x, TH.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
TH.s.sd.vec <- 0.05*TH.Sx

## create matrix
TH.popmat <- matrix(data = 0, nrow=TH.age.max+1, ncol=TH.age.max+1)
diag(TH.popmat[2:(TH.age.max+1),]) <- TH.Sx[-(TH.age.max+1)]
TH.popmat[TH.age.max+1,TH.age.max+1] <- 0 # catastrophic mortality at longevity TH.Sx[TH.age.max+1]
TH.popmat[1,] <- TH.pred.p.mm
colnames(TH.popmat) <- c(0:TH.age.max)
rownames(TH.popmat) <- c(0:TH.age.max)
TH.popmat.orig <- TH.popmat ## save original matrix

## matrix properties
max.lambda(TH.popmat.orig) ## 1-yr lambda
TH.lm.pred
max.r(TH.popmat.orig) # rate of population change, 1-yr
TH.ssd <- stable.stage.dist(TH.popmat.orig) ## stable stage distribution
plot(TH.age.vec, TH.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(TH.popmat.orig, TH.age.max) # reproductive value
TH.gen.l <- G.val(TH.popmat.orig, TH.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
TH.pop.found <- round(area*TH.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
TH.init.vec <- TH.ssd * TH.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*TH.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

TH.tot.F <- sum(TH.popmat.orig[1,])
TH.popmat <- TH.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
TH.n.mat <- matrix(0, nrow=TH.age.max+1,ncol=(t+1))
TH.n.mat[,1] <- TH.init.vec

## set up projection loop
for (i in 1:t) {
  TH.n.mat[,i+1] <- TH.popmat %*% TH.n.mat[,i]
}

TH.n.pred <- colSums(TH.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(TH.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
TH.K.max <- 1*TH.pop.found
TH.K.vec <- c(1, TH.K.max/2, 0.75*TH.K.max, TH.K.max) 
TH.red.vec <- c(1,0.94,0.71,0.645)
plot(TH.K.vec, TH.red.vec,pch=19,type="b")
TH.Kred.dat <- data.frame(TH.K.vec, TH.red.vec)

# logistic power function a/(1+(x/b)^c)
TH.param.init <- c(1, 2*TH.K.max, 2)
TH.fit.lp <- nls(TH.red.vec ~ a/(1+(TH.K.vec/b)^c), 
                 data = TH.Kred.dat,
                 algorithm = "port",
                 start = c(a = TH.param.init[1], b = TH.param.init[2], c = TH.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TH.fit.lp.summ <- summary(TH.fit.lp)
plot(TH.K.vec, TH.red.vec, pch=19,xlab="N",ylab="reduction factor")
TH.K.vec.cont <- seq(1,2*TH.pop.found,1)
TH.pred.lp.fx <- coef(TH.fit.lp)[1]/(1+(TH.K.vec.cont/coef(TH.fit.lp)[2])^coef(TH.fit.lp)[3])
lines(TH.K.vec.cont, TH.pred.lp.fx, lty=3,lwd=3,col="red")

TH.a.lp <- coef(TH.fit.lp)[1]
TH.b.lp <- coef(TH.fit.lp)[2]
TH.c.lp <- coef(TH.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
TH.n.mat <- matrix(0, nrow=TH.age.max+1, ncol=(t+1))
TH.n.mat[,1] <- TH.init.vec
TH.popmat <- TH.popmat.orig

## set up projection loop
for (i in 1:t) {
  TH.totN.i <- sum(TH.n.mat[,i])
  TH.pred.red <- as.numeric(TH.a.lp/(1+(TH.totN.i/TH.b.lp)^TH.c.lp))
  diag(TH.popmat[2:(TH.age.max+1),]) <- (TH.Sx[-(TH.age.max+1)])*TH.pred.red
  TH.popmat[TH.age.max+1,TH.age.max+1] <- 0
  TH.popmat[1,] <- TH.pred.p.mm
  TH.n.mat[,i+1] <- TH.popmat %*% TH.n.mat[,i]
}

TH.n.pred <- colSums(TH.n.mat)
plot(yrs, TH.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=TH.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

TH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TH.s.arr <- TH.m.arr <- array(data=NA, dim=c(t+1, TH.age.max+1, iter))
TH.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  TH.popmat <- TH.popmat.orig
  
  TH.n.mat <- matrix(0, nrow=TH.age.max+1,ncol=(t+1))
  TH.n.mat[,1] <- TH.init.vec
  TH.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    TH.s.alpha <- estBetaParams(TH.Sx, TH.s.sd.vec^2)$alpha
    TH.s.beta <- estBetaParams(TH.Sx, TH.s.sd.vec^2)$beta
    TH.s.stoch <- rbeta(length(TH.s.alpha), TH.s.alpha, TH.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    TH.fert.stch <- rnorm(length(TH.popmat[,1]), TH.pred.p.mm, TH.m.sd.vec)
    TH.m.arr[i,,e] <- ifelse(TH.fert.stch < 0, 0, TH.fert.stch)
    
    TH.totN.i <- sum(TH.n.mat[,i], na.rm=T)
    TH.pred.red.vec[i] <- TH.a.lp/(1+(TH.totN.i/TH.b.lp)^TH.c.lp)
    
    diag(TH.popmat[2:(TH.age.max+1),]) <- (TH.s.stoch[-(TH.age.max+1)])*TH.pred.red.vec[i]
    TH.popmat[TH.age.max+1,TH.age.max+1] <- 0
    TH.popmat[1,] <- TH.m.arr[i,,e]
    TH.n.mat[,i+1] <- TH.popmat %*% TH.n.mat[,i]

    TH.s.arr[i,,e] <- TH.s.stoch * TH.pred.red.vec[i]
    
  } # end i loop
  
  TH.n.sums.mat[e,] <- ((as.vector(colSums(TH.n.mat))/TH.pop.found))
  TH.pred.red.mn[e] <- mean(TH.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

TH.n.md <- apply(TH.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
TH.n.up <- apply(TH.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TH.n.lo <- apply(TH.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## SARCOPHILUS (harrisi) (SH)
## sources: Bradshaw & Brook (2005) Ecography

# mass
SH.mass <- 6.1 # Sarcophilus (female; Guiler 1978 in Bradshaw & Brook 2005)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
SH.rm.pred <- 10^(0.6914 - (0.2622*log10(SH.mass*1000)))
SH.rm.pred
SH.lm.pred <- exp(SH.rm.pred)
SH.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
## log(y) = β0 + β1 log(x) (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227)
## for carnivores: intercept: 1.77; slope = −1.02 (upper)
## for carnivores: intercept: -2.37; slope = −1.13 (upper)
SH.D.pred <- (10^(1.77 + (-1.02*log10(SH.mass))))/2 # divided by 2 for females only
SH.D.pred # animals/km2

# (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227); they said log10, but only makes sense at ln
lmu <- c(-2.4098957, -0.32328582, 2.5504842)
lDu <- c(4.3807473,2.2986877,-0.70201576)
lml <- c(-2.4042623,0.26574665,2.3397434)
lDl <- c(0.80104977, -1.3227539, -2.954604)
plot(lmu, lDu, pch=19, xlim=c(min(lml),max(lmu)), ylim=c(min(lDl),max(lDu)))
lmlDu.fit <- lm(lDu~lmu)
abline(lmlDu.fit, lty=2, col="red")
points(lml, lDl, pch=19)
lmlDl.fit <- lm(lDl~lml)
abline(lmlDl.fit, lty=2, col="red")
SH.D.pred.l <- exp(as.numeric(coef(lmlDl.fit)[1] + coef(lmlDl.fit)[2]*log(SH.mass)))/2
SH.D.pred.u <- exp(as.numeric(coef(lmlDu.fit)[1] + coef(lmlDu.fit)[2]*log(SH.mass)))/2
SH.D.pred.u

## 0.3–0.7 individuals per km² (http://www.environment.gov.au/cgi-bin/sprat/public/publicspecies.pl?taxon_id=299)
SH.D.pred <- SH.D.pred.u

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
SH.age.max <- round(10^(0.89 + (0.13*log10(SH.mass*1000))), 0)
SH.age.max
SH.age.max <- 5 # (downward for Dasyurid shorter lifespan)

## age vector
SH.age.vec <- 0:SH.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
SH.F.pred <- exp(2.719 - (0.211*log(SH.mass*1000)))/2 # divided by 2 for females
SH.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
SH.alpha <- ceiling(exp(-1.34 + (0.214*log(SH.mass*1000))))
SH.alpha

## define m function with age
## sex ratio
SH.sr <- 0.5362

## no. pouch young
SH.npy <- c(0.061, 2, 3.6, 3.6, 3.6, 2)

## pouch prop female
SH.pypf <- mean(c(0.48,0.58,0.57,0.56,0.5))

## prop females carrying young
SH.pfcy <- 0.8

## survivals
SH.Sx <- c(0.3983, 0.82, 0.82, 0.82, 0.82, 0.2711)
plot(SH.age.vec, SH.Sx, pch=19, type="b", xlab="age (years)", ylab="Sx")
SH.s.sd.vec <- 0.05*SH.Sx

## pre-breeding design with 0-1 survival in first row
SH.m.vec <- c(SH.Sx[1]*SH.npy[1], SH.Sx[1]*SH.npy[2]*SH.pypf*SH.pfcy, SH.Sx[1]*SH.npy[3]*SH.pypf*SH.pfcy, SH.Sx[1]*SH.npy[4]*SH.pypf*SH.pfcy, SH.Sx[1]*SH.npy[5]*SH.pypf*SH.pfcy, SH.Sx[1]*SH.npy[6]*SH.pypf*SH.pfcy)
SH.m.sd.vec <- 0.05*SH.m.vec
plot(SH.age.vec, SH.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

## create matrix
SH.popmat <- matrix(data = 0, nrow=SH.age.max+1, ncol=SH.age.max+1)
diag(SH.popmat[2:(SH.age.max+1),]) <- SH.Sx[-1]
SH.popmat[SH.age.max+1,SH.age.max+1] <- 0 # catastrophic mortality at longevity SH.Sx[SH.age.max+1]
SH.popmat[1,] <- SH.m.vec
colnames(SH.popmat) <- c(0:SH.age.max)
rownames(SH.popmat) <- c(0:SH.age.max)
SH.popmat.orig <- SH.popmat ## save original matrix

## matrix properties
max.lambda(SH.popmat.orig) ## 1-yr lambda
SH.lm.pred
max.r(SH.popmat.orig) # rate of population change, 1-yr
SH.ssd <- stable.stage.dist(SH.popmat.orig) ## stable stage distribution
plot(SH.age.vec, SH.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(SH.popmat.orig, SH.age.max) # reproductive value
SH.gen.l <- G.val(SH.popmat.orig, SH.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
SH.pop.found <- round(area*SH.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
SH.init.vec <- SH.ssd * SH.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*SH.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

SH.tot.F <- sum(SH.popmat.orig[1,])
SH.popmat <- SH.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
SH.n.mat <- matrix(0, nrow=SH.age.max+1,ncol=(t+1))
SH.n.mat[,1] <- SH.init.vec

## set up projection loop
for (i in 1:t) {
  SH.n.mat[,i+1] <- SH.popmat %*% SH.n.mat[,i]
}

SH.n.pred <- colSums(SH.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(SH.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
SH.K.max <- 1*SH.pop.found
SH.K.vec <- c(1, SH.K.max/2, 0.75*SH.K.max, SH.K.max) 
SH.red.vec <- c(1,0.96,0.943,0.871)
plot(SH.K.vec, SH.red.vec,pch=19,type="b")
SH.Kred.dat <- data.frame(SH.K.vec, SH.red.vec)

# logistic power function a/(1+(x/b)^c)
SH.param.init <- c(1, 2*SH.K.max, 2)
SH.fit.lp <- nls(SH.red.vec ~ a/(1+(SH.K.vec/b)^c), 
                 data = SH.Kred.dat,
                 algorithm = "port",
                 start = c(a = SH.param.init[1], b = SH.param.init[2], c = SH.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
SH.fit.lp.summ <- summary(SH.fit.lp)
plot(SH.K.vec, SH.red.vec, pch=19,xlab="N",ylab="reduction factor")
SH.K.vec.cont <- seq(1,2*SH.pop.found,1)
SH.pred.lp.fx <- coef(SH.fit.lp)[1]/(1+(SH.K.vec.cont/coef(SH.fit.lp)[2])^coef(SH.fit.lp)[3])
lines(SH.K.vec.cont, SH.pred.lp.fx, lty=3,lwd=3,col="red")

SH.a.lp <- coef(SH.fit.lp)[1]
SH.b.lp <- coef(SH.fit.lp)[2]
SH.c.lp <- coef(SH.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
SH.n.mat <- matrix(0, nrow=SH.age.max+1, ncol=(t+1))
SH.n.mat[,1] <- SH.init.vec
SH.popmat <- SH.popmat.orig

## set up projection loop
for (i in 1:t) {
  SH.totN.i <- sum(SH.n.mat[,i])
  SH.pred.red <- as.numeric(SH.a.lp/(1+(SH.totN.i/SH.b.lp)^SH.c.lp))
  diag(SH.popmat[2:(SH.age.max+1),]) <- (SH.Sx[-1])*SH.pred.red
  SH.popmat[SH.age.max+1,SH.age.max+1] <- 0
  SH.popmat[1,] <- SH.m.vec
  SH.n.mat[,i+1] <- SH.popmat %*% SH.n.mat[,i]
}

SH.n.pred <- colSums(SH.n.mat)
plot(yrs, SH.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=SH.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

SH.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
SH.s.arr <- SH.m.arr <- array(data=NA, dim=c(t+1, SH.age.max+1, iter))
SH.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  SH.popmat <- SH.popmat.orig
  
  SH.n.mat <- matrix(0, nrow=SH.age.max+1,ncol=(t+1))
  SH.n.mat[,1] <- SH.init.vec
  SH.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    SH.s.alpha <- estBetaParams(SH.Sx, SH.s.sd.vec^2)$alpha
    SH.s.beta <- estBetaParams(SH.Sx, SH.s.sd.vec^2)$beta
    SH.s.stoch <- rbeta(length(SH.s.alpha), SH.s.alpha, SH.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    SH.fert.stch <- rnorm(length(SH.popmat[,1]), SH.m.vec, SH.m.sd.vec)
    SH.m.arr[i,,e] <- ifelse(SH.fert.stch < 0, 0, SH.fert.stch)
    
    SH.totN.i <- sum(SH.n.mat[,i], na.rm=T)
    SH.pred.red.vec[i] <- SH.a.lp/(1+(SH.totN.i/SH.b.lp)^SH.c.lp)
    
    diag(SH.popmat[2:(SH.age.max+1),]) <- (SH.s.stoch[-1])*SH.pred.red.vec[i]
    SH.popmat[SH.age.max+1,SH.age.max+1] <- 0
    SH.popmat[1,] <- SH.m.arr[i,,e]
    SH.n.mat[,i+1] <- SH.popmat %*% SH.n.mat[,i]

    SH.s.arr[i,,e] <- SH.s.stoch * SH.pred.red.vec[i]
    
  } # end i loop
  
  SH.n.sums.mat[e,] <- ((as.vector(colSums(SH.n.mat))/SH.pop.found))
  SH.pred.red.mn[e] <- mean(SH.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

SH.n.md <- apply(SH.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
SH.n.up <- apply(SH.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
SH.n.lo <- apply(SH.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations



##############################
## DASYURUS (maculatus) (DM)

# mass
DM.mass <- 1.68 ## 1.7 kg (Glen 2008-AJZ)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
DM.rm.pred <- 10^(0.6914 - (0.2622*log10(DM.mass*1000)))
DM.rm.pred
DM.lm.pred <- exp(DM.rm.pred)
DM.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
## log(y) = β0 + β1 log(x) (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227)
## for carnivores: intercept: 1.77; slope = −1.02 (upper)
## for carnivores: intercept: -2.37; slope = −1.13 (lower)
DM.D.pred <- (10^(1.77 + (-1.02*log10(DM.mass))))/2 # divided by 2 for females only
DM.D.pred # animals/km2


# (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13227); they said log10, but only makes sense at ln
lmu <- c(-2.4098957, -0.32328582, 2.5504842)
lDu <- c(4.3807473,2.2986877,-0.70201576)
lml <- c(-2.4042623,0.26574665,2.3397434)
lDl <- c(0.80104977, -1.3227539, -2.954604)
plot(lmu, lDu, pch=19, xlim=c(min(lml),max(lmu)), ylim=c(min(lDl),max(lDu)))
lmlDu.fit <- lm(lDu~lmu)
abline(lmlDu.fit, lty=2, col="red")
points(lml, lDl, pch=19)
lmlDl.fit <- lm(lDl~lml)
abline(lmlDl.fit, lty=2, col="red")
DM.D.pred.l <- exp(as.numeric(coef(lmlDl.fit)[1] + coef(lmlDl.fit)[2]*log(DM.mass)))/2
DM.D.pred.u <- exp(as.numeric(coef(lmlDu.fit)[1] + coef(lmlDu.fit)[2]*log(DM.mass)))/2
DM.D.pred.u

## in an area of unusually high abundance in north-eastern New South Wales, where density is conservatively estimated at 0.3 km−2 (0.1-0.5)
## 3.5 km-2 (Johnson pers. obs. in Fisher et al. 2001)
## Glen, A.S. (2008) Population attributes of the spotted-tailed quoll (Dasyurus maculatus) in north-eastern New South Wales. Australian Journal of Zoology, 56, 137-142.
DM.D.pred <- DM.D.pred.u # better match to Johnson pers. obs.

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
DM.age.max <- round(10^(0.89 + (0.13*log10(DM.mass*1000))), 0)
DM.age.max # way, way overestimated
DM.age.max <- 4 # Glen & Dickman 2013; Moro et al. 2019 & Cremona et al. 2104 (D. hallucatus)
## changed to 4 to account for a few females persisting into fourth year, but very low survival

## age vector
DM.age.vec <- 0:DM.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
## mean litter size = 4.9 (Glen, A.S. 2008 Population attributes of the spotted-tailed quoll (Dasyurus maculatus) in north-eastern New South Wales. Australian Journal of Zoology, 56, 137-142)
## litter size = 5.5 in Tasmania (Green & Scarborough)
## take as '5'?
## % females breeding = .643 Glen & Dickman (or 0.91 for D maculatus - Cremona et al. 214)
5*.643/2
DM.F.pred <- exp(2.719 - (0.211*log(DM.mass*1000)))/2 # divided by 2 for females
DM.F.pred # (almost spot on)

## age at primiparity
DM.alpha <- 1 # Glen & Dickman 2013; Moro et al. 2019 (D. hallucatus)

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
DM.s.tran <- ln.a.s + b.s*log(DM.mass*1000) + log(1)
DM.s.ad.yr <- exp(-exp(DM.s.tran))
DM.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.99*DM.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.8 # rate of mortality decline (also known as bt)
a2 <- 1 - (0.99*DM.s.ad.yr) # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.6e-04 # initial adult mortality rate (also known as βt)
longev <- DM.age.max
b3 <- 2.5 # rate of mortality increase
x <- seq(0,longev,1) # age vector
plot(x,log(h.x),pch=19,type="l")
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
DM.Sx <- c(0.95*DM.s.ad.yr, 1 - qx)
plot(x, DM.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
DM.s.sd.vec <- 0.05*DM.Sx

## pre-breeding design with 0-1 survival in first row
DM.m.vec <- c(0.4*DM.F.pred, rep(DM.F.pred,3), 0.95*DM.F.pred) # allowed some breeding in < 1 year because sexual maturing < 12 months in D maculatus (Moro et al. 2019)
DM.m.sd.vec <- 0.05*DM.m.vec
plot(DM.age.vec, DM.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

## create matrix
DM.popmat <- matrix(data = 0, nrow=DM.age.max+1, ncol=DM.age.max+1)
diag(DM.popmat[2:(DM.age.max+1),]) <- DM.Sx[-1]
DM.popmat[DM.age.max+1,DM.age.max+1] <- 0 # catastrophic mortality at longevity DM.Sx[DM.age.max+1]
DM.popmat[1,] <- DM.m.vec
colnames(DM.popmat) <- c(0:DM.age.max)
rownames(DM.popmat) <- c(0:DM.age.max)
DM.popmat.orig <- DM.popmat ## save original matrix

## matrix properties
max.lambda(DM.popmat.orig) ## 1-yr lambda
DM.lm.pred
max.r(DM.popmat.orig) # rate of population change, 1-yr
DM.ssd <- stable.stage.dist(DM.popmat.orig) ## stable stage distribution
plot(DM.age.vec, DM.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(DM.popmat.orig, DM.age.max) # reproductive value
DM.gen.l <- G.val(DM.popmat.orig, DM.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
DM.pop.found <- round(area*DM.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
DM.init.vec <- DM.ssd * DM.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*DM.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

DM.tot.F <- sum(DM.popmat.orig[1,])
DM.popmat <- DM.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
DM.n.mat <- matrix(0, nrow=DM.age.max+1,ncol=(t+1))
DM.n.mat[,1] <- DM.init.vec

## set up projection loop
for (i in 1:t) {
  DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
}

DM.n.pred <- colSums(DM.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(DM.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
DM.K.max <- 1*DM.pop.found
DM.K.vec <- c(1, 0.25*DM.K.max, DM.K.max/2, 0.75*DM.K.max, DM.K.max) 
DM.red.vec <- c(1,0.96,0.76,0.45,0.20)
plot(DM.K.vec, DM.red.vec,pch=19,type="b")
DM.Kred.dat <- data.frame(DM.K.vec, DM.red.vec)

# logistic power function a/(1+(x/b)^c)
DM.param.init <- c(1, DM.K.max, 2)
DM.fit.lp <- nls(DM.red.vec ~ a/(1+(DM.K.vec/b)^c), 
                 data = DM.Kred.dat,
                 algorithm = "port",
                 start = c(a = DM.param.init[1], b = DM.param.init[2], c = DM.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DM.fit.lp.summ <- summary(DM.fit.lp)
plot(DM.K.vec, DM.red.vec, pch=19,xlab="N",ylab="reduction factor")
DM.K.vec.cont <- seq(1,2*DM.pop.found,1)
DM.pred.lp.fx <- coef(DM.fit.lp)[1]/(1+(DM.K.vec.cont/coef(DM.fit.lp)[2])^coef(DM.fit.lp)[3])
lines(DM.K.vec.cont, DM.pred.lp.fx, lty=3,lwd=3,col="red")

DM.a.lp <- coef(DM.fit.lp)[1]
DM.b.lp <- coef(DM.fit.lp)[2]
DM.c.lp <- coef(DM.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
DM.n.mat <- matrix(0, nrow=DM.age.max+1, ncol=(t+1))
DM.n.mat[,1] <- DM.init.vec
DM.popmat <- DM.popmat.orig

## set up projection loop
for (i in 1:t) {
  DM.totN.i <- sum(DM.n.mat[,i])
  DM.pred.red <- as.numeric(DM.a.lp/(1+(DM.totN.i/DM.b.lp)^DM.c.lp))
  diag(DM.popmat[2:(DM.age.max+1),]) <- (DM.Sx[-1])*DM.pred.red
  DM.popmat[DM.age.max+1,DM.age.max+1] <- 0
  DM.popmat[1,] <- DM.m.vec
  DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
}

DM.n.pred <- colSums(DM.n.mat)
plot(yrs, DM.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=DM.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

DM.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DM.s.arr <- DM.m.arr <- array(data=NA, dim=c(t+1, DM.age.max+1, iter))
DM.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  DM.popmat <- DM.popmat.orig
  
  DM.n.mat <- matrix(0, nrow=DM.age.max+1,ncol=(t+1))
  DM.n.mat[,1] <- DM.init.vec
  DM.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    DM.s.alpha <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$alpha
    DM.s.beta <- estBetaParams(DM.Sx, DM.s.sd.vec^2)$beta
    DM.s.stoch <- rbeta(length(DM.s.alpha), DM.s.alpha, DM.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    DM.fert.stch <- rnorm(length(DM.popmat[,1]), DM.m.vec, DM.m.sd.vec)
    DM.m.arr[i,,e] <- ifelse(DM.fert.stch < 0, 0, DM.fert.stch)
    
    DM.totN.i <- sum(DM.n.mat[,i], na.rm=T)
    DM.pred.red.vec[i] <- DM.a.lp/(1+(DM.totN.i/DM.b.lp)^DM.c.lp)
    
    diag(DM.popmat[2:(DM.age.max+1),]) <- (DM.s.stoch[-1])*DM.pred.red.vec[i]
    DM.popmat[DM.age.max+1,DM.age.max+1] <- 0
    DM.popmat[1,] <- DM.m.arr[i,,e]
    DM.n.mat[,i+1] <- DM.popmat %*% DM.n.mat[,i]
    
    DM.s.arr[i,,e] <- DM.s.stoch * DM.pred.red.vec[i]
    
  } # end i loop
  
  DM.n.sums.mat[e,] <- ((as.vector(colSums(DM.n.mat))/DM.pop.found))
  DM.pred.red.mn[e] <- mean(DM.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

DM.n.md <- apply(DM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
DM.n.up <- apply(DM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DM.n.lo <- apply(DM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## TACHYGLOSSUS (aculeatus) (TA)

# mass
TA.mass <- 4.0 # mean(c(3.8,3.4,3.7,(mean(c(3.9,7))))) # short-beaked echidna Tachyglossus aculeatus (Nicol & Andersen 2007)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
TA.rm.allom.pred <- 10^(0.6914 - (0.2622*log10(TA.mass*1000)))
TA.rm.allom.pred

# Metabolic method - Duncan et al. 2007
TA.BMR <- 2.327 # Watts
TA.BMR.mass <- 2909.0 # grams
TA.BMR.wt <- TA.BMR/TA.BMR.mass
BMR.exponent <- 0.725 # Duncan et al. 2007
BMRwt.vec <- c(2.896951E-4,9.03558E-4,0.0027028106,0.011487795,0.027647145,0.04964901)
rmB.vec <- c(0.23230834,0.53004533,1.1605488,3.2882311,6.167897,9.608828)
plot(log10(BMRwt.vec),log10(rmB.vec),pch=19)
fitrmB <- lm(log10(rmB.vec) ~ log10(BMRwt.vec))
abline(fitrmB, lty=2, col="red")
summary(fitrmB)
TA.lrm.pred <- as.numeric(coef(fitrmB)[1] + coef(fitrmB)[2]*log10(TA.BMR.wt))
TA.rm.pred <- 10^TA.lrm.pred

TA.rm.pred <- 0.40
# rm/year = 0.40 (Schmidt-Nielsen K, Dawson T J, Crawford EC (1966) Temperature regulation in the echidna (Tachyglossus aculeatus) J Cell Physiol 67 : 63-72)
TA.lm.pred <- exp(TA.rm.pred)
TA.lm.pred

## for mammal carnivores: intercept: 1.91; slope = −1.02 (upper)
## for mammal carnivores: intercept: -1.17; slope = −0.76 (lower)
TA.D.pred.up <- (10^(1.91 + (-1.02*log10(TA.mass))))/2 # divided by 2 for females only
TA.D.pred.lo <- (10^(-1.17 + (-0.76*log10(TA.mass))))/2 # divided by 2 for females only
# 87/1000 ha = 87/10 km2 = 8.7 km-2 on KI (Rismiller & McKelvey 2000) = 4.35 females/km2
# 168/1000 ha = 167/10 km2 = 16.7 km-2 in Tasmania (Nicol & Anderson 2007) = 8.35 females/km2
TA.D.pred <- TA.D.pred.up # animals/km2 (checks out)
TA.D.pred

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
TA.age.max.allom <- round(10^(0.89 + (0.13*log10(TA.mass*1000))), 0)
TA.age.max.allom # underestimated (torpor, hibernation, low BMR)
TA.age.max <- 45 # (Nicol & Andersen 2007)

## age vector
TA.age.vec <- 0:TA.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
TA.F.allom.pred <- exp(2.719 - (0.211*log(TA.mass*1000)))/2 # divided by 2 for females
TA.F.allom.pred
TA.F.egg <- 1/2 # one egg/year; /2 for daughters
TA.F.pbreed <- 0.55 # females reproductively active (Nicol & Morrow 2012)
# 17 females produced 22 young over 7 years: 22/17/7 = 0.1849, or 0.1849/2 = 0.0924 (Rismiller & McKelvey 2000)
TA.F.pred <- TA.F.egg*TA.F.pbreed

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
TA.alpha.allom <- ceiling(exp(-1.34 + (0.214*log(TA.mass*1000))))
TA.alpha.allom
# only adults breed (Rismiller & McKelvey 2000)
# adult mass reached 3-5 years (Nicol & Andersen 2007)
TA.alpha <- 3

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
TA.s.tran <- ln.a.s + b.s*log(TA.mass*1000) + log(1)
TA.s.ad.yr.allom <- exp(-exp(TA.s.tran))
TA.s.ad.yr.allom
TA.s.ad.yr <- mean(c(0.94,0.98))
TA.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.05*TA.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.9 # rate of mortality decline (also known as bt)
a2 <- 1 - (TA.s.ad.yr) # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.13 # rate of mortality increase
longev <- TA.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
TA.Sx <- c(0.99*TA.s.ad.yr, 1 - qx)
plot(x, TA.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
TA.s.sd.vec <- 0.05*TA.Sx

## pre-breeding design with 0-1 survival in first row
TA.m.vec <- c(0, 0, 0, 0.5*TA.F.pred, 0.75*TA.F.pred, rep(TA.F.pred,(TA.age.max-4))) # 
TA.m.sd.vec <- 0.05*TA.m.vec
plot(TA.age.vec, TA.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
TA.m.dat <- data.frame(TA.age.vec, TA.m.vec)
param.init <- c(0.5, 4, -4)
TA.fit.logp <- nls(TA.m.vec ~ a / (1+(TA.age.vec/b)^c), 
                   data = TA.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TA.fit.logp.summ <- summary(TA.fit.logp)
plot(TA.age.vec, TA.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
TA.age.vec.cont <- seq(0,max(TA.age.vec),1)
TA.pred.p.mm <- coef(TA.fit.logp)[1] / (1+(TA.age.vec.cont/coef(TA.fit.logp)[2])^coef(TA.fit.logp)[3])
#DN.pred.p.mm <- ifelse(DN.pred.p.m > 1, 1, DN.pred.p.m)
lines(TA.age.vec.cont, TA.pred.p.mm,lty=2,lwd=3,col="red")

## create matrix
TA.popmat <- matrix(data = 0, nrow=TA.age.max+1, ncol=TA.age.max+1)
diag(TA.popmat[2:(TA.age.max+1),]) <- TA.Sx[-1]
TA.popmat[TA.age.max+1,TA.age.max+1] <- TA.Sx[TA.age.max+1]
TA.popmat[1,] <- TA.pred.p.mm
colnames(TA.popmat) <- c(0:TA.age.max)
rownames(TA.popmat) <- c(0:TA.age.max)
TA.popmat.orig <- TA.popmat ## save original matrix

## matrix properties
max.lambda(TA.popmat.orig) ## 1-yr lambda
TA.lm.pred
max.r(TA.popmat.orig) # rate of population change, 1-yr
TA.ssd <- stable.stage.dist(TA.popmat.orig) ## stable stage distribution
plot(TA.age.vec, TA.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(TA.popmat.orig, TA.age.max) # reproductive value
TA.gen.l <- G.val(TA.popmat.orig, TA.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
TA.pop.found <- round(area*TA.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
TA.init.vec <- TA.ssd * TA.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*TA.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

TA.tot.F <- sum(TA.popmat.orig[1,])
TA.popmat <- TA.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
TA.n.mat <- matrix(0, nrow=TA.age.max+1,ncol=(t+1))
TA.n.mat[,1] <- TA.init.vec

## set up projection loop
for (i in 1:t) {
  TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
}

TA.n.pred <- colSums(TA.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(TA.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
TA.K.max <- 1*TA.pop.found
TA.K.vec <- c(1, 0.25*TA.K.max, TA.K.max/2, 0.75*TA.K.max, TA.K.max) 
TA.red.vec <- c(1,0.993,0.97,0.93,0.8845)
plot(TA.K.vec, TA.red.vec,pch=19,type="b")
TA.Kred.dat <- data.frame(TA.K.vec, TA.red.vec)

# logistic power function a/(1+(x/b)^c)
TA.param.init <- c(1, TA.K.max, 2)
TA.fit.lp <- nls(TA.red.vec ~ a/(1+(TA.K.vec/b)^c), 
                 data = TA.Kred.dat,
                 algorithm = "port",
                 start = c(a = TA.param.init[1], b = TA.param.init[2], c = TA.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
TA.fit.lp.summ <- summary(TA.fit.lp)
plot(TA.K.vec, TA.red.vec, pch=19,xlab="N",ylab="reduction factor")
TA.K.vec.cont <- seq(1,2*TA.pop.found,1)
TA.pred.lp.fx <- coef(TA.fit.lp)[1]/(1+(TA.K.vec.cont/coef(TA.fit.lp)[2])^coef(TA.fit.lp)[3])
lines(TA.K.vec.cont, TA.pred.lp.fx, lty=3,lwd=3,col="red")

TA.a.lp <- coef(TA.fit.lp)[1]
TA.b.lp <- coef(TA.fit.lp)[2]
TA.c.lp <- coef(TA.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
TA.n.mat <- matrix(0, nrow=TA.age.max+1, ncol=(t+1))
TA.n.mat[,1] <- TA.init.vec
TA.popmat <- TA.popmat.orig

## set up projection loop
for (i in 1:t) {
  TA.totN.i <- sum(TA.n.mat[,i])
  TA.pred.red <- as.numeric(TA.a.lp/(1+(TA.totN.i/TA.b.lp)^TA.c.lp))
  diag(TA.popmat[2:(TA.age.max+1),]) <- (TA.Sx[-1])*TA.pred.red
  TA.popmat[TA.age.max+1,TA.age.max+1] <- 0
  TA.popmat[1,] <- TA.m.vec
  TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
}

TA.n.pred <- colSums(TA.n.mat)
plot(yrs, TA.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=TA.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

TA.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
TA.s.arr <- TA.m.arr <- array(data=NA, dim=c(t+1, TA.age.max+1, iter))
TA.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  TA.popmat <- TA.popmat.orig
  
  TA.n.mat <- matrix(0, nrow=TA.age.max+1,ncol=(t+1))
  TA.n.mat[,1] <- TA.init.vec
  TA.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    TA.s.alpha <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$alpha
    TA.s.beta <- estBetaParams(TA.Sx, TA.s.sd.vec^2)$beta
    TA.s.stoch <- rbeta(length(TA.s.alpha), TA.s.alpha, TA.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    TA.fert.stch <- rnorm(length(TA.popmat[,1]), TA.m.vec, TA.m.sd.vec)
    TA.m.arr[i,,e] <- ifelse(TA.fert.stch < 0, 0, TA.fert.stch)
    
    TA.totN.i <- sum(TA.n.mat[,i], na.rm=T)
    TA.pred.red.vec[i] <- TA.a.lp/(1+(TA.totN.i/TA.b.lp)^TA.c.lp)
    
    diag(TA.popmat[2:(TA.age.max+1),]) <- (TA.s.stoch[-1])*TA.pred.red.vec[i]
    TA.popmat[TA.age.max+1,TA.age.max+1] <- 0
    TA.popmat[1,] <- TA.m.arr[i,,e]
    TA.n.mat[,i+1] <- TA.popmat %*% TA.n.mat[,i]
    
    TA.s.arr[i,,e] <- TA.s.stoch * TA.pred.red.vec[i]
    
  } # end i loop
  
  TA.n.sums.mat[e,] <- ((as.vector(colSums(TA.n.mat))/TA.pop.found))
  TA.pred.red.mn[e] <- mean(TA.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

TA.n.md <- apply(TA.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
TA.n.up <- apply(TA.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
TA.n.lo <- apply(TA.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations


##############################
## MEGALIBGWILIA (ramsayi) (MR)

# mass
MR.mass <- 11 # 10-12 kg Megalibgwilia ramsayi (Johnson 2006)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
MR.rm.allom.pred <- 10^(0.6914 - (0.2622*log10(MR.mass*1000)))
MR.rm.allom.pred
MR.rm.pred <- (TA.rm.pred/TA.rm.allom.pred) * MR.rm.allom.pred # correct for over-estimation based on TA data
MR.rm.pred
MR.lm.pred <- exp(MR.rm.pred)
MR.lm.pred

## for mammal carnivores: intercept: 1.91; slope = −1.02 (upper)
## for mammal carnivores: intercept: -1.17; slope = −0.76 (lower)
MR.D.pred.up <- (10^(1.91 + (-1.02*log10(MR.mass))))/2 # divided by 2 for females only
MR.D.pred.lo <- (10^(-1.17 + (-0.76*log10(MR.mass))))/2 # divided by 2 for females only
# 87/1000 ha = 87/10 km2 = 8.7 km-2 on KI (Rismiller & McKelvey 2000) = 4.35 females/km2
# 168/1000 ha = 167/10 km2 = 16.7 km-2 in Tasmania (Nicol & Anderson 2007) = 8.35 females/km2
MR.D.pred <- MR.D.pred.up # animals/km2 (checks out)
MR.D.pred

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
MR.age.max.allom <- round(10^(0.89 + (0.13*log10(MR.mass*1000))), 0)
MR.age.max.allom # underestimated (torpor, hibernation, low BMR)
MR.age.max <- round((TA.age.max/TA.age.max.allom) * MR.age.max.allom, 0) # correct for over-estimation based on TA data
MR.age.max

## age vector
MR.age.vec <- 0:MR.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
MR.F.allom.pred <- exp(2.719 - (0.211*log(MR.mass*1000)))/2 # divided by 2 for females
MR.F.allom.pred
MR.F.pred <- (TA.F.pred/(TA.F.allom.pred/2)) * (MR.F.allom.pred/2) # correct for over-estimation based on TA data
MR.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
MR.alpha.allom <- ceiling(exp(-1.34 + (0.214*log(MR.mass*1000))))
MR.alpha.allom
MR.alpha <- round((TA.alpha/TA.alpha.allom) * MR.alpha.allom, 0) # correct for over-estimation based on TA data
MR.alpha

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
MR.s.tran <- ln.a.s + b.s*log(MR.mass*1000) + log(1)
MR.s.ad.yr.allom <- exp(-exp(MR.s.tran))
MR.s.ad.yr.allom
MR.s.ad.yr <- (TA.s.ad.yr/TA.s.ad.yr.allom) * MR.s.ad.yr.allom # correct for over-estimation based on TA data
MR.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.01*MR.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.5 # rate of mortality decline (also known as bt)
a2 <- 1 - (MR.s.ad.yr) # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.1 # rate of mortality increase
longev <- MR.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
MR.Sx <- c(0.985*MR.s.ad.yr, 1 - qx)
plot(x, MR.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
MR.s.sd.vec <- 0.05*MR.Sx

## pre-breeding design with 0-1 survival in first row
MR.m.vec <- c(0, 0, 0, 0.5*MR.F.pred, 0.75*MR.F.pred, rep(MR.F.pred,(MR.age.max-4))) # 
MR.m.sd.vec <- 0.05*MR.m.vec
plot(MR.age.vec, MR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
MR.m.dat <- data.frame(MR.age.vec, MR.m.vec)
param.init <- c(0.5, 4, -4)
MR.fit.logp <- nls(MR.m.vec ~ a / (1+(MR.age.vec/b)^c), 
                   data = MR.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
MR.fit.logp.summ <- summary(MR.fit.logp)
plot(MR.age.vec, MR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
MR.age.vec.cont <- seq(0,max(MR.age.vec),1)
MR.pred.p.mm <- coef(MR.fit.logp)[1] / (1+(MR.age.vec.cont/coef(MR.fit.logp)[2])^coef(MR.fit.logp)[3])
#DN.pred.p.mm <- ifelse(DN.pred.p.m > 1, 1, DN.pred.p.m)
lines(MR.age.vec.cont, MR.pred.p.mm,lty=2,lwd=3,col="red")

## create matrix
MR.popmat <- matrix(data = 0, nrow=MR.age.max+1, ncol=MR.age.max+1)
diag(MR.popmat[2:(MR.age.max+1),]) <- MR.Sx[-1]
MR.popmat[MR.age.max+1,MR.age.max+1] <- MR.Sx[MR.age.max+1]
MR.popmat[1,] <- MR.pred.p.mm
colnames(MR.popmat) <- c(0:MR.age.max)
rownames(MR.popmat) <- c(0:MR.age.max)
MR.popmat.orig <- MR.popmat ## save original matrix

## matrix properties
max.lambda(MR.popmat.orig) ## 1-yr lambda
MR.lm.pred
max.r(MR.popmat.orig) # rate of population change, 1-yr
MR.ssd <- stable.stage.dist(MR.popmat.orig) ## stable stage distribution
plot(MR.age.vec, MR.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(MR.popmat.orig, MR.age.max) # reproductive value
MR.gen.l <- G.val(MR.popmat.orig, MR.age.max) # mean generation length
MR.gen.l

## initial population vector
area <- 500*500 # km × km
MR.pop.found <- round(area*MR.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
MR.init.vec <- MR.ssd * MR.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*MR.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

MR.tot.F <- sum(MR.popmat.orig[1,])
MR.popmat <- MR.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
MR.n.mat <- matrix(0, nrow=MR.age.max+1,ncol=(t+1))
MR.n.mat[,1] <- MR.init.vec

## set up projection loop
for (i in 1:t) {
  MR.n.mat[,i+1] <- MR.popmat %*% MR.n.mat[,i]
}

MR.n.pred <- colSums(MR.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(MR.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
MR.K.max <- 1*MR.pop.found
MR.K.vec <- c(1, 0.25*MR.K.max, MR.K.max/2, 0.75*MR.K.max, MR.K.max) 
MR.red.vec <- c(1,0.991,0.968,0.934,0.8897)
plot(MR.K.vec, MR.red.vec,pch=19,type="b")
MR.Kred.dat <- data.frame(MR.K.vec, MR.red.vec)

# logistic power function a/(1+(x/b)^c)
MR.param.init <- c(1, MR.K.max, 2)
MR.fit.lp <- nls(MR.red.vec ~ a/(1+(MR.K.vec/b)^c), 
                 data = MR.Kred.dat,
                 algorithm = "port",
                 start = c(a = MR.param.init[1], b = MR.param.init[2], c = MR.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
MR.fit.lp.summ <- summary(MR.fit.lp)
plot(MR.K.vec, MR.red.vec, pch=19,xlab="N",ylab="reduction factor")
MR.K.vec.cont <- seq(1,2*MR.pop.found,1)
MR.pred.lp.fx <- coef(MR.fit.lp)[1]/(1+(MR.K.vec.cont/coef(MR.fit.lp)[2])^coef(MR.fit.lp)[3])
lines(MR.K.vec.cont, MR.pred.lp.fx, lty=3,lwd=3,col="red")

MR.a.lp <- coef(MR.fit.lp)[1]
MR.b.lp <- coef(MR.fit.lp)[2]
MR.c.lp <- coef(MR.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
MR.n.mat <- matrix(0, nrow=MR.age.max+1, ncol=(t+1))
MR.n.mat[,1] <- MR.init.vec
MR.popmat <- MR.popmat.orig

## set up projection loop
for (i in 1:t) {
  MR.totN.i <- sum(MR.n.mat[,i])
  MR.pred.red <- as.numeric(MR.a.lp/(1+(MR.totN.i/MR.b.lp)^MR.c.lp))
  diag(MR.popmat[2:(MR.age.max+1),]) <- (MR.Sx[-1])*MR.pred.red
  MR.popmat[MR.age.max+1,MR.age.max+1] <- 0
  MR.popmat[1,] <- MR.m.vec
  MR.n.mat[,i+1] <- MR.popmat %*% MR.n.mat[,i]
}

MR.n.pred <- colSums(MR.n.mat)
plot(yrs, MR.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=MR.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 10000
itdiv <- iter/10

MR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
MR.s.arr <- MR.m.arr <- array(data=NA, dim=c(t+1, MR.age.max+1, iter))
MR.pred.red.mn <- rep(0,iter)

for (e in 1:iter) {
  MR.popmat <- MR.popmat.orig
  
  MR.n.mat <- matrix(0, nrow=MR.age.max+1,ncol=(t+1))
  MR.n.mat[,1] <- MR.init.vec
  MR.pred.red.vec <- rep(0,t)
  
  for (i in 1:t) {
    # stochastic survival values
    MR.s.alpha <- estBetaParams(MR.Sx, MR.s.sd.vec^2)$alpha
    MR.s.beta <- estBetaParams(MR.Sx, MR.s.sd.vec^2)$beta
    MR.s.stoch <- rbeta(length(MR.s.alpha), MR.s.alpha, MR.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    MR.fert.stch <- rnorm(length(MR.popmat[,1]), MR.m.vec, MR.m.sd.vec)
    MR.m.arr[i,,e] <- ifelse(MR.fert.stch < 0, 0, MR.fert.stch)
    
    MR.totN.i <- sum(MR.n.mat[,i], na.rm=T)
    MR.pred.red.vec[i] <- MR.a.lp/(1+(MR.totN.i/MR.b.lp)^MR.c.lp)
    
    diag(MR.popmat[2:(MR.age.max+1),]) <- (MR.s.stoch[-1])*MR.pred.red.vec[i]
    MR.popmat[MR.age.max+1,MR.age.max+1] <- 0
    MR.popmat[1,] <- MR.m.arr[i,,e]
    MR.n.mat[,i+1] <- MR.popmat %*% MR.n.mat[,i]
    
    MR.s.arr[i,,e] <- MR.s.stoch * MR.pred.red.vec[i]
    
  } # end i loop
  
  MR.n.sums.mat[e,] <- ((as.vector(colSums(MR.n.mat))/MR.pop.found))
  MR.pred.red.mn[e] <- mean(MR.pred.red.vec, na.rm=T)
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

MR.n.md <- apply(MR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
MR.n.up <- apply(MR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
MR.n.lo <- apply(MR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
