source("markDensSim.R")
library(np)

# Calculation of power for the 24 scenarios: 
# (AMP-B, AMP-C, pooled) x (IC50, IC80) x (0.45, 0.6, 0.75, 0.9)

dataDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP"

#========== Running simulations ================================================
# universal parameters
alphaLR <- 0.1
alphaWald <- 0.05
powerLR <- numeric(24)
powerH0 <- numeric(24)
powerH00 <- numeric(24)
simulations <- 10  # 10000

# for comparing 2-sided LR p-value with 1-sided LR p-value
power2sidedLR <- numeric(24)

#=========
# AMP-B:

Np = 900
np = 34
taumax = 80  # weeks
markVE <- 0.45
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))

ic50 <- data$ic50.geometric.mean.imputed
log10ic50 <- data$ic50.geometric.mean.imputed.log10

# nonparametric density estimation using the package 'np'
densbw <- npudensbw(~ log10ic50, ckertype="epanechnikov")  # bandwidth selection
dens <- npudens(densbw)

for (i in 1:simulations) {
  result <- simulOne(Np, np, markVE, taumax, dens, "log10ic50")
  
  # likelihood ratio test (serves as sanity check)
  if (result$lrBeta.pval <= alphaLR & (result$beta > 0)) { 
    powerLR[1] <- powerLR[1] + 1
  }
  # two-sided likelihood ratio p-value
  if (result$lrBeta.pval <= alphaLR) { 
    power2sidedLR[1] <- power2sidedLR[1] + 1
  }
  
  # 1-sided test of H0: VE(v)=VE vs. alternative that beta > 0
  if (result$waldH0.pval <= alphaWald) {
    powerH0[1] <- powerH0[1] + 1
  }
  
  # one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
  if (result$weighted.waldH00.pval) {
    powerH00[1] <- powerH0[1] + 1
  }

}

# # the below powers should be comparable
# mean((lrBeta.pvals <= 0.1) & (thetaHat[2] > 0))  # p-val is below 2-sided alpha and point estimate is in right direction
# mean(waldH00.pvals <= 0.05)  # p-val is below 1-sided alpha

powerLR[1] <- powerLR[1]/simulations
power2sidedLR[1] <- power2sidedLR[1]/simulations
powerH0[1] <- powerH0[1]/simulations
powerH00[1] <- powerH00[1]/simulations

#========================================== Miscellaneous==========================
result <- simulOne(Np, np, markVE, taumax, dens, "log10ic50")
##sanity check tat values of covEst are small
result$covEst

# this illustrates how the nonparametric density can be integrated over
# also a sanity check that this is a true density, i.e., the integral = 1
integrate(dPredict, lower=-Inf, upper=Inf, npdensityObject=dens, varName="log10ic50")

hist(log10ic50, breaks=13, freq=FALSE)  # histogram of density of log10ic50
x <- seq(-3.5, 3.5, by=0.1)
plot(x, dPredict(x, dens, "log10ic50"), col="red", type="l")  # predicted density 
lines(x, dPredict(x, dens, "log10ic50")*exp(alpha+beta*x), col="blue")

# creating the density for the mark variable in placebos
test <- seq(-100, 100, length=100)
dPredict(test, dens, "log10ic50")
b <- function(x, npdensityObject, varName, alpha, beta) {
  dPredict(x, npdensityObject, varName)*exp(alpha + beta*x)
}
b(test, dens, "log10ic50", alpha=-10, beta=beta)
# sanity check that this is a true density, i.e., the integral = 1
integrate(b, lower=-100, upper=100, npdensityObject=dens, varName="log10ic50", alpha=alpha, beta=beta, subdivisions=2000)$value

# # Maximum likelihood estimation for Exp(lambda) gives lambdaHat = 1/mean(data)
# lambdaVraw <- 1/mean(IC50subB$ic50.geometric.mean.imputed)
