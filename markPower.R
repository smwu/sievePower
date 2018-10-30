# Calculation of power for 60 scenarios: 
# (AMP-B, AMP-C, pooled) x (2:1, 1:1) x (IC50, IC80) x (0, 0.3, 0.5, 0.7, 0.9)

rm(list=ls(all=TRUE))

dataDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP"
# Ted, you obviously need to change 'codeDir' to wherever you have your local repo
codeDir <- "h:/SCHARP/HVTN703/sievePower"
outDir <- "h:/SCHARP/HVTN703/sievePower/Routput"

library(survival)
library(np)
source(file.path(codeDir, "markHR.R"))
source(file.path(codeDir, "covEst.R"))
source(file.path(codeDir, "markDensSim.R"))

dataB <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))
dataC <- read.csv(file.path(dataDir, "catnap_vrc01_neut_c.csv"))
dataBandC <- read.csv(file.path(dataDir, "catnap_vrc01_neut_all.csv"))

nMC <- 1000
results <- lapply(1:nMC, getInferenceOneMC, 
                  trial=c("704", "703", "704and703"), 
                  randRatio=c("2:1", "1:1"), 
                  IC=c("IC50", "IC80"),
                  VEcoord=c(0, 0.3, 0.5, 0.7, 0.9),
                  taumax=80,
                  alpha1sided=0.05,
                  dataB=dataB,
                  dataC=dataC,
                  dataBandC=dataBandC)

save(results, file=file.path(outDir, "results.RData"))

# load the output list named 'res' from the parallelized version of the above 'lapply'
load(file.path(outDir, "results.RData"))

# calculate power estimates from 'results' and store them in data frame 'power'
trial <- c("704", "703", "704and703")
randRatio <- c("2:1", "1:1")
IC <- c("IC50", "IC80")
VEcoord <- c(0, 0.3, 0.5, 0.7, 0.9)
# these are the scenario labels in 'results'
scenario <- c(outer(c(outer(c(outer(trial, randRatio, paste1)), IC, paste1)), VEcoord, paste1))

power <- as.data.frame(matrix(NA, nrow=length(scenario), ncol=3))
# use the same labels in the data frame with power estimates
rownames(power) <- scenario
colnames(power) <- c("H00wald1sidedPower", "H0wald1sidedPower", "H0lr1sidedPower")

# these are the p-value labels in 'res' corresponding to colnames(power)
pvalLabels <- c("H00waldP1sided", "H0waldP1sided", "H0lrP2sided")
for (i in 1:NROW(power)){
  for (j in 1:NCOL(power)){
    pvals <- sapply(res, function(resultsOneMC, scenarioLabel, pvalLabel){ 
      p <- resultsOneMC[[scenarioLabel]][[pvalLabel]]
      return(ifelse(is.null(p), NA, p))
      }, scenarioLabel=scenario[i], pvalLabel=pvalLabels[j])
    
    if (j<3){  # for the Wald tests
      power[i,j] <- mean(pvals<=0.05, na.rm=TRUE)
    } else {  # for the LR test
      # we also need to extract the point estimates of beta
      betaHats <- sapply(res, function(resultsOneMC, scenarioLabel){ 
        b <- resultsOneMC[[scenarioLabel]][["betaHat"]]
        return(ifelse(is.null(b), NA, b))
        }, scenarioLabel=scenario[i])
      power[i,j] <- mean((pvals<=0.1) & (betaHats>0), na.rm=TRUE)
    }
  }
}

write.table(power, file=file.path(outDir, "power.csv"), row.names=TRUE, col.names=TRUE, quote=FALSE)



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
