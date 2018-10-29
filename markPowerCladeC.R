source("markDensSim.R")
library(np)

# revised 'calcPower' function
calcPower <- function(index, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, alphaLR, alphaWald, randomRatio) {
  WaldH00 <- 0
  WaldH0 <- 0
  LR <- 0
  twosidedLR <- 0
  for (i in 1:simulations) {
    result <- simulOne(Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, randomRatio)
    
    # one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
    if (result$weighted.waldH00.pval <= alphaWald) {
      WaldH00 <- WaldH00 + 1
    }
    # 1-sided Wald-type test of H0: VE(v)=VE vs. alternative that beta > 0
    if (result$waldH0.pval <= alphaWald) {
      WaldH0 <- WaldH0 + 1
    }
    # 1-sided likelihood ratio test (serves as sanity check)
    if ((result$lrBeta.pval <= alphaLR) & (result$thetaHat[2] > 0)) { 
      LR <- LR + 1
    }
    # two-sided likelihood ratio p-value
    if (result$lrBeta.pval <= alphaLR) {
      twosidedLR <- twosidedLR + 1
    }
  }
  return(list(WaldH00 = WaldH00/simulations, WaldH0 = WaldH0/simulations, LR = LR/simulations, twosidedLR = twosidedLR/simulations))
}


# load data
dataDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP"

mascola <- read.csv(file.path(dataDir,"Mascola_Acute_Clade_C_VRC01.csv"), header=TRUE, stringsAsFactors = FALSE)
mascola$IC50[mascola$IC50==">10"] <- "20"
mascola$IC80[mascola$IC80==">10"] <- "20"

# set parameters
alphaLR <- 0.1
alphaWald <- 0.05
simulations <- 20  # 10000
markVE <- c(0, 0.3, 0.5, 0.7, 0.9)

Np = 634
np = 34
taumax = 80

# rate parameters for failure time T and censoring time C
lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
lambdaC <- 0.1*Np/np*lambdaT


### Analyses for ic50
log10ic50 <- log10(as.numeric(mascola$IC50))
densbw <- npudensbw(~ log10ic50, ckertype="epanechnikov")  # bandwidth selection
dens <- npudens(densbw)

# initialize dataframe for power calculations
powerIC50 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
powerIC50_1to1 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
colnames(powerIC50) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")  
colnames(powerIC50_1to1) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")  

for (i in 1:length(markVE)) {
  coeff <- getAlphaBetaGamma(markVE[i], dens, "log10ic50")
  alpha <- coeff$alpha
  beta <- coeff$beta
  gamma <- coeff$gamma
  
  calc <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "log10ic50", alphaLR, alphaWald, 2)
  powerIC50$WaldH00[i] <- calc$WaldH00
  powerIC50$WaldH0[i] <- calc$WaldH0
  powerIC50$LR[i] <- calc$LR
  powerIC50$twosidedLR[i] <- calc$twosidedLR
  
  calc1to1 <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "log10ic50", alphaLR, alphaWald, 1)
  powerIC50_1to1$WaldH00[i] <- calc1to1$WaldH00
  powerIC50_1to1$WaldH0[i] <- calc1to1$WaldH0
  powerIC50_1to1$LR[i] <- calc1to1$LR
  powerIC50_1to1$twosidedLR[i] <- calc1to1$twosidedLR
}
rownames(powerIC50) <- paste0("CladeC_",markVE,"_ic50_2to1")
rownames(powerIC50_1to1) <- paste0("CladeC_",markVE,"_ic50_1to1")


### Analyses for ic80
log10ic80 <- log10(as.numeric(mascola$IC80))
densbw <- npudensbw(~ log10ic80, ckertype="epanechnikov")  # bandwidth selection
dens <- npudens(densbw)

powerIC80 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
powerIC80_1to1 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
colnames(powerIC80) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")  
colnames(powerIC80_1to1) <- c("WaldH00", "WaldH0", "LR", "twosidedLR") 

for (i in 1:length(markVE)) {
  coeff <- getAlphaBetaGamma(markVE[i], dens, "log10ic80")
  alpha <- coeff$alpha
  beta <- coeff$beta
  gamma <- coeff$gamma
  
  calc <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "log10ic80", alphaLR, alphaWald, 2)
  powerIC80$WaldH00[i] <- calc$WaldH00
  powerIC80$WaldH0[i] <- calc$WaldH0
  powerIC80$LR[i] <- calc$LR
  powerIC80$twosidedLR[i] <- calc$twosidedLR
  
  calc1to1 <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "log10ic80", alphaLR, alphaWald, 1)
  powerIC80_1to1$WaldH00[i] <- calc1to1$WaldH00
  powerIC80_1to1$WaldH0[i] <- calc1to1$WaldH0
  powerIC80_1to1$LR[i] <- calc1to1$LR
  powerIC80_1to1$twosidedLR[i] <- calc1to1$twosidedLR
}
rownames(powerIC80) <- paste0("CladeC_",markVE,"_ic80_2to1")
rownames(powerIC80_1to1) <- paste0("CladeC_",markVE,"_ic80_1to1")

# create and save combined dataframe 
power <- rbind(powerIC50, powerIC80, powerIC50_1to1, powerIC80_1to1)
write.csv(power, "power_CladeC.csv")



  