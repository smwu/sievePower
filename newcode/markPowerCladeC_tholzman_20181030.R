library("rslurm")

library(np)
library(survival)

source("markDensSim_tholzman_20181029.R")
source("../markHR.R")
source("../covEst.R")

# load data
#dataDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP"
dataDir <- "../data"

mascola <- read.csv(file.path(dataDir,"Mascola_Acute_Clade_C_VRC01.csv"), header=TRUE, stringsAsFactors = FALSE)
mascola$IC50[mascola$IC50==">10"] <- 20
mascola$IC80[mascola$IC80==">10"] <- 20

#colnames(mascola) <-
#    c("ic50.geometric.mean.imputed.log10", "ic80.geometric.mean.imputed.log10")

#new line from Stephanie
colnames(mascola) <- c("ID", "ic50.geometric.mean.imputed.log10", "ic80.geometric.mean.imputed.log10")
mascola$ic50.geometric.mean.imputed.log10 <- as.numeric(mascola$ic50.geometric.mean.imputed.log10)
mascola$ic80.geometric.mean.imputed.log10 <- as.numeric(mascola$ic80.geometric.mean.imputed.log10)
#nMC <- 1

nMC <- 1000

dataC <- mascola
dataB <- mascola
dataBandC <- mascola

pars <- data.frame(
    seed=seq(1,nMC)
    )
pars$trial <- rep("703",times=nMC)
pars$randRatio <- rep(I(list(c("2:1", "1:1"))),times=nMC)
pars$IC <- rep(I(list(c("IC50", "IC80"))),times=nMC)
pars$VEcoord <- rep(I(list(c(0, 0.3, 0.5, 0.7, 0.9))),times=nMC)
pars$taumax <- rep(80,times=nMC)
pars$alpha1sided <- rep(0.05,times=nMC)

results <- slurm_apply(
    getInferenceOneMC,
    pars,
    add_objects=c("dataB","dataC","dataBandC"),
    jobname="markPowerCladeC",
    nodes=1000,
    cpus_per_node=1,
    submit=FALSE
    )


#results <- lapply(1:nMC, getInferenceOneMC, 
#                  trial="703", 
#                  randRatio=c("2:1", "1:1"), 
#                  IC=c("IC50", "IC80"),
#                  VEcoord=c(0, 0.3, 0.5, 0.7, 0.9),
#                  taumax=80,
#                  alpha1sided=0.05,
#                  dataC=mascola)
print(results)
save(results, file="resultsCladeC.RData")


# #===============================================================
# # revised 'calcPower' function
# calcPower <- function(index, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, alphaLR, alphaWald, randomRatio) {
#   WaldH00 <- 0
#   WaldH0 <- 0
#   LR <- 0
#   twosidedLR <- 0
#   for (i in 1:simulations) {
#     result <- simulOne(Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, randomRatio, 1)
#     
#     # one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
#     if (result$H00waldP1sided <= alphaWald) {
#       WaldH00 <- WaldH00 + 1
#     }
#     # 1-sided Wald-type test of H0: VE(v)=VE vs. alternative that beta > 0
#     if (result$H0waldP1sided <= alphaWald) {
#       WaldH0 <- WaldH0 + 1
#     }
#     # 1-sided likelihood ratio test (serves as sanity check)
#     if ((result$H0lrP2sided <= alphaLR) & (result$betaHat > 0)) { 
#       LR <- LR + 1
#     }
#     # two-sided likelihood ratio p-value
#     if (result$H0lrP2sided <= alphaLR) {
#       twosidedLR <- twosidedLR + 1
#     }
#   }
#   return(list(WaldH00 = WaldH00/simulations, WaldH0 = WaldH0/simulations, LR = LR/simulations, twosidedLR = twosidedLR/simulations))
# }      
# 
# 
# # set parameters
# alphaLR <- 0.1
# alphaWald <- 0.05
# simulations <- 10000
# markVE <- c(0, 0.3, 0.5, 0.7, 0.9)
# 
# Np = 634
# np = 34
# taumax = 80
# 
# # rate parameters for failure time T and censoring time C
# lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
# lambdaC <- 0.1*Np/np*lambdaT
# 
# ### Analyses for ic50
# mark <- log10(as.numeric(mascola$IC50))
# densbw <- npudensbw(~ mark, ckertype="epanechnikov")  # bandwidth selection
# dens <- npudens(densbw)
# 
# # initialize dataframe for power calculations
# ptm <- proc.time()
# powerIC50 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
# powerIC50_1to1 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
# colnames(powerIC50) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")  
# colnames(powerIC50_1to1) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")
# 
# 
# for (i in 1:length(markVE)) {
#   coeff <- getAlphaBetaGamma(markVE[i], dens, "mark")
#   alpha <- coeff$alpha
#   beta <- coeff$beta
#   gamma <- coeff$gamma
#   
#   calc <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "mark", alphaLR, alphaWald, 2)
#   powerIC50$WaldH00[i] <- calc$WaldH00
#   powerIC50$WaldH0[i] <- calc$WaldH0
#   powerIC50$LR[i] <- calc$LR
#   powerIC50$twosidedLR[i] <- calc$twosidedLR
#   
#   calc1to1 <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "mark", alphaLR, alphaWald, 1)
#   powerIC50_1to1$WaldH00[i] <- calc1to1$WaldH00
#   powerIC50_1to1$WaldH0[i] <- calc1to1$WaldH0
#   powerIC50_1to1$LR[i] <- calc1to1$LR
#   powerIC50_1to1$twosidedLR[i] <- calc1to1$twosidedLR
# }
# rownames(powerIC50) <- paste0("703CladeC_2:1_IC50_",markVE)
# rownames(powerIC50_1to1) <- paste0("703CladeC_1:1_IC50_",markVE)
# 
# 
# ### Analyses for ic80
# mark <- log10(as.numeric(mascola$IC80))
# densbw <- npudensbw(~ mark, ckertype="epanechnikov")  # bandwidth selection
# dens <- npudens(densbw)
# 
# powerIC80 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
# powerIC80_1to1 <- as.data.frame(matrix(0, nrow=length(markVE), ncol=4))
# colnames(powerIC80) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")  
# colnames(powerIC80_1to1) <- c("WaldH00", "WaldH0", "LR", "twosidedLR") 
# 
# for (i in 1:length(markVE)) {
#   coeff <- getAlphaBetaGamma(markVE[i], dens, "mark")
#   alpha <- coeff$alpha
#   beta <- coeff$beta
#   gamma <- coeff$gamma
#   
#   calc <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "mark", alphaLR, alphaWald, 2)
#   powerIC80$WaldH00[i] <- calc$WaldH00
#   powerIC80$WaldH0[i] <- calc$WaldH0
#   powerIC80$LR[i] <- calc$LR
#   powerIC80$twosidedLR[i] <- calc$twosidedLR
#   
#   calc1to1 <- calcPower(i, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, "mark", alphaLR, alphaWald, 1)
#   powerIC80_1to1$WaldH00[i] <- calc1to1$WaldH00
#   powerIC80_1to1$WaldH0[i] <- calc1to1$WaldH0
#   powerIC80_1to1$LR[i] <- calc1to1$LR
#   powerIC80_1to1$twosidedLR[i] <- calc1to1$twosidedLR
# }
# rownames(powerIC80) <- paste0("703CladeC_2:1_IC80_",markVE)
# rownames(powerIC80_1to1) <- paste0("703CladeC_1:1_IC80_",markVE)
# 
# # create and save combined dataframe 
# power <- rbind(powerIC50, powerIC80, powerIC50_1to1, powerIC80_1to1)
# write.csv(power, "power_CladeC.csv")
# proc.time()-ptm
# 
# 
# 
#   
