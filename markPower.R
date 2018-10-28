source("markDensSim.R")
library(np)

# Calculation of power for the 30 scenarios: 
# (AMP-B, AMP-C, pooled) x (IC50, IC80) x (0, 0.3, 0.5, 0.7, 0.9)

dataDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP"

#========== Running simulations ================================================
# universal parameters
alphaLR <- 0.1
alphaWald <- 0.05
simulations <- 50  # 10000
power <- as.data.frame(matrix(data=0,nrow=30, ncol=4))
rownames(power) <- c("B_0_ic50", "B_0.3_ic50", "B_0.5_ic50", "B_0.7_ic50", "B_0.9_ic50", 
                     "C_0_ic50", "C_0.3_ic50", "C_0.5_ic50", "C_0.7_ic50", "C_0.9_ic50",
                     "Pooled_0_ic50", "Pooled_0.3_ic50", "Pooled_0.5_ic50", "Pooled_0.7_ic50", "Pooled_0.9_ic50",
                     "B_0_ic80", "B_0.3_ic80", "B_0.5_ic80", "B_0.7_ic80", "B_0.9_ic80", 
                     "C_0_ic80", "C_0.3_ic80", "C_0.5_ic80", "C_0.7_ic80", "C_0.9_ic80",
                     "Pooled_0_ic80", "Pooled_0.3_ic80", "Pooled_0.5_ic80", "Pooled_0.7_ic80", "Pooled_0.9_ic80")
# twosidedLR is for comparing the 2-sided LR p-value with the 1-sided LR p-value
colnames(power) <- c("WaldH00", "WaldH0", "LR", "twosidedLR")  

#======================= AMP-B, IC50 ========================================

Np = 900     # number of subjects in placebo group
np = 34      # number of cases in placebo group
taumax = 80  # follow-up time in weeks
markVE <- 0  # cutoff mark VE value (i.e., VE value corresponding to a mark of 0.3)
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))
log10ic50 <- data$ic50.geometric.mean.imputed.log10

# nonparametric density estimation using the package 'np'
densbw <- npudensbw(~ log10ic50, ckertype="epanechnikov")  # bandwidth selection
dens <- npudens(densbw)

# power calculations saved in a repeatedly updated dataframe, 'power'
power <- calcPower(1, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.3
power <- calcPower(2, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.5
power <- calcPower(3, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.7
power <- calcPower(4, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.9
power <- calcPower(5, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)

#==================== AMP-C, IC50 =============================================
Np = 634
np = 34
taumax = 80
markVE <- 0
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_c.csv"))
log10ic50 <- data$ic50.geometric.mean.imputed.log10

densbw <- npudensbw(~ log10ic50, ckertype="epanechnikov") 
dens <- npudens(densbw)

power <- calcPower(6, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.3
power <- calcPower(7, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.5
power <- calcPower(8, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.7
power <- calcPower(9, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.9
power <- calcPower(10, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)

#===================== Pooled, IC50 ================================= 
Np = 1534
np = 68
taumax = 80 
markVE <- 0
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_all.csv"))
log10ic50 <- data$ic50.geometric.mean.imputed.log10

densbw <- npudensbw(~ log10ic50, ckertype="epanechnikov")  
dens <- npudens(densbw)

power <- calcPower(11, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.3
power <- calcPower(12, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.5
power <- calcPower(13, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.7
power <- calcPower(14, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)
markVE <- 0.9
power <- calcPower(15, simulations, Np, np, markVE, taumax, dens, "log10ic50", alphaLR, alphaWald, power)

#======================= AMP-B, IC80 =================================================

Np = 900
np = 34
taumax = 80 
markVE <- 0
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))
log10ic80 <- data$ic80.geometric.mean.imputed.log10
log10ic80 <- na.omit(log10ic80)

densbw <- npudensbw(~ log10ic80, ckertype="epanechnikov")  
dens <- npudens(densbw)

power <- calcPower(16, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.3
power <- calcPower(17, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.5
power <- calcPower(18, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.7
power <- calcPower(19, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.9
power <- calcPower(20, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)

#==================== AMP-C, CI80 ====================================================
Np = 634
np = 34
taumax = 80  
markVE <- 0
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_c.csv"))
log10ic80 <- data$ic80.geometric.mean.imputed.log10
log10ic80 <- na.omit(log10ic80)

densbw <- npudensbw(~ log10ic80, ckertype="epanechnikov")  
dens <- npudens(densbw)

power <- calcPower(21, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.3
power <- calcPower(22, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.5
power <- calcPower(23, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.7
power <- calcPower(24, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.9
power <- calcPower(25, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)

#===================== Pooled, IC80 =================================================== 
Np = 1534
np = 68
taumax = 80  
markVE <- 0
data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_all.csv"))
log10ic80 <- data$ic80.geometric.mean.imputed.log10
log10ic80 <- na.omit(log10ic80)

densbw <- npudensbw(~ log10ic80, ckertype="epanechnikov") 
dens <- npudens(densbw)

power <- calcPower(26, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.3
power <- calcPower(27, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.5
power <- calcPower(28, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.7
power <- calcPower(29, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)
markVE <- 0.9
power <- calcPower(30, simulations, Np, np, markVE, taumax, dens, "log10ic80", alphaLR, alphaWald, power)

#===================== Save power calculations ===================================
power <- power/simulations
write.csv(power, "power.csv")


# # the below powers should be comparable
# mean((lrBeta.pvals <= 0.1) & (thetaHat[2] > 0))  # p-val is below 2-sided alpha and point estimate is in right direction
# mean(waldH00.pvals <= 0.05)  # p-val is below 1-sided alpha



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
