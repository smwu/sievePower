library(np)

# 'dPredict' returns nonparametric density estimates at 'x'
# 'npdensityObject' is the output object from 'npudens'
# 'varName' is the variable name used in the formula in 'npudensbw'
dPredict <- function(x, npdensityObject, varName){
  newData <- data.frame(x)
  colnames(newData) <- varName
  return(predict(npdensityObject, newdata=newData))
}

dataDir <- "h:/SCHARP/HVTN703/sievePower/data"

data <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))

ic50 <- data$ic50.geometric.mean.imputed
log10ic50 <- data$ic50.geometric.mean.imputed.log10

# nonparametric density estimation using the package 'np'
densbw <- npudensbw(~ log10ic50, ckertype="epanechnikov")
dens <- npudens(densbw)

# this illustrates how the nonparametric density can be integrated over
# also a sanity check that this is a true density, i.e., the integral = 1
integrate(dPredict, lower=-Inf, upper=Inf, npdensityObject=dens, varName="log10ic50")

hist(log10ic50, breaks=13, freq=FALSE)
x <- seq(-2, 2, by=0.1)
lines(x, dPredict(x, dens, "log10ic50"), col="red")