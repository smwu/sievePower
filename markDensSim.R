source("markHR.R")
library(survival)

# Calculation of lambdaT and lambdaC
calcRates <- function(Np, np, taumax) {
  lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
  lambdaC <- 0.1*Np/np*lambdaT
  return(c(lambdaT, lambdaC))
}

# Calculation of true alpha, beta, and gamma parameters
calcTrueParams <- function (lambdaVraw, markVE) {
  beta <- -log(1-markVE)/45.5
  alpha <- log(1 - beta/lambdaVraw)  # alpha = log(1 + (log(1-markVE)/(45.5*lambdaVraw)))
  phi <- (log(1 - markVE) - 50.5*beta)/2  # phi = 48/45.5*log(1-markVE)
  gamma <- phi - alpha  # gamma = log(1-markVE)*48/45.5 - log(1 + log(1-markVE)/(45.5*lambdaVraw))
  return(c(beta, alpha, gamma))
}

# calculation of VE
ve <- function(v, a, b, g){ 1-exp(a+b*v+g) }

# Simulate data and run likelihood ratio test
# n = number of subjects in placebo group = Np
simulOne <- function(Np, np, lambdaVraw, markVE, taumax){
  
  lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
  lambdaC <- 0.1*Np/np*lambdaT
  
  beta <- -log(1-markVE)/45.5
  alpha <- log(1 - beta/lambdaVraw)  # alpha = log(1 + (log(1-markVE)/(45.5*lambdaVraw)))
  phi <- (log(1 - markVE) - 50.5*beta)/2  # phi = 48/45.5*log(1-markVE)
  gamma <- phi - alpha
  
  Z <- c(rep(0, Np), rep(1, 2*Np))   # treatment group
  T0 <- rexp(Np, lambdaT)    # failure time for placebo
  T1 <- rexp(2*Np, lambdaT*exp(gamma))  # failure time for vaccine
  T <- c(T0,T1)   # failure times
  C <- rexp(3*Np, lambdaC)   # censoring times
  X <- pmin(T,C, taumax)   # observed time is minimum of failure, centoring, and study time
  d <- ifelse(T <= pmin(C,taumax),1,0)   # failure indicator (0 if censored)
  nInf0 <- sum(d*(1-Z))  # number of infected in placebo group
  nInf1 <- sum(d*Z)      # number of infected in vaccine group
  Vraw0 <- rexp(Np, lambdaVraw)  # raw mark in placebos
  Vraw1 <- rexp(2*Np, (lambdaVraw - beta))  # raw mark in vaccinees
  Vraw <- c(Vraw0,Vraw1)
  V <- ifelse(Vraw < 0.00076, log10(0.00076), ifelse(Vraw <= 50, log10(Vraw), log10(50)))  # scaled mark
  dRatio <- densRatio(V[d==1],Z[d==1])  # mark is only observed for cases
  
  if (dRatio$conv) {
    
    phReg <- coxph(Surv(X,d) ~ Z)
    thetaHat <- dRatio$coef
    gammaHat <- phReg$coef
    
    VE <- ve(V,thetaHat[1],thetaHat[2],gammaHat)
    VEraw <- ve(Vraw, thetaHat[1], thetaHat[2], gammaHat)
    
    ### likelihood ratio test of the null hypothesis that beta=0
    lrBeta.pval <- LRtest(V[d==1],Z[d==1],thetaHat[-length(thetaHat)],thetaHat[length(thetaHat)])$pval
    
    return( list(beta = beta, alpha = alpha, gamma = gamma, lambdaT = lambdaT, lambdaC = lambdaC, 
                 nInf0 = nInf0, nInf1 = nInf1, V = V, Vraw = Vraw, VE = VE, VEraw = VEraw, mark = V[d==1], trt.id = Z[d==1],
                 lrBeta.pval = lrBeta.pval, thetaHat = thetaHat[-length(thetaHat)], gammaHat = gammaHat))
    
    # return( list(beta = beta, gamma = gamma, nInf0 = nInf0, nInf1 = nInf1, V = V[d==1], Vraw = Vraw[d==1], 
    #              distGridRaw = distGridRaw, distGrid = distGrid, VEraw = VEraw, VE = VE, Z = Z[d==1], 
    #              lrBeta.pval = lrBeta.pval, thetaHat = thetaHat[-length(thetaHat)], gammaHat = gammaHat) )
  }
}


# 'plotDistHazVE' generates a plot with point and interval estimates of VE(v), and scatter- and box-plots of univariate index values by treatment in the top panel
# p-values of tests evaluating H00 and H0 are reported
# the number of distances is reported
# 'nSites' is the number of AA sites in the vaccine insert sequence included in calculating the distance
# distData is a dataframe with columns tx (vaccine=1, placebo=0) and mark (positive, NA values omitted)
# 'out' is the output from distHazVE
# 'raw' specifies if the mark varialbe is on the raw scale (TRUE) or not (FALSE)
simulPlotDistHazVE <- function(result, nSites=NULL, xlab, ylab, title, raw){
  
  cexAxis <- 0.75
  cexLab <- 0.8
  cexTitle <- 1
  cexText <- 0.85
  cexLegend <- 0.8
  
  # par(mar=c(3,3,3,2.5), oma=c(0,0,0,0), cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)  
  par(cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)
  if (raw == TRUE) {
    plot(result$Vraw, result$VEraw, xlim=range(result$Vraw), ylim=c(0,1), type = "n", xlab="", ylab="", bty="l", main=title)  ###yLim changed
    mtext(xlab, side=1, line=2, cex=cexLab)
    mtext(ylab, side=2, line=2.15, las=3, cex=cexLab)
    abline(h=0, col="gray30", lty="dotted", lwd=2)
    points(result$Vraw, result$VEraw)
    ### continuous mark-specific treatment efficacy estimates
    text(0.2, 0.2, paste0("No. of Cases (V:P): ",sum(result$trt.id==1),":",sum(result$trt.id==0), "\nLR test of H0: \nPE(v) = PE: ",
                ifelse(result$lrBeta.pval<0.001,"p < 0.001", paste0("p = ",format(result$lrBeta.pval,digits=2)))), pos=4, cex=cexText)
    
  } else {
    plot(result$V, result$VE, xlim=range(result$V), ylim=c(0,1), type = "n", xaxt = "n", xlab="", ylab="", bty="l", main=title)  ###yLim changed
    grid <- seq(min(result$V), max(result$V), length = 6)
    axis(side=1, at=grid, labels=round(10^grid,3), las=1)
    mtext(xlab, side=1, line=2, cex=cexLab)
    mtext(ylab, side=2, line=2.15, las=3, cex=cexLab)
    abline(h=0, col="gray30", lty="dotted", lwd=2)
    points(result$V, result$VE)
    ### continuous mark-specific treatment efficacy estimates
    text(x=grid[2], 0.2, paste0("No. of Cases (V:P): ",sum(result$tx==1),":",sum(result$tx==0), "\nLR test of H0: \nPE(v) = PE: ",
                ifelse(result$lrBeta.pval<0.001,"p < 0.001", paste0("p = ",format(result$lrBeta.pval,digits=2)))), cex=cexText)
    
  }
  
  if (!is.null(nSites)){mtext(paste0("Number of vaccine insert AA sites = ",nSites), side=1, line=-1, cex=0.8, adj=1)}
}

