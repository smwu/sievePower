# 've' returns vaccine efficacy values given parameters 'alpha', 'beta', and 'gamma' 
ve <- function(v, a, b, g){ 1-exp(a+b*v+g) }

# 'dPredict' returns nonparametric density estimates at 'x'
# 'npdensityObject' is the output object from 'npudens'
# 'varName' is the variable name used in the formula in 'npudensbw'
dPredict <- function(x, npdensityObject, varName){
  newData <- data.frame(x)
  colnames(newData) <- varName
  return(predict(npdensityObject, newdata=newData))
}

# 'f0' is the density of the mark variable in the placebo group
# density is given for values of `v`
f0 <- function(v, dens, varName){
  dPredict(v, dens, varName)
}

# 'f1' is the density of the mark variable in the vaccine group
# density if given for values of `v`
f1 <- function(v, dens, varName, alpha, beta){
  f0(v, dens, varName)*exp(alpha + beta*v)
}

# 'intf1' is the equation that is solved to obtain parameter 'alpha'
intf1 <- function(alpha, dens, varName, beta){
  integrate(f1, lower=-100, upper=100, dens=dens, varName=varName, alpha=alpha, beta=beta, subdivisions=2000, stop.on.error = FALSE)$value - 1
}

# 'getAlphaBetaGamma' calculates the hazard ratio model coefficients alpha, beta, and gamma for
# a given value of markVE = VE(0.3) = 1 - exp(alpha + beta*0.3 + gamma)
getAlphaBetaGamma <- function(markVE, dens, varName){
  beta <- log(1 - markVE)/log10(0.3/50)
  alpha <- uniroot(intf1, interval=c(-20,20), dens=dens, varName=varName, beta=beta)$root
  # Stephanie, I don't see why 'phi' is equal to the below expression
  # phi <- (log(1 - markVE) - log10(0.3*50)*beta)/2
  # because phi + beta * log10(50) must be equal to 0, I suggest to calculate 'phi' as follows:
  phi <- -beta*log10(50)
  # phi is the sum of the density ratio intercept 'alpha' and the overall hazard ratio 'gamma'
  gamma <- phi - alpha
  return(list(alpha=alpha, beta=beta, gamma=gamma))
}

# `simulOne` simulates data and runs 1-sided Wald-type tests for H00: VE(v)=0 and H0: VE(v)=VE and a likelihood ratio test for H0: VE(v)=VE 
# The function returns a list containing the simulated mark variable, the p-values of the hypothesis tests, 
# and the parameters estimated in the mark-density ratio and the marginal hazard function 
  # 'Np' is number of subjects in placebo group
  # 'np' number of cases in placebo group (fixed at around 34)
  # 'markVE' is the upper vaccine efficacy cutoff, defined as the vaccine efficacy achieved when the mark variable equals 0.3
  # 'taumax' is the follow-up time (in weeks)
  # 'dens' is the output object from the 'npudens' function in the 'np' package
  # 'varName' is the variable name used in the formula in 'npudensbw'
  # 'randomRatio' is the randomization ratio of treatment to placebo (e.g. '2' for 2:1 treatment:placebo randomization)
simulOne <- function(Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, randomRatio, seed){
  Z <- c(rep(0, Np), rep(1, randomRatio*Np))        # treatment group
  T0 <- rexp(Np, lambdaT)                           # failure times for placebo
  T1 <- rexp(randomRatio*Np, lambdaT*exp(gamma))    # failure times for vaccine
  T <- c(T0,T1)                                     # failure times
  C <- rexp((randomRatio+1)*Np, lambdaC)            # censoring times
  X <- pmin(T,C, taumax)                            # observed time: minimum of failure, censoring, and study time
  d <- ifelse(T <= pmin(C,taumax),1,0)              # failure indicator (0 if censored)
  nInf0 <- sum(d*(1-Z))                             # number of infected in placebo group
  nInf1 <- sum(d*Z)                                 # number of infected in vaccine group
  
  Xpoints <- seq(-3,3,len=25000)                    # fine grid of points ranging from -3 to 3 to be sampled from
  prob0 <- f0(Xpoints, dens, varName)               # sampling probability for placebos, using nonparametric density estimates
  prob1 <- f1(Xpoints, dens, varName, alpha, beta)  # sampling probabiliy for vaccinees, using nonparametric density estimates 
  V0 <- sample(Xpoints, size=Np, prob=prob0)        # sample with replacement with probability prob0 to simulate mark in placebos 
  V1 <- sample(Xpoints, size=randomRatio*Np, prob=prob1)  # sample with replacement with probability prob1 to simulate mark in vaccinees
  V <- c(V0, V1)                                    # mark variable
  V <- ifelse(V < log10(0.00076), log10(0.00076), ifelse(V <= log10(50), V, log10(50)))  # mark variable with extreme values censored
    
    # Vraw0 <- rexp(Np, lambdaVraw)  # raw mark in placebos
    # Vraw1 <- rexp(2*Np, (lambdaVraw - beta))  # raw mark in vaccinees
    # Vraw <- c(Vraw0,Vraw1)
  
  # calculate mark density ratio
  # mark is only observed for cases
  dRatio <- densRatio(V[d==1],Z[d==1])  
  
  if (dRatio$conv) {
    
    phReg <- coxph(Surv(X,d) ~ Z)  # calculate marginal hazard ratio
    thetaHat <- dRatio$coef        # estimates for betaHat and alphaHat
    gammaHat <- phReg$coef   
    
    # vaccine efficacies calculated using estimated parameters
    VE <- ve(V,thetaHat[1],thetaHat[2],gammaHat) 
    
    # variance and covariance estimates
    vthetaHat <- dRatio$var[1:2,1:2]
    vgammaHat <- drop(phReg$var) 
    covThG <- covEst(X,d,V,Z,thetaHat[1:2],thetaHat[3],gammaHat)
    
    # one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
    weighted.waldH00 <- (thetaHat[2]/vthetaHat[2,2] - gammaHat/vgammaHat)/
      sqrt(1/vthetaHat[2,2] + 1/vgammaHat - 2*covThG[2]/(vthetaHat[2,2]*vgammaHat))
    weighted.waldH00.pval <- 1 - pnorm(weighted.waldH00)
    
    # 1-sided test of H0: VE(v)=VE vs. alternative that beta > 0
    waldH0 <- thetaHat[2]/sqrt(vthetaHat[2,2])
    waldH0.pval <- 1 - pnorm(waldH0)
    
    # 2-sided likelihood ratio test of the null hypothesis H0: beta=0 (constant VE curve)
    lrBeta.pval <- LRtest(V[d==1],Z[d==1],thetaHat[-length(thetaHat)],thetaHat[length(thetaHat)])$pval
    
    return( list(nInf0 = nInf0, nInf1 = nInf1, alpha = alpha, beta = beta, gamma = gamma, alphaHat=thetaHat[1], 
                 betaHat=thetaHat[2], gammaHat = gammaHat, H00waldP1sided = weighted.waldH00.pval, H0waldP1sided = waldH0.pval, 
                 H0lrP2sided = lrBeta.pval))
  }
}

paste1 <- function(x, y){
  return(paste0(x,"_",y))
}

# 'getInferenceOneMC' returns a matrix of p-values for every combination of the levels of (trial x randRatio x ICpcent x VEcoord) for a single MC iteration
# 'seed' sets a random seed for data generation in 'simulOne'
# 'trial' takes on character strings "704", "703", "704and703"
# 'randRatio' takes on "2:1" and "1:1"
# 'IC' takes on character strings "IC50" and "IC80"
# 'VEcoord' takes on numeric values in [0,1] representing VE(log10(0.3))
# 'taumax' is the follow-up time (in weeks)
# 'alpha1sided' is a 1-sided alpha level
# 'dataB', 'dataC', and 'dataBandC' are the data frames for 704, 703, and 704and703 calculations 
getInferenceOneMC <- function(seed, trial, randRatio, IC, VEcoord, taumax, alpha1sided, dataB, dataC, dataBandC){
  scenario <- c(outer(c(outer(c(outer(trial, randRatio, paste1)), IC, paste1)), VEcoord, paste1))
  out <- vector("list", length=length(scenario))
  names(out) <- scenario
  
  for (i in 1:length(trial)){
    if (trial[i]=="704"){ 
      data <- dataB
      # 'Np' is number of subjects in placebo group
      Np <- 900
      # 'np' is the expected number of cases in placebo group
      np <- 34
    } else if (trial[i]=="703"){ 
      data <- dataC
      Np <- 634
      np <- 34
    } else {
      data <- dataBandC
      Np <- 900 + 634
      np <- 68
    }
    
    # rate parameters for failure time T and censoring time C
    lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
    lambdaC <- 0.1*Np/np*lambdaT
    
    for (k in 1:length(IC)){
      mark <- data[,paste0(tolower(IC[k]),".geometric.mean.imputed.log10")]
      
      # nonparametric kernel density estimate of the mark variable
      densbw <- npudensbw(~ mark, ckertype="epanechnikov")
      dens <- npudens(densbw)
      
      for (l in 1:length(VEcoord)){
        coeff <- getAlphaBetaGamma(VEcoord[l], dens, "mark")
        
        for (j in 1:length(randRatio)){
          scenarioLabel <- paste0(trial[i],"_",randRatio[j],"_",IC[k],"_",VEcoord[l])
          out[[scenarioLabel]] <- simulOne(Np=Np, np=np, lambdaT=lambdaT, lambdaC=lambdaC, alpha=coeff$alpha, beta=coeff$beta, gamma=coeff$gamma, 
                                           taumax=taumax, dens=dens, varName="mark", randomRatio=ifelse(randRatio[j]=="2:1",2,1), seed=seed)
        }
      }
    }
  }
  
  return(out)
}

# 'calcPower' calculates power for the 1-sided Wald-type test of H00:VE(v)=0, the 1-sided Wald-type test of H0:VE(v)=VE,
# the 1-sided likelihood ratio test of H0:beta=0, H1:beta>0, and the 2-sided likelihood ratio test of H0:beta=0.
# Given a dataframe, 'power', with columns 'WaldH00', 'WaldH0', 'LR', and 'twosidedLR' and each row representing a 
# different scenario (e.g., AMP-B, VE(0.3)=0.7, IC50 could characterize one scenario), the function modifies the 
# dataframe and returns the modified dataframe. 
  # 'index' is the index of the row (and scenario) that will be modified and contain the new power calculations
  # 'alphaLR' is the type 1 error rate for the likelihood ratio tests
  # 'alphaWald' is the type 1 error rate for the Wald-type tests
calcPower <- function(index, simulations, Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, alphaLR, alphaWald, power, randomRatio) {
  for (i in 1:simulations) {
    result <- simulOne(Np, np, lambdaT, lambdaC, alpha, beta, gamma, taumax, dens, varName, randomRatio)
    
    # one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
    if (result$weighted.waldH00.pval <= alphaWald) {
      power$WaldH00[index] <- power$WaldH00[index] + 1
    }
    
    # 1-sided Wald-type test of H0: VE(v)=VE vs. alternative that beta > 0
    if (result$waldH0.pval <= alphaWald) {
      power$WaldH0[index] <- power$WaldH0[index] + 1
    }
    
    # 1-sided likelihood ratio test (serves as sanity check)
    if ((result$lrBeta.pval <= alphaLR) & (result$betaHat > 0)) { 
      power$LR[index] <- power$LR[index] + 1
    }
    
    # two-sided likelihood ratio p-value
    if (result$lrBeta.pval <= alphaLR) {
      power$twosidedLR[index] <- power$twosidedLR[index] + 1
    }
  }
  
  return(power)
}


#=================================Extraneous functions=================================================================
# # 'plotDistHazVE' generates a plot with point and interval estimates of VE(v), and scatter- and box-plots of univariate index values by treatment in the top panel
# # p-values of tests evaluating H00 and H0 are reported
# # the number of distances is reported
# # 'nSites' is the number of AA sites in the vaccine insert sequence included in calculating the distance
# # distData is a dataframe with columns tx (vaccine=1, placebo=0) and mark (positive, NA values omitted)
# # 'out' is the output from distHazVE
# # 'raw' specifies if the mark varialbe is on the raw scale (TRUE) or not (FALSE)
# simulPlotDistHazVE <- function(result, nSites=NULL, xlab, ylab, title, raw){
#   
#   cexAxis <- 0.75
#   cexLab <- 0.8
#   cexTitle <- 1
#   cexText <- 0.85
#   cexLegend <- 0.8
#   
#   # par(mar=c(3,3,3,2.5), oma=c(0,0,0,0), cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)  
#   par(cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)
#   if (raw == TRUE) {
#     plot(result$Vraw, result$VEraw, xlim=range(result$Vraw), ylim=c(0,1), type = "n", xlab="", ylab="", bty="l", main=title)  ###yLim changed
#     mtext(xlab, side=1, line=2, cex=cexLab)
#     mtext(ylab, side=2, line=2.15, las=3, cex=cexLab)
#     abline(h=0, col="gray30", lty="dotted", lwd=2)
#     points(result$Vraw, result$VEraw)
#     ### continuous mark-specific treatment efficacy estimates
#     text(0.2, 0.2, paste0("No. of Cases (V:P): ",sum(result$trt.id==1),":",sum(result$trt.id==0), "\nLR test of H0: \nPE(v) = PE: ",
#                 ifelse(result$lrBeta.pval<0.001,"p < 0.001", paste0("p = ",format(result$lrBeta.pval,digits=2)))), pos=4, cex=cexText)
#     
#   } else {
#     plot(result$V, result$VE, xlim=range(result$V), ylim=c(0,1), type = "n", xaxt = "n", xlab="", ylab="", bty="l", main=title)  ###yLim changed
#     grid <- seq(min(result$V), max(result$V), length = 6)
#     axis(side=1, at=grid, labels=round(10^grid,3), las=1)
#     mtext(xlab, side=1, line=2, cex=cexLab)
#     mtext(ylab, side=2, line=2.15, las=3, cex=cexLab)
#     abline(h=0, col="gray30", lty="dotted", lwd=2)
#     points(result$V, result$VE)
#     ### continuous mark-specific treatment efficacy estimates
#     text(x=grid[2], 0.2, paste0("No. of Cases (V:P): ",sum(result$tx==1),":",sum(result$tx==0), "\nLR test of H0: \nPE(v) = PE: ",
#                 ifelse(result$lrBeta.pval<0.001,"p < 0.001", paste0("p = ",format(result$lrBeta.pval,digits=2)))), cex=cexText)
#     
#   }
#   
#   if (!is.null(nSites)){mtext(paste0("Number of vaccine insert AA sites = ",nSites), side=1, line=-1, cex=0.8, adj=1)}
# }


# # Calculation of lambdaT and lambdaC
# calcRates <- function(Np, np, taumax) {
#   lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
#   lambdaC <- 0.1*Np/np*lambdaT
#   return(c(lambdaT, lambdaC))
# }
# 
# # Calculation of true alpha, beta, and gamma parameters
# calcTrueParams <- function (lambdaVraw, markVE) {
#   beta <- -log(1-markVE)/49.5
#   alpha <- log(1 - beta/lambdaVraw)  # alpha = log(1 + (log(1-markVE)/(49.5*lambdaVraw)))
#   phi <- (log(1 - markVE) - 50.5*beta)/2  # phi = 48/49.5*log(1-markVE)
#   gamma <- phi - alpha  # gamma = log(1-markVE)*48/49.5 - log(1 + log(1-markVE)/(49.5*lambdaVraw))
#   return(c(beta, alpha, gamma))
# }
