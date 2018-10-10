source("markHR.R")

library(survival)

VE <- function(v, a, b, g){ 1-exp(a+b*v+g) }
dVE <- function(v, a, b, g){ -exp(a+b*v+g)*c(1,v,1) }
seVE <- function(v, Var, a, b, g){
  sapply(v, function(mark){ drop(sqrt(t(dVE(mark,a,b,g)) %*% Var %*% dVE(mark,a,b,g))) })
}
Ymat <- function(X,u,ssize) drop(X>=matrix(u,ssize,length(u),byrow=TRUE))
g <- function(theta,V){ exp(drop(V %*% theta)) }
dG <- function(theta,V){ t(g(theta,V) * V) }

# Calculation of lambdaT and lambdaC:
# # E.g. for 704:
# Np = 900
# np = 34
# taumax = 80  # weeks
calcRates <- function(Np, np, taumax) {
  lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
  lambdaC <- 0.1*Np/np*lambdaT
  return(c(lambdaT, lambdaC))
}

# n = number of subjects in vaccine group
simulOne <- function(beta,gamma,n){
  Z <- c(rep(0, n), rep(1, 2*n))   # treatment group
  T0 <- rexp(n, lambdaT)    # failure time for placebo
  T1 <- rexp(2*n, lambdaT*exp(gamma))  # failure time for vaccine
  T <- c(T0,T1)   # failure times
  C <- rexp(3*n, lambdaC)   # censoring times
  X <- pmin(T,C,3)   # observed time is minimum of failure, centoring, and study time
  d <- ifelse(T<=pmin(C,3),1,0)   # failure indicator (0 if censored)
  nInf0 <- sum(d*(1-Z))  # number of infected in placebo group
  nInf1 <- sum(d*Z)      # number of infected in vaccine group
  V0 <- log(1-(1-exp(-2))*runif(n))/(-2)
  V1 <- log(1+(beta-2)*(1-exp(-2))*runif(n)/(2*exp(alpha(beta))))/(beta-2)
  V <- c(V0,V1)
  dRatio <- densRatio(V[d==1],Z[d==1])
  
  if (dRatio$conv){
    
    phReg <- coxph(Surv(X,d)~Z)
    thetaHat <- dRatio$coef
    vthetaHat <- dRatio$var[1:2,1:2]
    gammaHat <- phReg$coef
    vgammaHat <- drop(phReg$var)
    
    ve <- VE(v,thetaHat[1],thetaHat[2],gammaHat)
    
    ### likelihood ratio test of the null hypothesis that beta=0
    lrBeta.pval <- LRtest(V[d==1],Z[d==1],thetaHat[-length(thetaHat)],thetaHat[length(thetaHat)])$pval
    
    return( c(beta, gamma, nInf0, nInf1, ve, lrBeta.pval, thetaHat[-length(thetaHat)], gammaHat) )
    
  }
}