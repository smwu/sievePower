source("markHR.R")
#source("covEst.R")

library(survival)

VE <- function(v, a, b, g){ 1-exp(a+b*v+g) }
dVE <- function(v, a, b, g){ -exp(a+b*v+g)*c(1,v,1) }
seVE <- function(v, Var, a, b, g){
  sapply(v, function(mark){ drop(sqrt(t(dVE(mark,a,b,g)) %*% Var %*% dVE(mark,a,b,g))) })
}
Ymat <- function(X,u,ssize) drop(X>=matrix(u,ssize,length(u),byrow=TRUE))
g <- function(theta,V){ exp(drop(V %*% theta)) }
dG <- function(theta,V){ t(g(theta,V) * V) }

simulOne <- function(beta,gamma,n){
  ssize <- 2*n
  Z <- rep(0:1, each=n)
  T0 <- rexp(n,lambda)    
  T1 <- rexp(n,lambda*exp(gamma))
  T <- c(T0,T1)
  C <- runif(ssize, 0, 15)
  X <- pmin(T,C,3)
  d <- ifelse(T<=pmin(C,3),1,0)
  nInf0 <- sum(d*(1-Z))
  nInf1 <- sum(d*Z)
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
    covThG <- covEst(X,d,V,Z,thetaHat[1:2],thetaHat[3],gammaHat)
    covThG <- numeric(2)
    
    ve <- VE(v,thetaHat[1],thetaHat[2],gammaHat)
    Sigma <- cbind(rbind(vthetaHat,covThG), c(covThG,vgammaHat))
    se <- seVE(v,Sigma,thetaHat[1],thetaHat[2],gammaHat)
    ci <- c(ve + qnorm(0.975)*se %o% c(-1,1))
    
    ### Wald test of H0: VE(v)=VE
    waldH0 <- thetaHat[2]/sqrt(vthetaHat[2,2])
    
    ### Wald test of H00: VE(v)=0
    waldH00 <- drop(t(c(thetaHat[2], gammaHat)) %*% solve(Sigma[2:3,2:3]) %*% c(thetaHat[2], gammaHat))
    
    ### one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
    weighted.waldH00 <- (thetaHat[2]/vthetaHat[2,2] - gammaHat/vgammaHat)/
      sqrt(1/vthetaHat[2,2] + 1/vgammaHat - 2*covThG[2]/(vthetaHat[2,2]*vgammaHat))
    
    ### log-rank test of the equality of time-to-infection distributions (test statistic and p-value)
    logrank <- survdiff(Surv(X,d)~Z)$chisq
    logrank.pval <- 1-pchisq(logrank,1)
    
    ### likelihood ratio test of the null hypothesis that beta=0
    lrBeta.pval <- LRtest(V[d==1],Z[d==1],thetaHat[-length(thetaHat)],thetaHat[length(thetaHat)])$pval
    
    ### Kolmogorov-Smirnov-type test of conditional independence between T and V given Z
    # indTV.pval <- bootpval(cbind(X,V,d))
    
    return( c(beta, gamma, nInf0, nInf1, ve, se, ci, waldH0, waldH00, weighted.waldH00, logrank.pval, lrBeta.pval, thetaHat[-length(thetaHat)], gammaHat) )
    
  }
}