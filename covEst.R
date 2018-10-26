### estimator for \cov(\hat{\phi},\hat{\lambda}) using Theorem 1 in Juraska and Gilbert (2013, Biometrics)

library(survival)

covEst <- function(time, find, mark, txind, phi.hat, lambda.hat, gamma.hat){
	n <- length(time)
	m <- sum(find)
	time.f <- time[find==1]
	V.f <- cbind(1,mark[find==1])
	txind.f <- txind[find==1]
	time.fM <- matrix(time.f, nrow=n, ncol=m, byrow=TRUE)
	VV <- apply(V.f,1,tcrossprod)
	nmark <- NCOL(V.f)

	g <- function(phi){ exp(drop(V.f %*% phi)) }
	dG <- function(phi){ t(g(phi) * V.f) }
	d2G <- function(phi){ array(t(t(VV)*g(phi)),dim=c(nmark,nmark,m)) }
	dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m)) }

	score1.vect <- function(phi, lambda){
		t((-lambda/(1+lambda*(g(phi)-1)) + txind.f/g(phi)) * t(dG(phi)))
	}
	xi <- function(gamma){ crossprod(time>=time.fM, txind*exp(gamma*txind)) }
	zeta <- function(gamma){ crossprod(time>=time.fM, exp(gamma*txind)) }
	eta <- drop(xi(gamma.hat)/zeta(gamma.hat))             
	score3.vect <- function(gamma){ txind.f-eta }
	l.vect <- function(gamma){
		survprob.vect <- c(1, summary(survfit(Surv(time,find)~1), times=sort(time.f))$surv)
		surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
		time.fMsq <- time.fM[1:m,]
		crossprod(time.f>=time.fMsq, surv.increm*(txind.f*exp(gamma*txind.f) - eta*exp(gamma*txind.f))/zeta(gamma))
	}
	score1 <- function(phi, lambda){
		drop(-lambda * dG(phi) %*% (1/(1+lambda*(g(phi)-1))) + dG(phi) %*% (txind/g(phi)))
	}
	score2 <- function(phi, lambda){
		-sum((g(phi)-1)/(1+lambda*(g(phi)-1)))
	}
	score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
	jack11 <- function(phi, lambda){
		d2Gperm <- aperm(d2G(phi), c(3,1,2))
		dGdGperm <- aperm(dGdG(phi), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(1+lambda*(g(phi)-1))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(1+lambda*(g(phi)-1))^2), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(txind.f/g(phi)), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(txind.f/g(phi)^2), c(2,3,1)),c(1,2),sum)
		-lambda*(term1 - lambda*term2) + term3 - term4
	}
	jack21 <- function(phi, lambda){
		drop(-dG(phi) %*% (1/(1+lambda*(g(phi)-1))^2))
	}
	jack22 <- function(phi, lambda){
		sum(((g(phi)-1)/(1+lambda*(g(phi)-1)))^2)
	}
	jack <- function(phi, lambda){
		j21 <- jack21(phi,lambda)
		(cbind(rbind(jack11(phi,lambda),j21),c(j21,jack22(phi,lambda))))/n
	}
	jack33 <- sum(eta*(eta-1))/n
	
	p <- mean(find)
	omega <- drop(score1.vect(phi.hat,lambda.hat) %*% (score3.vect(gamma.hat) + p*l.vect(gamma.hat))/n - 
	sum(score3.vect(gamma.hat) + p*l.vect(gamma.hat))*apply(score1.vect(phi.hat,lambda.hat),1,sum)/(n^2))   # a vector with 2 components
	drop(solve(jack(phi.hat,lambda.hat))[1:2,1:2] %*% omega)/(n*jack33)
}

### estimator for \cov(\hat{\phi}_{ipw},\hat{\lambda}) based on Theorem 1 of JG (2013)
covEstIPW <- function(time, find, mark, txind, aux.miss, phi.hat, lambda.hat, gamma.hat){
	n <- length(time)
	m.complete <- sum(find==1 & !is.na(mark))
	m.f <- sum(find==1)
	time.f <- time[find==1]
	time.fM <- matrix(time.f, nrow=n, ncol=m.f, byrow=TRUE)
	time.complete <- time[find==1 & !is.na(mark)]
	V.f <- cbind(1,mark[find==1])
	V.complete <- na.omit(V.f)
	na.idx <- attr(V.complete,"na.action")	
	txind.complete <- txind[find==1 & !is.na(mark)]
	txind.f <- txind[find==1]
	aux.miss.f <- aux.miss[find==1]
	time.completeM <- matrix(time.complete, nrow=n, ncol=m.complete, byrow=TRUE)
	VV.complete <- apply(V.complete,1,tcrossprod)
	nmark <- NCOL(V.complete)

	g <- function(phi){ exp(drop(V.complete %*% phi)) }
	dG <- function(phi){ t(g(phi) * V.complete) }
	d2G <- function(phi){ array(t(t(VV.complete)*g(phi)),dim=c(nmark,nmark,m.complete)) }
	dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m.complete)) }

	score1.vect <- function(phi, lambda){
		vect <- matrix(0, nrow=nmark, ncol=m.f)
		vect[,-na.idx] <- t((-lambda/(pi*(1+lambda*(g(phi)-1))) + txind.complete/(pi*g(phi))) * t(dG(phi)))
		vect
	}
	xi <- function(gamma){ crossprod(time>=time.fM, txind*exp(gamma*txind)) }
	zeta <- function(gamma){ crossprod(time>=time.fM, exp(gamma*txind)) }
	eta <- drop(xi(gamma.hat)/zeta(gamma.hat))             
	score3.vect <- function(gamma){ txind.f-eta }
	l.vect <- function(gamma){
		survprob.vect <- c(1, summary(survfit(Surv(time,find)~1))$surv)
		surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
		time.fMsq <- time.fM[1:m.f,]
		crossprod(time.f>=time.fMsq, surv.increm*(txind.f*exp(gamma*txind.f) - eta*exp(gamma*txind.f))/zeta(gamma))
	}
	score1 <- function(phi, lambda){
		drop(-lambda * dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1)))) + dG(phi) %*% (txind.complete/(pi*g(phi))))
	}
	score2 <- function(phi, lambda){
		-sum((g(phi)-1)/(pi*(1+lambda*(g(phi)-1))))
	}
	score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
	jack11 <- function(phi, lambda){
		d2Gperm <- aperm(d2G(phi), c(3,1,2))
		dGdGperm <- aperm(dGdG(phi), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(phi)-1)))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(phi)-1))^2)), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(txind.complete/(pi*g(phi))), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(txind.complete/(pi*g(phi)^2)), c(2,3,1)),c(1,2),sum)
		-lambda*(term1 - lambda*term2) + term3 - term4
	}
	jack21 <- function(phi, lambda){
		drop(-dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1))^2)))
	}
	jack22 <- function(phi, lambda){
		sum(((g(phi)-1)^2)/(pi*(1+lambda*(g(phi)-1))^2))
	}
	jack <- function(phi, lambda){
		j21 <- jack21(phi,lambda)
		(cbind(rbind(jack11(phi,lambda),j21),c(j21,jack22(phi,lambda))))/n
	}
	jack33 <- sum(eta*(eta-1))/n

	r <- apply(V.f, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
	pi.all <- glm(r ~ txind.f*aux.miss.f, family=binomial)$fitted
	if (!is.null(na.idx)){
		pi <- pi.all[-na.idx]
	} else {
		pi <- pi.all
	}
	if (any(pi<0.005)){ stop("Selection probabilities not bounded away from 0.") }

	p <- mean(find==1)
	omega <- drop(score1.vect(phi.hat,lambda.hat) %*% (score3.vect(gamma.hat) + p*l.vect(gamma.hat))/n - 
	sum(score3.vect(gamma.hat) + p*l.vect(gamma.hat))*apply(score1.vect(phi.hat,lambda.hat),1,sum)/(n^2))   # a vector with 2 components
	drop(solve(jack(phi.hat,lambda.hat))[1:2,1:2] %*% omega)/(n*jack33)
}

### estimator for \cov(\hat{\phi}_{aug},\hat{\lambda}) based on Theorem 1 of JG (2013)
covEstAUG <- function(time, find, mark, txind, aux.miss, aux, phi.hat, lambda.hat, gamma.hat){
	n <- length(time)
	m.complete <- sum(find==1 & !is.na(mark))
	m.f <- sum(find==1)
	time.f <- time[find==1]
	time.fM <- matrix(time.f, nrow=n, ncol=m.f, byrow=TRUE)
	time.complete <- time[find==1 & !is.na(mark)]
	V.f <- cbind(1,mark[find==1])
	V.complete <- na.omit(V.f)
	na.idx <- attr(V.complete,"na.action")	
	txind.complete <- txind[find==1 & !is.na(mark)]
	txind.f <- txind[find==1]
	aux.complete <- aux[find==1 & !is.na(mark)]
	aux.f <- aux[find==1]
	aux.miss.f <- aux.miss[find==1]
	time.completeM <- matrix(time.complete, nrow=n, ncol=m.complete, byrow=TRUE)
	VV.complete <- apply(V.complete,1,tcrossprod)
	nmark <- NCOL(V.complete)

	g <- function(phi){ exp(drop(V.complete %*% phi)) }
	dG <- function(phi){ t(g(phi) * V.complete) }
	d2G <- function(phi){ array(t(t(VV.complete)*g(phi)),dim=c(nmark,nmark,m.complete)) }
	dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m.complete)) }

	score1.complete.vect <- function(phi, lambda){
		(-lambda/(1+lambda*(g(phi)-1)) + txind.complete/g(phi)) * t(dG(phi))
	}
	score2.complete.vect <- function(phi, lambda){
		-(g(phi)-1)/(1+lambda*(g(phi)-1))
	}
	aug.mean1 <- function(phi, lambda){
		U <- score1.complete.vect(phi, lambda)
		predicted.vals <- sapply(1:NCOL(U), function(col){
			fit <- lm(U[,col] ~ txind.complete*aux.complete + I(aux.complete^2))
			predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
		})
		predicted.vals
	}
	aug.mean2 <- function(phi, lambda){
		U <- score2.complete.vect(phi, lambda)
		fit <- lm(U ~ txind.complete*aux.complete + I(aux.complete^2))
		predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
	}
	score1.vect <- function(phi, lambda){
		vect <- matrix(0, nrow=m.f, ncol=nmark)
		vect[-na.idx,] <- (-lambda/(pi*(1+lambda*(g(phi)-1))) + txind.complete/(pi*g(phi))) * t(dG(phi))
		vect <- vect + (aug.mean1(phi, lambda) * (1-r/pi.all))
		t(vect)
	}
	xi <- function(gamma){ crossprod(time>=time.fM, txind*exp(gamma*txind)) }
	zeta <- function(gamma){ crossprod(time>=time.fM, exp(gamma*txind)) }
	eta <- drop(xi(gamma.hat)/zeta(gamma.hat))             
	score3.vect <- function(gamma){ txind.f-eta }
	l.vect <- function(gamma){
		survprob.vect <- c(1, summary(survfit(Surv(time,find)~1))$surv)
		surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
		time.fMsq <- time.fM[1:m.f,]
		crossprod(time.f>=time.fMsq, surv.increm*(txind.f*exp(gamma*txind.f) - eta*exp(gamma*txind.f))/zeta(gamma))
	}
	score1 <- function(phi, lambda){
		drop(-lambda * dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1)))) + dG(phi) %*% (txind.complete/(pi*g(phi))) +
		t(aug.mean1(phi, lambda)) %*% (1-r/pi.all))
	}
	score2 <- function(phi, lambda){
		-sum((g(phi)-1)/(pi*(1+lambda*(g(phi)-1)))) + sum(aug.mean2(phi, lambda)*(1-r/pi.all))
	}
	score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
	jack11 <- function(phi, lambda){
		d2Gperm <- aperm(d2G(phi), c(3,1,2))
		dGdGperm <- aperm(dGdG(phi), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(phi)-1)))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(phi)-1))^2)), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(txind.complete/(pi*g(phi))), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(txind.complete/(pi*g(phi)^2)), c(2,3,1)),c(1,2),sum)

		d2U.phi1 <- aperm(d2Gperm*(1/(1+lambda*(g(phi)-1))), c(2,3,1))
		d2U.phi2 <- aperm(dGdGperm*(1/(1+lambda*(g(phi)-1))^2), c(2,3,1))
		d2U.phi3 <- aperm(d2Gperm*(txind.complete/g(phi)), c(2,3,1))
		d2U.phi4 <- aperm(dGdGperm*(txind.complete/g(phi)^2), c(2,3,1))
		d2U.phi <- -lambda*(d2U.phi1 - lambda*d2U.phi2) + d2U.phi3 - d2U.phi4
		predicted.vals.jack <- array(0, dim=c(nmark,nmark,m.f))
		for (i in 1:nmark){
		for (j in 1:nmark){
			resp <- d2U.phi[i,j,]
			fit <- lm(resp ~ txind.complete*aux.complete + I(aux.complete^2))
			predicted.vals.jack[i,j,] <- predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
		}}
		weighted.predicted.vals.jack <- apply(aperm(aperm(predicted.vals.jack, c(3,1,2))*(1-r/pi.all), c(2,3,1)), c(1,2), sum)
		-lambda*(term1 - lambda*term2) + term3 - term4 + weighted.predicted.vals.jack
	}
	jack21 <- function(phi, lambda){
		d2U.phi.lambda <- -t(dG(phi)) * (1/(1+lambda*(g(phi)-1))^2)
		predicted.vals.jack <- matrix(0,nrow=m.f,ncol=nmark)
		for (i in 1:nmark){
			resp <- d2U.phi.lambda[,i]
			fit <- lm(resp ~ txind.complete*aux.complete + I(aux.complete^2))
			predicted.vals.jack[,i] <- predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
		}
		weighted.predicted.vals.jack <- colSums(predicted.vals.jack*(1-r/pi.all))
		drop(-dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1))^2))) + weighted.predicted.vals.jack
	}
	jack22 <- function(phi, lambda){
		d2U.lambda <- ((g(phi)-1)/(1+lambda*(g(phi)-1)))^2
		fit <- lm(d2U.lambda ~ txind.complete*aux.complete + I(aux.complete^2))
		predicted.vals.jack <- predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
		weighted.predicted.vals.jack <- sum(predicted.vals.jack*(1-r/pi.all))
		sum(((g(phi)-1)^2)/(pi*(1+lambda*(g(phi)-1))^2)) + weighted.predicted.vals.jack
	}
	jack <- function(phi, lambda){
		j21 <- jack21(phi,lambda)
		(cbind(rbind(jack11(phi,lambda),j21),c(j21,jack22(phi,lambda))))/n
	}
	jack33 <- sum(eta*(eta-1))/n

	r <- apply(V.f, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
	pi.all <- glm(r ~ txind.f*aux.miss.f, family=binomial)$fitted
	if (!is.null(na.idx)){
		pi <- pi.all[-na.idx]
	} else {
		pi <- pi.all
	}
	if (any(pi<0.005)){ stop("Selection probabilities not bounded away from 0.") }

	p <- mean(find==1)
	omega <- drop(score1.vect(phi.hat,lambda.hat) %*% (score3.vect(gamma.hat) + p*l.vect(gamma.hat))/n - 
	sum(score3.vect(gamma.hat) + p*l.vect(gamma.hat))*apply(score1.vect(phi.hat,lambda.hat),1,sum)/(n^2))   # a vector with 2 components
	drop(solve(jack(phi.hat,lambda.hat))[1:2,1:2] %*% omega)/(n*jack33)
}

### Usage:
# covEst(X,d,V,Z,dRatio$coef[1:2],dRatio$coef[3],gammaHat)
# covEstIPW(X,d,V,Z,A.miss,dRatio$coef[1:2],dRatio$coef[3],gammaHat)
# covEstAUG(X,d,V,Z,A.miss,A.miss,dRatio$coef[1:2],dRatio$coef[3],gammaHat)
		








