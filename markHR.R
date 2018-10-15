densRatio <- function(mark, trt.id){
	V <- cbind(1,mark)
	z <- trt.id
	nmark <- NCOL(V)
	ninf <- NROW(V)
	VV <- apply(V,1,tcrossprod)
		
	g <- function(theta){ exp(drop(V %*% theta)) }
	dG <- function(theta){ t(g(theta) * V) }
	d2G <- function(theta){ array(t(t(VV)*g(theta)),dim=c(nmark,nmark,ninf)) }
	dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

	score1 <- function(theta, lambda){
		drop(-lambda * dG(theta) %*% (1/(1+lambda*(g(theta)-1))) + dG(theta) %*% (z/g(theta)))
	}
	score2 <- function(theta, lambda){
		-sum((g(theta)-1)/(1+lambda*(g(theta)-1)))
	}
	score <- function(theta, lambda){ c(score1(theta,lambda),score2(theta,lambda)) }
	jack11 <- function(theta, lambda){
		d2Gperm <- aperm(d2G(theta), c(3,1,2))
		dGdGperm <- aperm(dGdG(theta), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(1+lambda*(g(theta)-1))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(1+lambda*(g(theta)-1))^2), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(z/g(theta)), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(z/g(theta)^2), c(2,3,1)),c(1,2),sum)
		-lambda*(term1 - lambda*term2) + term3 - term4
	}
	jack21 <- function(theta, lambda){
		drop(-dG(theta) %*% (1/(1+lambda*(g(theta)-1))^2))
	}
	jack22 <- function(theta, lambda){
		sum(((g(theta)-1)/(1+lambda*(g(theta)-1)))^2)
	}
	jack <- function(theta, lambda){
		j21 <- jack21(theta,lambda)
		cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
	}

	param.old <- numeric(nmark+1)
	param.new <- c(numeric(nmark),0.5)
	while (sum((param.new - param.old)^2)>1e-8){
		param.old <- param.new
		jackInv <- try(solve(jack(param.old[-(nmark+1)],param.old[nmark+1])), silent=TRUE)
		if (class(jackInv)!="try-error"){
			param.new <- param.old - drop(jackInv %*% score(param.old[-(nmark+1)],
			param.old[nmark+1]))
		}
		if (sum(is.nan(param.new))>0) break
	}
	theta.new <- param.new[-(nmark+1)]
	lambda.new <- param.new[nmark+1]
	
	SigmaHat <- function(theta, lambda){
		L <- -lambda * t(dG(theta)) * (1/(1+lambda*(g(theta)-1))) + t(dG(theta)) * (z/g(theta))
		L <- cbind(L, (g(theta)-1)/(1+lambda*(g(theta)-1)))
		crossprod(L)/ninf
	}

	JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
	if (class(JackInv)!="try-error"){
		Var <- ninf * JackInv %*% SigmaHat(theta.new,lambda.new) %*% JackInv
		names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",
		paste("beta",1:(nmark-1),sep=""),"lambda")
	} else {
		Var <- NULL
	}

	list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new),
	conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}

### goodness-of-fit test of the validity of a density ratio model (Qin and Zhang, 1997)
### can handle univariate and bivariate marks
### B = number of bootstrap iterations
densRatioGOFtest <- function(mark, trt.id, theta=NULL, lambda=NULL, B=1000){
	ninf <- length(trt.id)   ## number of infections
	n0 <- sum(1-trt.id)      ## number of placebo infections
	n1 <- sum(trt.id)        ## number of vaccine infections
	if (sum(is.null(theta),is.null(lambda))>0){
		param <- densRatio(mark, trt.id)$coef
		theta <- param[-length(param)]
		lambda <- param[length(param)]
	}

	g <- function(mark, theta){ exp(drop(cbind(1,mark) %*% theta)) }
	p <- function(mark, theta, lambda, m){ 1/(m*(1+lambda*(g(mark,theta)-1))) }
	F0np <- function(mark, trt.id){
		if (NCOL(as.matrix(mark))==1){
			F0np.vector <- sapply(mark, function(mark.i){ mean(mark[trt.id==0]<=mark.i) })
		} else {
			F0np.vector <- apply(mark, 1, function(mark.i){ mean(mark[trt.id==0,1]<=mark.i[1] & mark[trt.id==0,2]<=mark.i[2]) })
		}
		return( F0np.vector )
	}
	F0sp <- function(mark, prob){
		if (NCOL(as.matrix(mark))==1){
			F0sp.vector <- sapply(mark, function(mark.i){ sum(prob*ifelse(mark<=mark.i,1,0)) })
		} else {
			F0sp.vector <- apply(mark, 1, function(mark.i){ sum(prob*ifelse(mark[,1]<=mark.i[1] & mark[,2]<=mark.i[2],1,0)) })
		}
		return( F0sp.vector )
	}

	f0 <- p(mark,theta,lambda,ninf)
	delta <- sqrt(ninf)*max(abs(F0np(mark,trt.id) - F0sp(mark,f0)))

	N0.ast <- rmultinom(B,n0,f0)
	N1.ast <- rmultinom(B,n1,f0*g(mark,theta))
	v0.ast <- lapply(1:B, function(iter){
		if (NCOL(as.matrix(mark))==1){
			out <- rep(mark,N0.ast[,iter])
		} else {
			out <- cbind(rep(mark[,1],N0.ast[,iter]), rep(mark[,2],N0.ast[,iter]))
		}
		return( out )
	})
	v1.ast <- lapply(1:B, function(iter){
		if (NCOL(as.matrix(mark))==1){
			out <- rep(mark,N1.ast[,iter])
		} else {
			out <- cbind(rep(mark[,1],N1.ast[,iter]), rep(mark[,2],N1.ast[,iter]))
		}
		return( out )
	})
	z.ast <- rep(0:1,c(n0,n1))

	teststat <- sapply(1:B, function(iter){
		v.ast.iter <- rbind(as.matrix(v0.ast[[iter]]),as.matrix(v1.ast[[iter]]))
		param <- densRatio(v.ast.iter, z.ast)$coef
		theta.ast <- param[-length(param)]
		lambda.ast <- param[length(param)]
		sqrt(ninf)*max(abs(F0np(v.ast.iter,z.ast) - F0sp(v.ast.iter,p(v.ast.iter,theta.ast,lambda.ast,ninf))))
	})
	pval <- mean(teststat>=delta)
	list(teststat=delta, pval=pval, theta=theta, lambda=lambda)
}


LRtest <- function(mark, trt.id, theta.hat, lambda.hat){
	V <- cbind(1,mark)
	z <- trt.id
	nmark <- NCOL(V)
		
	g <- function(theta){ exp(drop(V %*% theta)) }
	loglik <- function(theta, lambda){ -sum(log(1 + lambda*(g(theta)-1))) + sum(z*log(g(theta))) }

	teststat <- 2*(loglik(theta.hat, lambda.hat) - loglik(rep(0,nmark), 0))
	pval <- 1-pchisq(teststat,nmark-1)
	list(teststat=teststat, pval=pval)
}

densRatioIPW <- function(mark, trt.id, aux){
	V <- cbind(1,mark)
	V.complete <- na.omit(V)
	z <- trt.id
	na.idx <- attr(V.complete,"na.action")
	if (!is.null(na.idx)){
		z.complete <- z[-na.idx]
	} else {
		z.complete <- z
	}
	nmark <- NCOL(V.complete)
	ninf <- NROW(V.complete)
	VV.complete <- apply(V.complete,1,tcrossprod)
		
	g <- function(theta){ exp(drop(V.complete %*% theta)) }
	dG <- function(theta){ t(g(theta) * V.complete) }
	d2G <- function(theta){ array(t(t(VV.complete)*g(theta)),dim=c(nmark,nmark,ninf)) }
	dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

	score1 <- function(theta, lambda){
		drop(-lambda * dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1)))) + dG(theta) %*% (z.complete/(pi*g(theta))))
	}
	score2 <- function(theta, lambda){
		-sum((g(theta)-1)/(pi*(1+lambda*(g(theta)-1))))
	}
	score <- function(theta, lambda){ c(score1(theta,lambda),score2(theta,lambda)) }
	jack11 <- function(theta, lambda){
		d2Gperm <- aperm(d2G(theta), c(3,1,2))
		dGdGperm <- aperm(dGdG(theta), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(theta)-1)))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(theta)-1))^2)), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(z.complete/(pi*g(theta))), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(z.complete/(pi*g(theta)^2)), c(2,3,1)),c(1,2),sum)
		-lambda*(term1 - lambda*term2) + term3 - term4
	}
	jack21 <- function(theta, lambda){
		drop(-dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1))^2)))
	}
	jack22 <- function(theta, lambda){
		sum(((g(theta)-1)^2)/(pi*(1+lambda*(g(theta)-1))^2))
	}
	jack <- function(theta, lambda){
		j21 <- jack21(theta,lambda)
		cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
	}

	r <- apply(V, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
	# pi.all <- glm(r ~ z + aux + z*aux, family=binomial)$fitted
	pi.all <- glm(r ~ z, family=binomial)$fitted
	if (!is.null(na.idx)){
		pi <- pi.all[-na.idx]
	} else {
		pi <- pi.all
	}
	if (any(pi<0.005)){ stop("Selection probabilities not bounded away from 0.") }

	param.old <- numeric(nmark+1)
	param.new <- c(numeric(nmark),0.5)
	while (sum((param.new - param.old)^2)>1e-8){
		param.old <- param.new
		jackInv <- try(solve(jack(param.old[-(nmark+1)],param.old[nmark+1])), silent=TRUE)
		if (class(jackInv)!="try-error"){
			param.new <- param.old - drop(jackInv %*% score(param.old[-(nmark+1)],param.old[nmark+1]))
		}
		if (sum(is.nan(param.new))>0) break
	}
	theta.new <- param.new[-(nmark+1)]
	lambda.new <- param.new[nmark+1]
	
	Resid <- function(theta, lambda){
		U <- matrix(0,nrow=length(z),ncol=nmark+1)
		U[-na.idx,1:nmark] <- -lambda * t(dG(theta)) * (1/(pi*(1+lambda*(g(theta)-1)))) + t(dG(theta)) * (z.complete/(pi*g(theta)))
		U[-na.idx,nmark+1] <- (g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))
		# S <- (r-pi.all) * cbind(1,z,aux,z*aux)
		S <- (r-pi.all) * cbind(1,z)
		# resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
		resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2])$resid })
		Resids <- do.call("cbind",resids)
		crossprod(Resids)/ninf
	}

	JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
	if (class(JackInv)!="try-error"){
		Var <- ninf * JackInv %*% Resid(theta.new,lambda.new) %*% JackInv
		names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",paste("beta",1:(nmark-1),sep=""),"lambda")
	} else {
		Var <- NULL
	}

	list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), probs=pi, conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}

densRatioAUG <- function(mark, trt.id, aux.miss, aux){
	V <- cbind(1,mark)
	V.complete <- na.omit(V)
	z <- trt.id
	na.idx <- attr(V.complete,"na.action")
	if (!is.null(na.idx)){
		z.complete <- z[-na.idx]
		aux.complete <- aux[-na.idx]
	} else {
		z.complete <- z
		aux.complete <- aux
	}
	nmark <- NCOL(V.complete)
	ninf <- NROW(V.complete)
	VV.complete <- apply(V.complete,1,tcrossprod)
		
	g <- function(theta){ exp(drop(V.complete %*% theta)) }
	dG <- function(theta){ t(g(theta) * V.complete) }
	d2G <- function(theta){ array(t(t(VV.complete)*g(theta)),dim=c(nmark,nmark,ninf)) }
	dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

	fscore.i1 <- function(theta, lambda){
		-lambda * t(dG(theta)) * (1/(1+lambda*(g(theta)-1))) + t(dG(theta)) * (z.complete/g(theta))
	}
	fscore.i2 <- function(theta, lambda){
		-(g(theta)-1)/(1+lambda*(g(theta)-1))
	}
	fscore.i <- function(theta, lambda){
		cbind(fscore.i1(theta, lambda), fscore.i2(theta, lambda))
	}
	aug.mean1 <- function(theta, lambda){
		U <- fscore.i1(theta, lambda)
		predicted.vals <- sapply(1:NCOL(U), function(col){
			fit <- lm(U[,col] ~ z.complete*aux.complete + I(aux.complete^2))
			predict(fit, data.frame(z.complete=z, aux.complete=aux))
		})
		return( predicted.vals )
	}
	aug.mean2 <- function(theta, lambda){
		U <- fscore.i2(theta, lambda)
		fit <- lm(U ~ z.complete*aux.complete + I(aux.complete^2))
		predict(fit, data.frame(z.complete=z, aux.complete=aux))
	}

	score1 <- function(theta, lambda){
		drop(-lambda * dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1)))) + dG(theta) %*% (z.complete/(pi*g(theta))) +
		t(aug.mean1(theta, lambda)) %*% (1-r/pi.all))
	}
	score2 <- function(theta, lambda){
		-sum((g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))) + sum(aug.mean2(theta, lambda)*(1-r/pi.all))
	}
	score <- function(theta, lambda){ c(score1(theta,lambda), score2(theta,lambda)) }
	jack11 <- function(theta, lambda){
		d2Gperm <- aperm(d2G(theta), c(3,1,2))
		dGdGperm <- aperm(dGdG(theta), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(theta)-1)))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(theta)-1))^2)), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(z.complete/(pi*g(theta))), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(z.complete/(pi*g(theta)^2)), c(2,3,1)),c(1,2),sum)

		d2U.theta1 <- aperm(d2Gperm*(1/(1+lambda*(g(theta)-1))), c(2,3,1))
		d2U.theta2 <- aperm(dGdGperm*(1/(1+lambda*(g(theta)-1))^2), c(2,3,1))
		d2U.theta3 <- aperm(d2Gperm*(z.complete/g(theta)), c(2,3,1))
		d2U.theta4 <- aperm(dGdGperm*(z.complete/g(theta)^2), c(2,3,1))
		d2U.theta <- -lambda*(d2U.theta1 - lambda*d2U.theta2) + d2U.theta3 - d2U.theta4
		predicted.vals.jack <- array(0, dim=c(nmark,nmark,length(z)))
		for (i in 1:nmark){
		for (j in 1:nmark){
			resp <- d2U.theta[i,j,]
			fit <- lm(resp ~ z.complete*aux.complete + I(aux.complete^2))
			predicted.vals.jack[i,j,] <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
		}}
		weighted.predicted.vals.jack <- apply(aperm(aperm(predicted.vals.jack, c(3,1,2))*(1-r/pi.all), c(2,3,1)), c(1,2), sum)

		return( -lambda*(term1 - lambda*term2) + term3 - term4 + weighted.predicted.vals.jack )
	}
	jack21 <- function(theta, lambda){
		d2U.theta.lambda <- -t(dG(theta)) * (1/(1+lambda*(g(theta)-1))^2)
		predicted.vals.jack <- matrix(0,nrow=length(z),ncol=nmark)
		for (i in 1:nmark){
			resp <- d2U.theta.lambda[,i]
			fit <- lm(resp ~ z.complete*aux.complete + I(aux.complete^2))
			predicted.vals.jack[,i] <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
		}
		weighted.predicted.vals.jack <- colSums(predicted.vals.jack*(1-r/pi.all))
		return( drop(-dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1))^2))) + weighted.predicted.vals.jack )
	}
	jack22 <- function(theta, lambda){
		d2U.lambda <- ((g(theta)-1)/(1+lambda*(g(theta)-1)))^2
		fit <- lm(d2U.lambda ~ z.complete*aux.complete + I(aux.complete^2))
		predicted.vals.jack <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
		weighted.predicted.vals.jack <- sum(predicted.vals.jack*(1-r/pi.all))
		return( sum(((g(theta)-1)^2)/(pi*(1+lambda*(g(theta)-1))^2)) + weighted.predicted.vals.jack )
	}
	jack <- function(theta, lambda){
		j21 <- jack21(theta,lambda)
		cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
	}

	r <- apply(V, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
	# pi.all <- glm(r ~ z + aux.miss + z*aux.miss, family=binomial)$fitted
	pi.all <- glm(r ~ z, family=binomial)$fitted
	if (!is.null(na.idx)){
		pi <- pi.all[-na.idx]
	} else {
		pi <- pi.all
	}
	if (any(pi==0)){ stop("Selection probabilities not bounded away from 0.") }

	param.old <- numeric(nmark+1)
	param.new <- c(numeric(nmark),0.5)
#	param.old <- densRatioIPW(mark, trt.id, aux)$coef
#	param.new <- param.old + c(numeric(nmark),1e-2)
	while (sum((param.new - param.old)^2)>1e-8){
		param.old <- param.new
		jackInv <- try(solve(jack(param.old[-(nmark+1)],param.old[nmark+1])), silent=TRUE)
		if (class(jackInv)!="try-error"){
			param.new <- param.old - drop(jackInv %*% score(param.old[-(nmark+1)],param.old[nmark+1]))
		}
		if (sum(is.nan(param.new))>0) break
	}
	theta.new <- param.new[-(nmark+1)]
	lambda.new <- param.new[nmark+1]
	
	Resid <- function(theta, lambda){
		U <- matrix(0,nrow=length(z),ncol=nmark+1)
		U[-na.idx,1:nmark] <- -lambda * t(dG(theta)) * (1/(pi*(1+lambda*(g(theta)-1)))) + t(dG(theta)) * (z.complete/(pi*g(theta)))
		U[,1:nmark] <- U[,1:nmark] + aug.mean1(theta, lambda) * (1-r/pi.all)
		U[-na.idx,nmark+1] <- (g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))
		U[,nmark+1] <- U[,nmark+1] + aug.mean2(theta, lambda) * (1-r/pi.all)
		# S <- (r-pi.all) * cbind(1,z,aux.miss,z*aux.miss)
		S <- (r-pi.all) * cbind(1,z)
		# resids <- sapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
		resids <- sapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2])$resid })
		crossprod(resids)/ninf
	}

	JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
	if (class(JackInv)!="try-error"){
		Var <- ninf * JackInv %*% Resid(theta.new,lambda.new) %*% JackInv
		names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",paste("beta",1:(nmark-1),sep=""),"lambda")
	} else {
		Var <- NULL
	}

	list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}

### the original version of the AUG estimating equation with the augmented term
### based on the conditional cdf of the mark
densRatioAUGorg <- function(mark, trt.id, aux, f.aux.misspecified=FALSE){
	V <- cbind(1,mark)
	V.complete <- na.omit(V)
	z <- trt.id
	na.idx <- attr(V.complete,"na.action")
	if (!is.null(na.idx)){
		z.complete <- z[-na.idx]
		aux.complete <- aux[-na.idx]
	} else {
		z.complete <- z
		aux.complete <- aux
	}
	nmark <- NCOL(V.complete)
	ninf <- NROW(V.complete)
	VV.complete <- apply(V.complete,1,tcrossprod)
		
	g <- function(theta){ exp(drop(V.complete %*% theta)) }
	dG <- function(theta){ t(g(theta) * V.complete) }
	d2G <- function(theta){ array(t(t(VV.complete)*g(theta)),dim=c(nmark,nmark,ninf)) }
	dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

	score1 <- function(theta, lambda, pi, aux1, aux2){
		drop(-lambda * dG(theta) %*% ((1/(1+lambda*(g(theta)-1)))*(1/pi + aux1)) + dG(theta) %*% ((1/g(theta))*(z.complete/pi + aux2)) )
	}
	score2 <- function(theta, lambda, pi, aux1){
		-sum((g(theta)-1)*(1/(1+lambda*(g(theta)-1)))*(1/pi + aux1) )
	}
	score <- function(theta, lambda, pi, aux1, aux2){ c(score1(theta,lambda,pi,aux1,aux2),score2(theta,lambda,pi,aux1)) }
	jack11 <- function(theta, lambda, pi, aux1, aux2){
		d2Gperm <- aperm(d2G(theta), c(3,1,2))
		dGdGperm <- aperm(dGdG(theta), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(1+lambda*(g(theta)-1)))*(1/pi + aux1), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(1+lambda*(g(theta)-1))^2)*(1/pi + aux1), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(1/g(theta))*(z.complete/pi + aux2), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(1/g(theta)^2)*(z.complete/pi + aux2), c(2,3,1)),c(1,2),sum)
		-lambda*(term1 - lambda*term2) + term3 - term4
	}
	jack21 <- function(theta, lambda, pi, aux1){
		drop(-dG(theta) %*% ((1/(1+lambda*(g(theta)-1))^2)*(1/pi + aux1)) )
	}
	jack22 <- function(theta, lambda, pi, aux1){
		sum(((g(theta)-1)^2)*(1/(1+lambda*(g(theta)-1))^2)*(1/pi + aux1) )
	}
	jack <- function(theta, lambda, pi, aux1, aux2){
		j21 <- jack21(theta,lambda,pi,aux1)
		cbind(rbind(jack11(theta,lambda,pi,aux1,aux2),j21),c(j21,jack22(theta,lambda,pi,aux1)))
	}

	r <- apply(V, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
	pi.all <- glm(r ~ z + aux + z*aux, family=binomial)$fitted
	if (!is.null(na.idx)){
		pi <- pi.all[-na.idx]
	} else {
		pi <- pi.all
	}
	if (any(pi==0)){ stop("Selection probabilities not bounded away from 0.") }

	coef.ipw <- densRatioIPW(mark, trt.id, aux)$coef
	theta.ipw <- coef.ipw[-(nmark+1)]
	lambda.ipw <- coef.ipw[nmark+1]
	# f.mark <- 1/(sum(r/pi.all)*(1 + lambda.ipw*(g(theta.ipw) - 1)))
	f.mark <- 1/(ninf*(1 + lambda.ipw*(g(theta.ipw) - 1)))
	kappa.mle <- max(pmax(V.complete[,2]/aux.complete, (1-V.complete[,2])/(1-aux.complete))) - 1

	# 'mark.vect' has the length of 'V.complete'
	f.aux <- function(mark.vect, aux.val){
		((1+kappa.mle)/kappa.mle)*ifelse((mark.vect/(1+kappa.mle) <= aux.val) & (aux.val <= (mark.vect+kappa.mle)/(1+kappa.mle)),1,0)
	}
	r.hat <- function(mark.vect, z.val, aux.val){
		r.hat.num <- f.aux(mark.vect, aux.val) * (z.val*g(theta.ipw) + (1-z.val)) * f.mark
		return( r.hat.num/sum(r.hat.num) )
	}
	r.aux <- sapply(1:ninf, function(i){ r.hat(V.complete[,2], z.complete[i], aux.complete[i]) })
	weighted.r.aux.sums <- drop(r.aux %*% (1-1/pi))
	z.weighted.r.aux.sums <- drop(r.aux %*% (z.complete * (1-1/pi)))

#	if (f.aux.misspecified){
#		f.aux <- ((1+kappa.mle)/kappa.mle)*ifelse((V.complete[,2]/(1+kappa.mle) <= aux.complete) & (aux.complete <= (V.complete[,2]+kappa.mle)/(1+kappa.mle)),1,0)
#	} else {
#		f.aux <- ((1+kappa.mle)/kappa.mle)*rep(1,ninf)
#	}
#	r.hat.num <- f.aux * (z.complete*g(theta.ipw) + (1-z.complete)) * f.mark
#	r.hat <- r.hat.num/sum(r.hat.num)

	param.old <- numeric(nmark+1)
	param.new <- c(numeric(nmark),0.5)
	while (sum((param.new - param.old)^2)>1e-8){
		param.old <- param.new
		jackInv <- try(solve(jack(param.old[-(nmark+1)],param.old[nmark+1],pi=pi,aux1=weighted.r.aux.sums,aux2=z.weighted.r.aux.sums)), silent=TRUE)
		if (class(jackInv)!="try-error"){
			param.new <- param.old - drop(jackInv %*% score(param.old[-(nmark+1)],param.old[nmark+1],pi=pi,aux1=weighted.r.aux.sums,aux2=z.weighted.r.aux.sums))
		}
		if (sum(is.nan(param.new))>0) break
	}
	theta.new <- param.new[-(nmark+1)]
	lambda.new <- param.new[nmark+1]
	
	Resid <- function(theta, lambda, pi, pi.all, r.aux){
		U <- matrix(0,nrow=length(z),ncol=nmark+1)
		U[-na.idx,1:nmark] <- -lambda * t(dG(theta)) * (1/(1+lambda*(g(theta)-1)))*(1/pi) + t(-lambda * dG(theta) %*% ((1/(1+lambda*(g(theta)-1))) * t(t(r.aux) * (1-1/pi)))) + t(dG(theta)) * (z.complete/g(theta))*(1/pi) + t(dG(theta) %*% ((1/g(theta)) * t(t(r.aux) * z.complete *(1-1/pi)))) 
		U[-na.idx,nmark+1] <- (g(theta)-1)*(1/(1+lambda*(g(theta)-1)))*(1/pi) + drop((t(r.aux) * (1-1/pi)) %*% ((g(theta)-1)/(1+lambda*(g(theta)-1))))
		S <- (r-pi.all) * cbind(1,z,aux,z*aux)
		resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
		Resids <- do.call("cbind",resids)
		crossprod(Resids)/ninf
	}

	JackInv <- try(solve(jack(theta.new,lambda.new,pi=pi,aux1=weighted.r.aux.sums,aux2=z.weighted.r.aux.sums)), silent=TRUE)
	if (class(JackInv)!="try-error"){
		Var <- ninf * JackInv %*% Resid(theta.new,lambda.new,pi=pi,pi.all=pi.all,r.aux=r.aux) %*% JackInv
		names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",paste("beta",1:(nmark-1),sep=""),"lambda")
	} else {
		Var <- NULL
	}

	list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new,pi,aux1=weighted.r.aux.sums,aux2=z.weighted.r.aux.sums), probs=pi, conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}


