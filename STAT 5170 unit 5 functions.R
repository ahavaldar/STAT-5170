# =================================================================================
# STAT 5170: Applied Time Series
# Specialized functons for unit 5
# =================================================================================
#
# The four functions defined below are described as follows.

# run.dsamp.reg.wn: This function is used for de-trending a time series by regression, under the assumption that the residuals are a white-noise time series. The core algorithm used in this function is essentially a copy of the algorithm we explored in the previous learning unit for regression analysis. The sample from the posterior distribution is simulated by direct sampling within a Gibbs strategy. The algorithm itself executes quite quickly.

# run.mcmc.reg.ar1: This function is used for de-trending a time series by regression, under the assumption that the residuals are an AR(1) time series. For a given correlation matrix R, the posterior distribution reflects that obtained under Bayesian regression as usual, except that the least squares statistics that define its mean vector and covariance matrix (of the regression-coefficient posterior distribution), and sum-of-squared error statisitcs (of the white-noise variance posterior distribution) incorporate the correlation matrix in the same manner as in the modified regression summary function defined above. The model does not assume a fixed starting value of the time-series of residual deviations, which distinguishes it from the AR model and computational algorithm we examined in the previous learning unit. The sample from the posterior distribution is simulated using a Gibbs strategy within an MCMC algorithm, using direct sampling to simulate the regression coefficients and white-noise variance, and an Metropolis-Hastings step to simulate the autoregressive parameter. The algorithm requires a burn in phase and monitoring in order to specify its tuning parameters.

# run.mcmc.diff.wn: This function is used for differencing analysis of a time series, under the assumption that the time series of differences is a white-noise time series with a non-zero mean. The model assumes fixed initial value x0 of the original time series, before differencing, which is treated as an additional parameter. The sample from the posterior distribution is simulated using a Gibbs strategy within an MCMC algorithm, using direct sampling to simulate the regression coefficients and white-noise variance, and an Metropolis-Hastings step to simulate the initial value x0. The algorithm requires a burn in phase and monitoring in order to specify its tuning parameters.

# run.mcmc.diff.ar1: This function is used for differencing analysis of a time series, under the assumption that the time series of differences is an AR(1) time series, possibly with a non-zero mean. The posterior distribution of the time series of differences matches that of an intercept-only regression analysis (with no regressor variables) in order to model the non-zero mean. The model does not assume a fixed starting value of the difference time series. However, it does assume a fixed inital value x0 of the original time series, before differencing. The sample from the posterior distribution is simulated using a Gibbs strategy within an MCMC algorithm, using direct sampling to simulate the regression coefficients and white-noise variance, and an Metropolis-Hastings step to simulate the initial value x0 and autoregressive parameter. The algorithm requires a burn in phase and monitoring in order to specify its tuning parameters.

# All of these functions produce similated samples of individual parameters, of the log-data-generating density evaluated at those paramters, of measured and predicted test statistics used to check for autocorrelation in relevant components of the model. The functons also calculate the posterior predictive p-values associated with the predicted autocorrelations, and the predictive diagnostics DIC and WAIC. The functions "run.mcmc.diff.wn" and "run.mcmc.diff.ar1" also produce samples of m-step-ahead forecasts.

# The input to every one of these functions consists of two lists. In what follows, the name of these lists are taken to be "model.parms" and "comp.parms". Output is returned in the form of a list. In what follows, the name of this list is taken to be "result". Descriptions of the variables they contain are as follows:

# ***** run.dsamp.reg.wn *****

# model.parms$x.vect = the time series to be analyzed stored as a response vector
# model.parms$Z.mat = the associated matrix of regressors

# comp.parms$n.samp = the number of iterations the agorithm is to be run
# comp.parms$max.lag = the maximum lag autotocorrelation at which predicted autocorrelations are calcuated

# result$samp$gamma0: simulated sample of the lag-0 autocorrelation
# result$samp$beta: simulated sample of the regression coefficients
# result$samp$logdgdens: simulated sample of the log-data-generating density evaluated at the model parameters

# result$stat$DIC: Deviance information statistic
# result$stat$p.DIC: effective number of parameters assocaited with the DIC statistic
# result$stat$WAIC: Watanabe Akiake information statistic
# result$stat$p.WAIC: effective number of parameters assocaited with the WAIC statistic

# result$pval$indiv: posterior predictive p-values associated with individual absolute-values of autocorrelations in the predicted residual deviations.
# result$pval$max: posterior predictive p-values associated with the maximum of absolute-value autocorrelations up to a lag
# result$pval$sum: posterior predictive p-values associated with the cumulative sum of absolute-value autocorrelations up to a lag

# ***** run.mcmc.reg.ar1 *****

# Input variables are as in "run.dsamp.reg.wn" with the addition of the following

# model.parms$phi1 = an initial value for the autoregression parameter.
# comp.parms$phi1.rad = the radius used to simulate proposed values of the autoregression parameter using a uniform shift with reflecting boundaries 
# comp.parms$max.kappa = a maximum allowed "condition number" used to check the stability of the autocorrelation matrix. Set this value to 1000.

# The variable "comp.parms$max.lag" doubles in specifying the maximum lag of bands that are specified in the sparse autocorrelation matrix, R

# The output variables are as in "run.dsamp.reg.wn" with the addition of the following

# result$rate$phi1: acceptance rate of proposed values of the autoregressive parameter
# result$rate$poorcond: the rate at which a proposed candidate value produced a poorly conditioned autocorrelation matrix, hence was rejected.

# ***** run.mcmc.diff.wn *****

# model.parms$x = the time series to be analyzed stored as a sequence of numeric values
# model.parms$x0 = an initial value for time series's initial value, x0

# comp.parms$n.samp = the number of iterations the agorithm is to be run
# comp.parms$max.lag = the maximum lag autotocorrelation at which predicted autocorrelations are calcuated
# comp.parms$x0.rad = the radius used to simulate proposed values of the time series's initial value using a uniform shift
# comp.parms$m = the maximum number steps ahead in m-step-ahead predictions

# result$samp$gamma0: simulated sample of the lag-0 autocorrelation
# result$samp$x0: simulated sample of the time series's initial value
# result$samp$delta: simulated sample of the mean of the time series of differences
# result$samp$logdgdens: simulated sample of the log-data-generating density evaluated at the model parameters
# result$samp$fore: simulated sample of the predicted m-step ahead forecasts

# result$stat$DIC: Deviance information statistic
# result$stat$p.DIC: effective number of parameters assocaited with the DIC statistic
# result$stat$WAIC: Watanabe Akiake information statistic
# result$stat$p.WAIC: effective number of parameters assocaited with the WAIC statistic

# result$pval$indiv: posterior predictive p-values associated with individual absolute-values of autocorrelations in the predicted time series of differences
# result$pval$max: posterior predictive p-values associated with the maximum of absolute-value autocorrelations up to a lag
# result$pval$sum: posterior predictive p-values associated with the cumulative sum of absolute-value autocorrelations up to a lag

# result$rate$x0: acceptance rate of proposed values of the time series's initial value

# ***** run.mcmc.diff.ar1 *****

# Input variables are as in "run.mcmc.diff.wn" with the addition of the following

# model.parms$phi1 = an initial value for the autoregression parameter.
# comp.parms$phi1.rad = the radius used to simulate proposed values of the autoregression parameter using a uniform shift with reflecting boundaries 
# comp.parms$max.kappa = a maximum allowed "condition number" used to check the stability of the autocorrelation matrix. Set this value to 1000.

# The variable "comp.parms$max.lag" doubles in specifying the maximum lag of bands that are specified in the sparse autocorrelation matrix, R

# The output variables are as in "run.mcmc.diff.wn" with the addition of the following

# result$rate$phi1: acceptance rate of proposed values of the autoregressive parameter
# result$rate$poorcond: the rate at which a proposed candidate value produced a poorly conditioned autocorrelation matrix, hence was rejected.

# =================================================================================
# Bayesian regression analyses for detrending
# =================================================================================
#
# ---------------------------------------------------------------------------------
# Residuals are a white noise time series
# ---------------------------------------------------------------------------------

run.dsamp.reg.wn <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	reg.summary <- function(model.parms) {
		n <- dim(model.parms$Z.mat)[1]
		r <- dim(model.parms$Z.mat)[2]
		ZpZ <- t(model.parms$Z.mat) %*% model.parms$Z.mat
		ZpX <- t(model.parms$Z.mat) %*% model.parms$x.vect
		ZpZ.inv <- solve(ZpZ)
		sqrt.ZpZ.inv <- sqrtm(ZpZ.inv)
		b.hat <- ZpZ.inv %*% ZpX
		res <- model.parms$x.vect - model.parms$Z.mat %*% b.hat
		SSE <- sum(res^2)
		result <- list(ZpZ.inv=ZpZ.inv, sqrt.ZpZ.inv=sqrt.ZpZ.inv, b.hat=b.hat, res=res, SSE=SSE)
		return(result)
	}
	get.logdgendens <- function(model.parms, ls.stat) {
		n <- dim(model.parms$Z.mat)[1]
		result <- list()
		result$pts <- -0.5*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$res^2/model.parms$gamma0
		result$all <- -0.5*n*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$SSE/model.parms$gamma0
		return(result)
	}
	sim.regen.ts <- function(model.parms) {
		n <- dim(model.parms$Z.mat)[1]
		regen.model.parms <- model.parms
		regen.model.parms$res <- rnorm(n=n, mean=0, sd=sqrt(model.parms$gamma0))
		return(regen.model.parms)
	}
	pred.acf <- function(model.parms, comp.parms) {
		n <- dim(model.parms$Z.mat)[1]
		acf.seq <- numeric(length=comp.parms$max.lag)
		for (lag in 1:comp.parms$max.lag) {
			acf.seq[lag] <- sum(model.parms$res[(1+lag):n]*model.parms$res[1:(n-lag)])/((n-lag)*model.parms$gamma0)
		}
		return(acf.seq)
	}
	# ----- ----- --- main --- ----- -----
	n <- dim(model.parms$Z.mat)[1]
	r <- dim(model.parms$Z.mat)[2]
	nu <- n - r
	# Initial iteration
	curr.model.parms <- model.parms
	curr.ls.stat <- reg.summary(curr.model.parms)
	# Set up space for storing results
	post.gamma0.samp <- numeric(length=comp.parms$n.samp)
	post.beta.samp <- matrix(nrow=r, ncol=comp.parms$n.samp)
	logdgdens.samp <- list()
	logdgdens.samp$pts <- matrix(nrow=n, ncol=comp.parms$n.samp)
	logdgdens.samp$all <- numeric(length=comp.parms$n.samp)
	acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	# Implement remaining simulation
	for (i.samp in 1:comp.parms$n.samp) {
		# Update gamma0 and beta
		curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
		curr.model.parms$beta <- curr.ls.stat$b.hat + curr.ls.stat$sqrt.ZpZ.inv %*% as.matrix(rnorm(n=r, mean=0, sd=sqrt(curr.model.parms$gamma0)))
		curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
		# Record sample values
		post.gamma0.samp[i.samp] <- curr.model.parms$gamma0
		post.beta.samp[, i.samp] <- curr.model.parms$beta
		logdgdens.samp$pts[, i.samp] <- curr.model.parms$logpdens$pts
		logdgdens.samp$all[i.samp] <- curr.model.parms$logpdens$all
		# Regenerate the time series, summarize, and record summaries
		curr.model.parms$res <- model.parms$x.vect - model.parms$Z.mat %*% curr.model.parms$beta
		acf.measr.samp[, i.samp] <- pred.acf(curr.model.parms, comp.parms)
		regen.model.parms <- sim.regen.ts(curr.model.parms)
		acf.regen.samp[, i.samp] <- pred.acf(regen.model.parms, comp.parms)
		for (lag in 1:comp.parms$max.lag) {
			max.acf.measr.samp[lag, i.samp] <- max(abs(acf.measr.samp[1:lag, i.samp]))
			max.acf.regen.samp[lag, i.samp] <- max(abs(acf.regen.samp[1:lag, i.samp]))
			sum.acf.measr.samp[lag, i.samp] <- sum(abs(acf.measr.samp[1:lag, i.samp]))
			sum.acf.regen.samp[lag, i.samp] <- sum(abs(acf.regen.samp[1:lag, i.samp]))
		}
	}
	# Calculate DIC
	Bayes.model.parms <- model.parms
	Bayes.model.parms$gamma0 <- exp(mean(log(post.gamma0.samp)))
	Bayes.model.parms$beta <- rowSums(post.beta.samp) / comp.parms$n.samp
	Bayes.ls.stat <- reg.summary(Bayes.model.parms)
	Bayes.model.parms$logpdens <- get.logdgendens(Bayes.model.parms, Bayes.ls.stat)
	p.DIC <- 2*(Bayes.model.parms$logpdens$all - mean(logdgdens.samp$all))
	DIC <- -2*Bayes.model.parms$logpdens$all + 2*p.DIC
	# Calculate WAIC
	mean.dgdens <- rowSums(exp(logdgdens.samp$pts)) / comp.parms$n.samp
	sum.log.mean.dgdens <- sum(log(mean.dgdens))
	mean.logdgdens <- rowSums(logdgdens.samp$pts) / comp.parms$n.samp
	p.WAIC <- 2*(sum.log.mean.dgdens - sum(mean.logdgdens))
	WAIC <- -2*sum.log.mean.dgdens + 2*p.WAIC
	# Calculate posterior predictive p-values
	pval.indiv.acf <- numeric(length=comp.parms$max.lag)
	pval.max.acf <- numeric(length=comp.parms$max.lag)
	pval.sum.acf <- numeric(length=comp.parms$max.lag)
	for (lag in 1:comp.parms$max.lag) {
		pval.indiv.acf[lag] <- length(which(abs(acf.regen.samp[lag,]) > abs(acf.measr.samp[lag,]))) / comp.parms$n.samp
		pval.max.acf[lag] <- length(which(max.acf.regen.samp[lag,] > max.acf.measr.samp[lag,])) / comp.parms$n.samp
		pval.sum.acf[lag] <- length(which(sum.acf.regen.samp[lag,] > sum.acf.measr.samp[lag,])) / comp.parms$n.samp
	}
	# Record results
	samp <- list(gamma0=post.gamma0.samp, beta=post.beta.samp, logdgdens=logdgdens.samp)
	stat <- list(DIC=DIC, p.DIC=p.DIC, WAIC=WAIC, p.WAIC=p.WAIC)
	pval <- list(indiv=pval.indiv.acf, max=pval.max.acf, sum=pval.sum.acf)
	result <- list(samp=samp, stat=stat, pval=pval)
	return(result)
}

# ---------------------------------------------------------------------------------
# Residuals are an AR(1) time series
# ---------------------------------------------------------------------------------

run.mcmc.reg.ar1 <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	get.acf.ar1 <- function(phi1, max.lag) {
		acf.seq <- phi1^(0:max.lag)
		return(acf.seq)
	}
	reg.summary.acorr <- function(model.parms) {
		n <- dim(model.parms$Z.mat)[1]
		r <- dim(model.parms$Z.mat)[2]
		max.lagp1 <- length(model.parms$acf.seq)
		max.lag <- max.lagp1-1
		bands <- t(matrix(data=model.parms$acf.seq[1:max.lagp1], nrow= max.lagp1, ncol=n))
		R.mat <- bandSparse(n=n, k=0:max.lag, diag=bands, symm=TRUE) # The bandSparse() function is specialized syntax for banded matrices
		kappa.val <- kappa(R.mat) # Check for poor conditioning
		detR <- det(R.mat)
		if (detR <= 0) {
			kappa.val <- Inf
		}
		logdetR <- log(abs(detR))
		R.inv <- solve(R.mat)
		ZpRZ <- as.matrix(t(model.parms$Z.mat) %*% R.inv %*% model.parms$Z.mat)
		ZpRX <- as.matrix(t(model.parms$Z.mat) %*% R.inv %*% model.parms$x.vect)
		ZpRZ.inv <- solve(ZpRZ)
		sqrt.ZpRZ.inv <- sqrtm(ZpRZ.inv)
		b.hat <- ZpRZ.inv %*% ZpRX
		res <- model.parms$x.vect - model.parms$Z.mat %*% b.hat
		SSE <- as.numeric(t(res) %*% R.inv %*% res)
		if (SSE < 0) {
			kappa.val <- Inf
		}
		result <- list(R.inv=R.inv, ZpRZ.inv=ZpRZ.inv, sqrt.ZpRZ.inv=sqrt.ZpRZ.inv, b.hat=b.hat, res=res, SSE=SSE, logdetR=logdetR, R.inv=R.inv, kappa=kappa.val)
		return(result)
	}
	get.logdgendens <- function(model.parms, ls.stat) {
		n <- dim(model.parms$Z.mat)[1]
		result <- list()
		result$pts <- -0.5*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$res^2/model.parms$gamma0
		result$all <- -0.5*n*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$logdetR - 0.5*ls.stat$SSE/model.parms$gamma0
		return(result)
	}
	propose.update.phi1 <- function(curr.phi1, radius) {
		result <- NA
		cand.phi1 <- curr.phi1 + runif(n=1, min=-radius, max=radius)
		hi <- 1
		lo <- -1
		if (cand.phi1 < lo) {
			cand.phi1 <- lo + (lo - cand.phi1)
		} else if (cand.phi1 > hi) {
			cand.phi1 <- hi - (cand.phi1 - hi)
		}
		if ((cand.phi1 > lo) & (cand.phi1 < hi)) {
			result <- cand.phi1
		}
		return(result)
	}
	sim.regen.ts <- function(model.parms) {
		n <- dim(model.parms$Z.mat)[1]
		ar1.samp <- numeric(length=n)
		ar1.samp[1] <- rnorm(n=1, mean=0, sd=sqrt(model.parms$gamma0))
		sigw <- sqrt((1-model.parms$phi1^2)*model.parms$gamma0)
		wn.samp <- rnorm(n=n, mean=0, sd=sigw)
		pastx <- ar1.samp[1]
		for (t in 2:n) {
			ar1.samp[t] <- sum(model.parms$phi1*pastx) + wn.samp[t]
			pastx <- ar1.samp[t]
		}
		regen.model.parms <- model.parms
		regen.model.parms$res <- ar1.samp
		return(regen.model.parms)
	}
	pred.acf <- function(model.parms, comp.parms) {
		n <- dim(model.parms$Z.mat)[1]
		acf.seq <- numeric(length=comp.parms$max.lag)
		for (lag in 1:comp.parms$max.lag) {
			acf.seq[lag] <- sum(model.parms$res[(1+lag):n]*model.parms$res[1:(n-lag)])/((n-lag)*model.parms$gamma0)
		}
		return(acf.seq)
	}
	# ----- ----- --- main --- ----- -----
	n <- dim(model.parms$Z.mat)[1]
	r <- dim(model.parms$Z.mat)[2]
	nu <- n - r
	# Initial iteration
	curr.model.parms <- model.parms
	curr.model.parms$acf.seq <- get.acf.ar1(model.parms$phi1, comp.parms$max.lag)
	curr.ls.stat <- reg.summary.acorr(curr.model.parms)
	curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
	curr.model.parms$beta <- curr.ls.stat$b.hat + curr.ls.stat$sqrt.ZpRZ.inv %*% as.matrix(rnorm(n=r, mean=0, sd=sqrt(curr.model.parms$gamma0)))
	curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
	cand.model.parms <- curr.model.parms
	# Set up space for storing results
	post.gamma0.samp <- numeric(length=comp.parms$n.samp)
	post.beta.samp <- matrix(nrow=r, ncol=comp.parms$n.samp)
	post.phi1.samp <- numeric(length=comp.parms$n.samp)
	accept.phi1.cnt <- 0
	poorcond.cnt <- 0
	logdgdens.samp <- list()
	logdgdens.samp$pts <- matrix(nrow=n, ncol=comp.parms$n.samp)
	logdgdens.samp$all <- numeric(length=comp.parms$n.samp)
	acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	# Implement remaining simulation
	for (i.samp in 1:comp.parms$n.samp) {
		# Update phi1
		cand.model.parms$phi1 <- propose.update.phi1(curr.model.parms$phi1, comp.parms$phi1.rad)
		cand.model.parms$acf.seq <- get.acf.ar1(cand.model.parms$phi1, comp.parms$max.lag)
		cand.ls.stat <- reg.summary.acorr(cand.model.parms)
		if (cand.ls.stat$kappa < comp.parms$max.kappa) {
			cand.model.parms$logpdens <- get.logdgendens(cand.model.parms, cand.ls.stat)
			log.accept.prob <- cand.model.parms$logpdens$all - curr.model.parms$logpdens$all
		} else {
			poorcond.cnt <- poorcond.cnt + 1
			log.accept.prob <- -Inf
		}
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.phi1.cnt <- accept.phi1.cnt + 1
			curr.model.parms <- cand.model.parms
			curr.ls.stat <- cand.ls.stat
		} else {
			cand.model.parms <- curr.model.parms
			cand.ls.stat <- curr.ls.stat
		}
		# Update gamma0 and beta
		curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
		curr.model.parms$beta <- curr.ls.stat$b.hat + curr.ls.stat$sqrt.ZpRZ.inv %*% as.matrix(rnorm(n=r, mean=0, sd=sqrt(curr.model.parms$gamma0)))
		curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
		# Record sample values
		post.phi1.samp[i.samp] <- curr.model.parms$phi1
		post.gamma0.samp[i.samp] <- curr.model.parms$gamma0
		post.beta.samp[, i.samp] <- curr.model.parms$beta
		logdgdens.samp$pts[, i.samp] <- curr.model.parms$logpdens$pts
		logdgdens.samp$all[i.samp] <- curr.model.parms$logpdens$all
		# Regenerate the time series, summarize, and record summaries
		curr.model.parms$res <- model.parms$x.vect - model.parms$Z.mat %*% curr.model.parms$beta
		acf.measr.samp[, i.samp] <- pred.acf(curr.model.parms, comp.parms)
		regen.model.parms <- sim.regen.ts(curr.model.parms)
		acf.regen.samp[, i.samp] <- pred.acf(regen.model.parms, comp.parms)
		for (lag in 1:comp.parms$max.lag) {
			max.acf.measr.samp[lag, i.samp] <- max(abs(acf.measr.samp[1:lag, i.samp]))
			max.acf.regen.samp[lag, i.samp] <- max(abs(acf.regen.samp[1:lag, i.samp]))
			sum.acf.measr.samp[lag, i.samp] <- sum(abs(acf.measr.samp[1:lag, i.samp]))
			sum.acf.regen.samp[lag, i.samp] <- sum(abs(acf.regen.samp[1:lag, i.samp]))
		}
	}
	# Calculate DIC
	Bayes.model.parms <- model.parms
	Bayes.model.parms$phi1 <- mean(post.phi1.samp)
	Bayes.model.parms$gamma0 <- exp(mean(log(post.gamma0.samp)))
	Bayes.model.parms$beta <- rowSums(post.beta.samp) / comp.parms$n.samp
	Bayes.model.parms$acf.seq <- get.acf.ar1(Bayes.model.parms$phi1, comp.parms$max.lag)
	Bayes.ls.stat <- reg.summary.acorr(Bayes.model.parms)
	Bayes.model.parms$logpdens <- get.logdgendens(Bayes.model.parms, Bayes.ls.stat)
	p.DIC <- 2*(Bayes.model.parms$logpdens$all - mean(logdgdens.samp$all))
	DIC <- -2*Bayes.model.parms$logpdens$all + 2*p.DIC
	# Calculate WAIC
	mean.dgdens <- rowSums(exp(logdgdens.samp$pts)) / comp.parms$n.samp
	sum.log.mean.dgdens <- sum(log(mean.dgdens))
	mean.logdgdens <- rowSums(logdgdens.samp$pts) / comp.parms$n.samp
	p.WAIC <- 2*(sum.log.mean.dgdens - sum(mean.logdgdens))
	WAIC <- -2*sum.log.mean.dgdens + 2*p.WAIC
	# Calculate posterior predictive p-values
	pval.indiv.acf <- numeric(length=comp.parms$max.lag)
	pval.max.acf <- numeric(length=comp.parms$max.lag)
	pval.sum.acf <- numeric(length=comp.parms$max.lag)
	for (lag in 1:comp.parms$max.lag) {
		pval.indiv.acf[lag] <- length(which(abs(acf.regen.samp[lag,]) > abs(acf.measr.samp[lag,]))) / comp.parms$n.samp
		pval.max.acf[lag] <- length(which(max.acf.regen.samp[lag,] > max.acf.measr.samp[lag,])) / comp.parms$n.samp
		pval.sum.acf[lag] <- length(which(sum.acf.regen.samp[lag,] > sum.acf.measr.samp[lag,])) / comp.parms$n.samp
	}
	# Record results
	samp <- list(gamma0=post.gamma0.samp, beta=post.beta.samp, phi1=post.phi1.samp, logdgdens=logdgdens.samp)
	rate <- list(phi1=accept.phi1.cnt/comp.parms$n.samp, poorcond=poorcond.cnt/comp.parms$n.samp)
	stat <- list(DIC=DIC, p.DIC=p.DIC, WAIC=WAIC, p.WAIC=p.WAIC)
	pval <- list(indiv=pval.indiv.acf, max=pval.max.acf, sum=pval.sum.acf)
	result <- list(samp=samp, rate=rate, stat=stat, pval=pval)
	return(result)
}

# =================================================================================
# Bayesian differencing analyses
# =================================================================================

# ---------------------------------------------------------------------------------
# Differences are a shifted white noise time series
# ---------------------------------------------------------------------------------

run.mcmc.diff.wn <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	reg.summary <- function(model.parms) {
		n <- length(model.parms$x)
		sqrt.v.inv <- 1/sqrt(n)
		d.hat <- mean(model.parms$x.diff)
		res <- model.parms$x.diff - d.hat
		SSE <- sum(res^2)
		result <- list(sqrt.v.inv=sqrt.v.inv, d.hat=d.hat, res=res, SSE=SSE)
		return(result)
	}
	get.logdgendens <- function(model.parms, ls.stat) {
		n <- length(model.parms$x)
		result <- list()
		result$pts <- -0.5*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$res^2/model.parms$gamma0
		result$all <- -0.5*n*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$SSE/model.parms$gamma0
		return(result)
	}
	propose.update.x0 <- function(curr.x0, radius) {
		cand.x0 <- curr.x0 + runif(n=1, min=-radius, max=radius)
		return(cand.x0)
	}
	sim.regen.ts <- function(model.parms) {
		n <- length(model.parms$x)
		regen.model.parms <- model.parms
		regen.model.parms$x.diff <- rnorm(n=n, mean=model.parms$delta, sd=sqrt(model.parms$gamma0))
		pastx <- model.parms$x0
		for (t in 1:n) {
			regen.model.parms$x[t] <- regen.model.parms$x.diff[t] + pastx
			pastx <- regen.model.parms$x[t]
		}
		return(regen.model.parms)
	}
	sim.future.ts <- function(model.parms, regen.model.parms, comp.parms) {
		n <- length(model.parms$x)
		fore.model.parms <- model.parms
		fore.model.parms$x.diff <- rnorm(n=comp.parms$m, mean=model.parms$delta, sd=sqrt(model.parms$gamma0))
		fore.model.parms$x <- numeric(length=comp.parms$m)
		pastx <- model.parms$x[n]
		for (t in 1:comp.parms$m) {
			fore.model.parms$x[t] <- fore.model.parms$x.diff[t] + pastx
			pastx <- fore.model.parms$x[t]
		}
		return(fore.model.parms)
	}
	pred.acf <- function(model.parms, comp.parms) {
		n <- length(model.parms$x)
		acf.seq <- numeric(length=comp.parms$max.lag)
		for (lag in 1:comp.parms$max.lag) {
			acf.seq[lag] <- sum(model.parms$x.diff[(1+lag):n]*model.parms$x.diff[1:(n-lag)])/((n-lag)*model.parms$gamma0)
		}
		return(acf.seq)
	}
	# ----- ----- --- main --- ----- -----
	n <- length(model.parms$x)
	nu <- n - 1
	# Initial iteration
	curr.model.parms <- model.parms
	curr.model.parms$x.diff <- c(model.parms$x[1] - model.parms$x0, diff(model.parms$x))
	curr.ls.stat <- reg.summary(curr.model.parms)
	curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
	curr.model.parms$delta <- curr.ls.stat$d.hat + curr.ls.stat$sqrt.v.inv*rnorm(n=1, mean=0, sd=sqrt(curr.model.parms$gamma0))
	curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
	cand.model.parms <- curr.model.parms
	# Set up space for storing results
	post.gamma0.samp <- numeric(length=comp.parms$n.samp)
	post.delta.samp <- numeric(length=comp.parms$n.samp)
	post.x0.samp <- numeric(length=comp.parms$n.samp)
	accept.x0.cnt <- 0
	logdgdens.samp <- list()
	logdgdens.samp$pts <- matrix(nrow=n, ncol=comp.parms$n.samp)
	logdgdens.samp$all <- numeric(length=comp.parms$n.samp)
	acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	post.fore.samp <- matrix(nrow=comp.parms$m, ncol=comp.parms$n.samp)
	# Implement remaining simulation
	for (i.samp in 1:comp.parms$n.samp) {
		# Update x0
		cand.model.parms$x0 <- propose.update.x0(curr.model.parms$x0, comp.parms$x0.rad)
		cand.model.parms$x.diff[1] <- model.parms$x[1] - cand.model.parms$x0
		cand.ls.stat <- reg.summary(cand.model.parms)
		cand.model.parms$logpdens <- get.logdgendens(cand.model.parms, cand.ls.stat)
		log.accept.prob <- cand.model.parms$logpdens$all - curr.model.parms$logpdens$all
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.x0.cnt <- accept.x0.cnt + 1
			curr.model.parms <- cand.model.parms
			curr.ls.stat <- cand.ls.stat
		} else {
			cand.model.parms <- curr.model.parms
			cand.ls.stat <- curr.ls.stat
		}
		# Update gamma0 and beta
		curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
		curr.model.parms$delta <- curr.ls.stat$d.hat + curr.ls.stat$sqrt.v.inv*rnorm(n=1, mean=0, sd=sqrt(curr.model.parms$gamma0))
		curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
		# Record sample values
		post.x0.samp[i.samp] <- curr.model.parms$x0
		post.gamma0.samp[i.samp] <- curr.model.parms$gamma0
		post.delta.samp[i.samp] <- curr.model.parms$delta
		logdgdens.samp$pts[, i.samp] <- curr.model.parms$logpdens$pts
		logdgdens.samp$all[i.samp] <- curr.model.parms$logpdens$all
		# Regenerate the time series, summarize, and record summaries
		acf.measr.samp[, i.samp] <- pred.acf(curr.model.parms, comp.parms)
		regen.model.parms <- sim.regen.ts(curr.model.parms)
		fore.model.parms <- sim.future.ts(curr.model.parms, regen.model.parms, comp.parms)
		post.fore.samp[, i.samp] <- fore.model.parms$x
		acf.regen.samp[, i.samp] <- pred.acf(regen.model.parms, comp.parms)
		for (lag in 1:comp.parms$max.lag) {
			max.acf.measr.samp[lag, i.samp] <- max(abs(acf.measr.samp[1:lag, i.samp]))
			max.acf.regen.samp[lag, i.samp] <- max(abs(acf.regen.samp[1:lag, i.samp]))
			sum.acf.measr.samp[lag, i.samp] <- sum(abs(acf.measr.samp[1:lag, i.samp]))
			sum.acf.regen.samp[lag, i.samp] <- sum(abs(acf.regen.samp[1:lag, i.samp]))
		}
	}
	# Calculate DIC
	Bayes.model.parms <- model.parms
	Bayes.model.parms$x0 <- mean(post.x0.samp)
	Bayes.model.parms$x.diff <- c(model.parms$x[1] - Bayes.model.parms$x0, diff(model.parms$x))
	Bayes.model.parms$gamma0 <- exp(mean(log(post.gamma0.samp)))
	Bayes.model.parms$delta <- mean(post.delta.samp)
	Bayes.ls.stat <- reg.summary(Bayes.model.parms)
	Bayes.model.parms$logpdens <- get.logdgendens(Bayes.model.parms, Bayes.ls.stat)
	p.DIC <- 2*(Bayes.model.parms$logpdens$all - mean(logdgdens.samp$all))
	DIC <- -2*Bayes.model.parms$logpdens$all + 2*p.DIC
	# Calculate WAIC
	mean.dgdens <- rowSums(exp(logdgdens.samp$pts)) / comp.parms$n.samp
	sum.log.mean.dgdens <- sum(log(mean.dgdens))
	mean.logdgdens <- rowSums(logdgdens.samp$pts) / comp.parms$n.samp
	p.WAIC <- 2*(sum.log.mean.dgdens - sum(mean.logdgdens))
	WAIC <- -2*sum.log.mean.dgdens + 2*p.WAIC
	# Calculate posterior predictive p-values
	pval.indiv.acf <- numeric(length=comp.parms$max.lag)
	pval.max.acf <- numeric(length=comp.parms$max.lag)
	pval.sum.acf <- numeric(length=comp.parms$max.lag)
	for (lag in 1:comp.parms$max.lag) {
		pval.indiv.acf[lag] <- length(which(abs(acf.regen.samp[lag,]) > abs(acf.measr.samp[lag,]))) / comp.parms$n.samp
		pval.max.acf[lag] <- length(which(max.acf.regen.samp[lag,] > max.acf.measr.samp[lag,])) / comp.parms$n.samp
		pval.sum.acf[lag] <- length(which(sum.acf.regen.samp[lag,] > sum.acf.measr.samp[lag,])) / comp.parms$n.samp
	}
	# Record results
	samp <- list(x0=post.x0.samp, gamma0=post.gamma0.samp, delta=post.delta.samp, logdgdens=logdgdens.samp, fore=post.fore.samp)
	rate <- list(x0=accept.x0.cnt/comp.parms$n.samp)
	stat <- list(DIC=DIC, p.DIC=p.DIC, WAIC=WAIC, p.WAIC=p.WAIC)
	pval <- list(indiv=pval.indiv.acf, max=pval.max.acf, sum=pval.sum.acf)
	result <- list(samp=samp, rate=rate, stat=stat, pval=pval)
	return(result)
}

# ---------------------------------------------------------------------------------
# Differences are an AR(1) time series, possibly with a non-zero mean
# ---------------------------------------------------------------------------------

run.mcmc.diff.ar1 <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	get.acf.ar1 <- function(phi1, max.lag) {
		acf.seq <- phi1^(0:max.lag)
		return(acf.seq)
	}
	reg.summary.acorr <- function(model.parms) {
		n <- length(model.parms$x)
		max.lagp1 <- length(model.parms$acf.seq)
		max.lag <- max.lagp1-1
		bands <- t(matrix(data=model.parms$acf.seq[1:max.lagp1], nrow= max.lagp1, ncol=n))
		R.mat <- bandSparse(n=n, k=0:max.lag, diag=bands, symm=TRUE) # The bandSparse() function is specialized syntax for banded matrices
		kappa.val <- kappa(R.mat) # Check for poor conditioning
		detR <- det(R.mat)
		if (detR <= 0) {
			kappa.val <- Inf
		}
		logdetR <- log(abs(det(R.mat)))
		R.inv <- solve(R.mat)
		v.scl <- sum(R.inv)
		sqrt.v.inv <- 1/sqrt(abs(v.scl))
		d.hat <- sum(R.inv %*% as.matrix(model.parms$x.diff)) / v.scl
		res <- model.parms$x.diff - d.hat
		SSE <- as.numeric(t(res) %*% R.inv %*% res)
		if (SSE < 0) {
			kappa.val <- Inf
		}
		result <- list(R.inv=R.inv, v.scl=v.scl, sqrt.v.inv=sqrt.v.inv, d.hat=d.hat, res=res, SSE=SSE, logdetR=logdetR, kappa=kappa.val)
		return(result)
	}
	reg.recalc.acorr <- function(model.parms, ls.stat) {
		n <- length(model.parms$x)
		result <- ls.stat
		result$d.hat <- sum(ls.stat$R.inv %*% as.matrix(model.parms$x.diff)) / ls.stat$v.scl
		result$res <- model.parms$x.diff - result$d.hat
		result$SSE <- as.numeric(t(result$res) %*% ls.stat$R.inv %*% result$res)
		return(result)
	}
	get.logdgendens <- function(model.parms, ls.stat) {
		n <- length(model.parms$x)
		result <- list()
		result$pts <- -0.5*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$res^2/model.parms$gamma0
		result$all <- -0.5*n*log(2*pi*model.parms$gamma0) - 0.5*ls.stat$logdetR - 0.5*ls.stat$SSE/model.parms$gamma0
		return(result)
	}
	propose.update.x0 <- function(curr.x0, radius) {
		cand.x0 <- curr.x0 + runif(n=1, min=-radius, max=radius)
		return(cand.x0)
	}
	propose.update.phi1 <- function(curr.phi1, radius) {
		result <- NA
		cand.phi1 <- curr.phi1 + runif(n=1, min=-radius, max=radius[1])
		hi <- 1
		lo <- -1
		if (cand.phi1 < lo) {
			cand.phi1 <- lo + (lo - cand.phi1)
		} else if (cand.phi1 > hi) {
			cand.phi1 <- hi - (cand.phi1 - hi)
		}
		if ((cand.phi1 > lo) & (cand.phi1 < hi)) {
			result <- cand.phi1
		}
		return(result)
	}
	sim.regen.ts <- function(model.parms) {
		n <- length(model.parms$x)
		sigw <- sqrt((1-model.parms$phi1^2)*model.parms$gamma0)
		alpha.val <- (1-model.parms$phi1)*model.parms$delta
		wn.samp <- rnorm(n=n, mean=0, sd=sigw)
		regen.model.parms <- model.parms
		regen.model.parms$x.diff <- numeric(length=n)
		regen.model.parms$x.diff[1] <- rnorm(n=1, mean=0, sd=sqrt(model.parms$gamma0))
		regen.model.parms$x[1] <- regen.model.parms$x.diff[1] + model.parms$x0
		pastx.diff <- regen.model.parms$x.diff[1]
		pastx <- regen.model.parms$x[1]
		for (t in 2:n) {
			regen.model.parms$x.diff[t] <- alpha.val + model.parms$phi1*pastx.diff + wn.samp[t]
			regen.model.parms$x[t] <- regen.model.parms$x.diff[t] + pastx
			pastx.diff <- regen.model.parms$x.diff[t]
			pastx <- regen.model.parms$x[t]
		}
		return(regen.model.parms)
	}
	sim.future.ts <- function(model.parms, regen.model.parms, comp.parms) {
		n <- length(model.parms$x)
		sigw <- sqrt((1-model.parms$phi1^2)*model.parms$gamma0)
		alpha.val <- (1-model.parms$phi1)*model.parms$delta
		wn.samp <- rnorm(n=comp.parms$m, mean=0, sd=sigw)
		fore.model.parms <- model.parms
		fore.model.parms$x.diff <- numeric(length=comp.parms$m)
		fore.model.parms$x <- numeric(length=comp.parms$m)
		pastx.diff <- regen.model.parms$x.diff[n]
		pastx <- model.parms$x[n]
		for (t in 1:comp.parms$m) {
			fore.model.parms$x.diff[t] <- alpha.val + model.parms$phi1*pastx.diff + wn.samp[t]
			fore.model.parms$x[t] <- fore.model.parms$x.diff[t] + pastx
			pastx.diff <- fore.model.parms$x.diff[t]
			pastx <- fore.model.parms$x[t]
		}
		return(fore.model.parms)
	}
	pred.acf <- function(model.parms, comp.parms) {
		n <- length(model.parms$x)
		acf.seq <- numeric(length=comp.parms$max.lag)
		for (lag in 1:comp.parms$max.lag) {
			acf.seq[lag] <- sum(model.parms$x.diff[(1+lag):n]*model.parms$x.diff[1:(n-lag)])/((n-lag)*model.parms$gamma0)
		}
		return(acf.seq)
	}
	# ----- ----- --- main --- ----- -----
	n <- length(model.parms$x)
	nu <- n - 1
	# Initial iteration
	curr.model.parms <- model.parms
	curr.model.parms$x.diff <- c(model.parms$x[1] - model.parms$x0, diff(model.parms$x))
	curr.model.parms$acf.seq <- get.acf.ar1(model.parms$phi1, comp.parms$max.lag)
	curr.ls.stat <- reg.summary.acorr(curr.model.parms)
	curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
	curr.model.parms$delta <- curr.ls.stat$d.hat + curr.ls.stat$sqrt.v.inv*rnorm(n=1, mean=0, sd=sqrt(curr.model.parms$gamma0))
	curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
	cand.model.parms <- curr.model.parms
	# Set up space for storing results
	post.gamma0.samp <- numeric(length=comp.parms$n.samp)
	post.delta.samp <- numeric(length=comp.parms$n.samp)
	post.x0.samp <- numeric(length=comp.parms$n.samp)
	post.phi1.samp <- numeric(length=comp.parms$n.samp)
	accept.x0.cnt <- 0
	accept.phi1.cnt <- 0
	poorcond.cnt <- 0
	logdgdens.samp <- list()
	logdgdens.samp$pts <- matrix(nrow=n, ncol=comp.parms$n.samp)
	logdgdens.samp$all <- numeric(length=comp.parms$n.samp)
	acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	max.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.measr.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	sum.acf.regen.samp <- matrix(nrow=comp.parms$max.lag, ncol=comp.parms$n.samp)
	post.fore.samp <- matrix(nrow=comp.parms$m, ncol=comp.parms$n.samp)
	# Implement remaining simulation
	for (i.samp in 1:comp.parms$n.samp) {
		# Update x0
		cand.model.parms$x0 <- propose.update.x0(curr.model.parms$x0, comp.parms$x0.rad)
		cand.model.parms$x.diff[1] <- model.parms$x[1] - cand.model.parms$x0
		cand.ls.stat <- reg.recalc.acorr(cand.model.parms, curr.ls.stat)
		cand.model.parms$logpdens <- get.logdgendens(cand.model.parms, cand.ls.stat)
		log.accept.prob <- cand.model.parms$logpdens$all - curr.model.parms$logpdens$all
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.x0.cnt <- accept.x0.cnt + 1
			curr.model.parms <- cand.model.parms
			curr.ls.stat <- cand.ls.stat
		} else {
			cand.model.parms <- curr.model.parms
			cand.ls.stat <- curr.ls.stat
		}
		# Update phi1
		cand.model.parms$phi1 <- propose.update.phi1(curr.model.parms$phi1, comp.parms$phi1.rad)
		cand.model.parms$acf.seq <- get.acf.ar1(cand.model.parms$phi1, comp.parms$max.lag)
		cand.ls.stat <- reg.summary.acorr(cand.model.parms)
		if (cand.ls.stat$kappa < comp.parms$max.kappa) {
			cand.model.parms$logpdens <- get.logdgendens(cand.model.parms, cand.ls.stat)
			log.accept.prob <- cand.model.parms$logpdens$all - curr.model.parms$logpdens$all
		} else {
			poorcond.cnt <- poorcond.cnt + 1
			log.accept.prob <- -Inf
		}
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.phi1.cnt <- accept.phi1.cnt + 1
			curr.model.parms <- cand.model.parms
			curr.ls.stat <- cand.ls.stat
		} else {
			cand.model.parms <- curr.model.parms
			cand.ls.stat <- curr.ls.stat
		}
		# Update gamma0 and beta
		curr.model.parms$gamma0 <- curr.ls.stat$SSE / rchisq(n=1, df=nu)
		curr.model.parms$delta <- curr.ls.stat$d.hat + curr.ls.stat$sqrt.v.inv*rnorm(n=1, mean=0, sd=sqrt(curr.model.parms$gamma0))
		curr.model.parms$logpdens <- get.logdgendens(curr.model.parms, curr.ls.stat)
		# Record sample values
		post.x0.samp[i.samp] <- curr.model.parms$x0
		post.phi1.samp[i.samp] <- curr.model.parms$phi1
		post.gamma0.samp[i.samp] <- curr.model.parms$gamma0
		post.delta.samp[i.samp] <- curr.model.parms$delta
		logdgdens.samp$pts[, i.samp] <- curr.model.parms$logpdens$pts
		logdgdens.samp$all[i.samp] <- curr.model.parms$logpdens$all
		# Regenerate the time series, summarize, and record summaries
		acf.measr.samp[, i.samp] <- pred.acf(curr.model.parms, comp.parms)
		regen.model.parms <- sim.regen.ts(curr.model.parms)
		fore.model.parms <- sim.future.ts(curr.model.parms, regen.model.parms, comp.parms)
		post.fore.samp[, i.samp] <- fore.model.parms$x
		acf.regen.samp[, i.samp] <- pred.acf(regen.model.parms, comp.parms)
		for (lag in 1:comp.parms$max.lag) {
			max.acf.measr.samp[lag, i.samp] <- max(abs(acf.measr.samp[1:lag, i.samp]))
			max.acf.regen.samp[lag, i.samp] <- max(abs(acf.regen.samp[1:lag, i.samp]))
			sum.acf.measr.samp[lag, i.samp] <- sum(abs(acf.measr.samp[1:lag, i.samp]))
			sum.acf.regen.samp[lag, i.samp] <- sum(abs(acf.regen.samp[1:lag, i.samp]))
		}
	}
	# Calculate DIC
	Bayes.model.parms <- model.parms
	Bayes.model.parms$x0 <- mean(post.x0.samp)
	Bayes.model.parms$x.diff <- c(model.parms$x[1] - Bayes.model.parms$x0, diff(model.parms$x))
	Bayes.model.parms$phi1 <- mean(post.phi1.samp)
	Bayes.model.parms$gamma0 <- exp(mean(log(post.gamma0.samp)))
	Bayes.model.parms$delta <- mean(post.delta.samp)
	Bayes.model.parms$acf.seq <- get.acf.ar1(Bayes.model.parms$phi1, comp.parms$max.lag)
	Bayes.ls.stat <- reg.summary.acorr(Bayes.model.parms)
	Bayes.model.parms$logpdens <- get.logdgendens(Bayes.model.parms, Bayes.ls.stat)
	p.DIC <- 2*(Bayes.model.parms$logpdens$all - mean(logdgdens.samp$all))
	DIC <- -2*Bayes.model.parms$logpdens$all + 2*p.DIC
	# Calculate WAIC
	mean.dgdens <- rowSums(exp(logdgdens.samp$pts)) / comp.parms$n.samp
	sum.log.mean.dgdens <- sum(log(mean.dgdens))
	mean.logdgdens <- rowSums(logdgdens.samp$pts) / comp.parms$n.samp
	p.WAIC <- 2*(sum.log.mean.dgdens - sum(mean.logdgdens))
	WAIC <- -2*sum.log.mean.dgdens + 2*p.WAIC
	# Calculate posterior predictive p-values
	pval.indiv.acf <- numeric(length=comp.parms$max.lag)
	pval.max.acf <- numeric(length=comp.parms$max.lag)
	pval.sum.acf <- numeric(length=comp.parms$max.lag)
	for (lag in 1:comp.parms$max.lag) {
		pval.indiv.acf[lag] <- length(which(abs(acf.regen.samp[lag,]) > abs(acf.measr.samp[lag,]))) / comp.parms$n.samp
		pval.max.acf[lag] <- length(which(max.acf.regen.samp[lag,] > max.acf.measr.samp[lag,])) / comp.parms$n.samp
		pval.sum.acf[lag] <- length(which(sum.acf.regen.samp[lag,] > sum.acf.measr.samp[lag,])) / comp.parms$n.samp
	}
	# Record results
	samp <- list(x0=post.x0.samp, phi1=post.phi1.samp, gamma0=post.gamma0.samp, delta=post.delta.samp, logdgdens=logdgdens.samp, fore=post.fore.samp)
	rate <- list(x0=accept.x0.cnt/comp.parms$n.samp, phi1=accept.phi1.cnt/comp.parms$n.samp, poorcond=poorcond.cnt/comp.parms$n.samp)
	stat <- list(DIC=DIC, p.DIC=p.DIC, WAIC=WAIC, p.WAIC=p.WAIC)
	pval <- list(indiv=pval.indiv.acf, max=pval.max.acf, sum=pval.sum.acf)
	result <- list(samp=samp, rate=rate, stat=stat, pval=pval)
	return(result)
}
