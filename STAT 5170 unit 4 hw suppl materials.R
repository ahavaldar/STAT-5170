# =================================================================================
# STAT 5170: Applied Time Series
# Suppleental materials for learning unit 4 homework
# =================================================================================
#
# ---------------------------------------------------------------------------------
# Function that calculates the switch-count of a zero-one sequence
# ---------------------------------------------------------------------------------

test.stat <- function(x.seq) {
	n <- length(x.seq)
	stat <- sum(abs(x.seq[2:n] - x.seq[1:(n-1)]))
	return(stat)
}

# ---------------------------------------------------------------------------------
# Function that calculates the mean of adjacent-residual absolute sums
# ---------------------------------------------------------------------------------

test.stat <- function(res) {
	n <- length(res)
	stat <- mean(abs(res[1:n-1] + res[2:n]))
	return(stat)
}

# ---------------------------------------------------------------------------------
# Function that implements MCMC calculation for Bayesian analysis of a zero-mean AR(1) times series with initial value zero
# ---------------------------------------------------------------------------------

run.mcmc.m0.ar1 <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	log.dgen.dens <- function(x, sigw2, x0, phi) {
		n <- length(x)
		p <- length(phi)
		x.ext <- c(x0, x)
		Usum <- 0
		for (t in (p+1):(n+p)) {
			Usum <- Usum + (x.ext[t] - sum(phi*x.ext[seq(from=t-1, to=t-p)]))^2
		}
		val <- -0.5*n*log(2*pi*sigw2) - 0.5*Usum/sigw2
		return(val)
	}
	log.prior.dens <- function(sigw2, x0, phi) {
		val <- -log(sigw2)
		return(val)
	}
	propose.update.A <- function(curr.parm, radius) {
		cand.parm <- curr.parm + runif(n=1, min=-radius, max=radius)
		return(cand.parm)
	}
	propose.update.B <- function(curr.parm, radius) {
		cand.parm <- abs(curr.parm + runif(n=1, min=-radius, max=radius))
		return(cand.parm)
	}
	propose.update.C <- function(curr.parm, radius) {
		cand.parm <- curr.parm + runif(n=1, min=-radius, max=radius)
		if (cand.parm < -1) {
			cand.parm <- -1 - (cand.parm + 1)
		} else if (cand.parm > 1) {
			cand.parm <- 1 - (cand.parm - 1)
		}
		return(cand.parm)
	}
	# ----- ----- --- main --- ----- -----
	x <- model.parms$x
	n <- length(model.parms$x)
	p <- length(model.parms$phi)
	m <- model.parms$m
	curr.sigw2 <- model.parms$sigw2
	curr.x0 <- 0
	curr.phi <- model.parms$phi
	n.samp <- comp.parms$n.samp
	sigw2.rad <- comp.parms$sigw2.rad
	phi.rad <- comp.parms$phi.rad
	cand.phi <- curr.phi
	accept.sigw2.cnt <- 0
	accept.phi.cnt <- rep(x=0, times=p)
	post.sigw2.samp <- numeric(length=n.samp)
	post.phi.samp <- matrix(data=NA, nrow=p, ncol=n.samp)
	post.predx.samp <- matrix(data=NA, nrow=m, ncol=n.samp)
	curr.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.x0, curr.phi)
	curr.log.prior.dens <- log.prior.dens(curr.sigw2, curr.x0, curr.phi)
	curr.log.post.dens <- curr.log.dgen.dens + curr.log.prior.dens
	for (i.samp in 1:n.samp) {
		# Update sigw2
		cand.sigw2 <- propose.update.B(curr.sigw2, sigw2.rad)
		cand.log.dgen.dens <- log.dgen.dens(x, cand.sigw2, curr.x0, curr.phi)
		cand.log.prior.dens <- log.prior.dens(cand.sigw2, curr.x0, curr.phi)
		cand.log.post.dens <- cand.log.dgen.dens + cand.log.prior.dens
		log.accept.prob <- cand.log.post.dens - curr.log.post.dens
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.sigw2.cnt <- accept.sigw2.cnt + 1
			curr.sigw2 <- cand.sigw2
			curr.log.dgen.dens <- cand.log.dgen.dens
			curr.log.prior.dens <- cand.log.prior.dens
			curr.log.post.dens <- cand.log.post.dens
		}
		# Update phi
		for (j in 1:p) {
			cand.phi[j] <- propose.update.C(curr.phi[j], phi.rad[j])
			cand.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.x0, cand.phi)
			cand.log.prior.dens <- log.prior.dens(curr.sigw2, curr.x0, cand.phi)
			cand.log.post.dens <- cand.log.dgen.dens + cand.log.prior.dens
			log.accept.prob <- cand.log.post.dens - curr.log.post.dens
			unif <- runif(n=1, min=0, max=1)
			if (unif <= exp(log.accept.prob)) {
				accept.phi.cnt[j] <- accept.phi.cnt[j] + 1
				curr.phi[j] <- cand.phi[j]
				curr.log.dgen.dens <- cand.log.dgen.dens
				curr.log.prior.dens <- cand.log.prior.dens
				curr.log.post.dens <- cand.log.post.dens
			}
		} #for (j in 1:p
		# Simulate future values
		predx <- numeric(length=m)
		pastx <- x[seq(from=n, to=n-p+1)]
		for (t in 1:m) {
			predx.mean <- sum(curr.phi*pastx)
			predx.sd <- sqrt(curr.sigw2)
			predx[t] <- rnorm(n=1, mean=predx.mean, sd=predx.sd)
			pastx <- c(predx[t], pastx[1:p-1])
		}
		# Record sample values
		post.sigw2.samp[i.samp] <- curr.sigw2
		post.phi.samp[,i.samp] <- curr.phi
		post.predx.samp[,i.samp] <- predx
	}
	samp <- list(sigw2=post.sigw2.samp, phi=post.phi.samp, predx=post.predx.samp)
	rate <- list(sigw2=accept.sigw2.cnt/n.samp, phi=accept.phi.cnt/n.samp)
	result <- list(samp=samp, rate=rate)
	return(result)
}

