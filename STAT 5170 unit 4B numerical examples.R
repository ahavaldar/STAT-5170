# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part B of learning unit 4
# =================================================================================
#
# This set of numerical examples expands our exploration of Bayesian inference in time series analysis by introducing Markov chain Monte Carlo (MCMC) computational techiques, which can be especially helpful for estimating the parameters of time series models. We also explore a number of classical estimation techniques, some of which can be helpful for configuring Bayesian computational algorithms to inncrease efficiency.

# =================================================================================
# Analysis of binomial count data by MCMC
# =================================================================================
#
# In a previous numerical example, we implemented an example analysis of binimial count data based on a beta prior distribution. Because the setup is mathematically easy to manage, it was not necessary to implement a sophisticaed compuational algorithm: the posterior distribution is also beta, and may be accesed directly using R's built-in functionality for working with beta distributions.

# However, to familiarize ourselves with MCMC algorithms, it is worthwhile to consider how MCMC calculations could play out in such an analysis.

# Recall the setup of this example: the binomial count variable x has data-generating distribution x|theta ~ binomial(n, theta), and the prior distribution has theta ~ beta(alpha, beta). We know that the posterior distribution has theta|x ~ beta(alpha+x, beta+n-x), but we will pretend that information is not available.

# Example values of the relevant quantities are specified as follows.

n <- 10
x <- 3
alpha <- 1
beta <- 1

# The following user-defined function implements a Metropolis-Hastings algorithm to simulate a sample from the posterior distribtion. Observe that it makes use of several subfunctions to help with the more technical calculations. 

# The first subfunction generates a candidate parameter using a uniform shift with reflecting boundaries.

# The second subfunction evaluates the posterior density up to proportional equivalence. A logarithmic transformation is applied to manage potential computational imprecision.

# This function also makes use fo R's list functionality for passing arguements. Observe that the function has just two arguments, model.parms and comp.parms. Each of these, however, stores several varibles, the first associated with the data and model parameters, and the second associated with computational and tuning parameters that have to do with controlling the simulation algorithm. Variables within these lists are accessed using by stating the list name, followed by the dollar sign symbol ($), and then the relevant variable name. The output to the function is also a list.

run.mcmc.binom <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	propose.update <- function(curr.parm, radius) {
		cand.parm <- curr.parm + runif(n=1, min=-radius, max=radius)
		if (cand.parm < 0) {
			cand.parm <- -cand.parm
		} else if (cand.parm > 1) {
			cand.parm <- 1 - (cand.parm - 1)
		}
		return(cand.parm)
	}
	log.post.dens <- function(theta, x, n, alpha, beta) {
		val <- dbinom(x=x, size=n, prob=theta, log=TRUE) + dbeta(x=theta, shape1=alpha, shape2=beta, log=TRUE)
		return(val)
	}
	# ----- ----- --- main --- ----- -----
	x <- model.parms$x
	n <- model.parms$n
	curr.theta <- model.parms$theta
	alpha <- model.parms$alpha
	beta <- model.parms$beta
	n.samp <- comp.parms$n.samp
	theta.rad <- comp.parms$theta.rad
	accept.theta.cnt <- 0
	post.theta.samp <- numeric(length=n.samp)
	curr.log.post.dens <- log.post.dens(curr.theta, x, n, alpha, beta)
	for (i.samp in 1:n.samp) {
		# Update theta
		cand.theta <- propose.update(curr.theta, theta.rad)
		cand.log.post.dens <- log.post.dens(cand.theta, x, n, alpha, beta)
		log.accept.prob <- cand.log.post.dens - curr.log.post.dens
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.theta.cnt <- accept.theta.cnt + 1
			curr.theta <- cand.theta
			curr.log.post.dens <- cand.log.post.dens
		}
		post.theta.samp[i.samp] <- curr.theta
	}
	result <- list(samp=post.theta.samp, rate=accept.theta.cnt/n.samp)
	return(result)
}

# The following code takes the MCMC algorithm through a burn in phase. The model and computational parameters are specified as follows. Starting with the model parameters, the first several variables in that list store the data measurements and prior parameters.

model.parms <- list()
model.parms$x <- x
model.parms$n <- n
model.parms$alpha <- alpha
model.parms$beta <- beta

# The next variable stores the current value of the parameter. Since we are just setting up the algorihtm this value would be an initial value of theta. Let's try a random initial value.

model.parms$theta <- runif(n=1, min=0, max=1)
model.parms$theta

# Now we specify the computational parameters. Only two variables need to be specified in this list. The first is the algorithm's number of iterations. In the burn in phase, this number does not need to be very large--it will be set much larger in the implementation phase. The second parameter is the "radius" of the proposal distribution. This is the value c that specifies the width of the uniform shift with reflecting boundaries. In this numerical example, a large-ish value of this radius works quite well.

comp.parms <- list()
comp.parms$n.samp <- 100
comp.parms$theta.rad <- 0.95

# The following code carries out the burh in phase.

result <- run.mcmc.binom(model.parms, comp.parms)

# The list that is returned contains two variables. The first stores the simulated posterior sample. This is displayed as a histogram as follows

hist(result$samp, xlim=c(0, 1), xlab="theta", main="simulated posterior sample (burn in)")

# In the burn-in phase, this histogram is not necessary reliable as representing the posterior distribion. This is because a poor starting value may lead the sample to include a handful of values that are far from parameters with high posterior probabilitly.

# Plotting the trajectory of the sample values across iterations offers a useful diagnostic for checking whether the algorithm has moved into higher posterior probability regions, and would thereafter produce simulated samples that are representative of the posterior distribution.

# The following code produces a graph of this trajectory. The initial value is marked with a circle.

plot(1:comp.parms$n.samp, result$samp, type="l", lty=1, lwd=1, col="black", xlab="iteration", ylab="theta", main="")
points(1, result$samp[1], pch=1, cex=1)

# Here we see that very few values are simulated near the starting value. In other words, the algoritm has very quickly passed through the burn in stage. Some would say that it has "converged" to the posterior distribution. In any case, we are ready to implement the algorithm in its full capacity. The following code replaces the intial value with the last value simulated in the burn-in phase, increases the number of samples to be simulated, and runs the algorithm.

model.parms$theta <- result$samp[length(result$samp)]
comp.parms$n.samp <- 10000
result <- run.mcmc.binom(model.parms, comp.parms)

# A histogram of the simulated posterior sample is generated by the code below. Several posterior quantiles are also calculated and indicated in the plot. The posterior mean and standard deviation is also calculated.

qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
n.qnt <- length(qprobs)
post.quant <- as.numeric(quantile(result$samp, probs=qprobs))
for (i.qnt in 1:n.qnt) {
	print(paste(sprintf("%4.1f", 100*qprobs[i.qnt]), "% quantile: ", sprintf("%6.4f", post.quant[i.qnt]), sep=""))
}
print(paste("mean: ", sprintf("%6.4f", mean(result$samp)), sep=""))
print(paste("sdev: ", sprintf("%6.4f", sd(result$samp)), sep=""))
yref <- 1400
hist(result$samp, xlim=c(0, 1), xlab="theta", main="simulated posterior sample")
for (i.qnt in 1:n.qnt) {
	lines(post.quant[i.qnt]*c(1,1), yref*c(0.05,0.5), lty=1, lwd=3, col="red")
}

# The output to the MCMC function also includes a rate at which the Metropolis-Hastings step accepted the proposed simulated value across all iterations. This is converted into a rate by the following code.

print(paste("Acceptance rate for theta: ", sprintf("%6.4f", result$rate), sep=""))

# A common guideline is that the proposal distribtuon would be configured so that the acceptance rate is between 25% and 40%. We see that the rate is a little high in this example, that that is largely reflective of the problem's simplicity. In any case, the stated rates are only guidelines.

# We can begin to understand the relationship between the acceptance rate and proposal distribution by plotting the trajectory of the sample values across iterations and see how its patterns change in response to changing the radius of the proposal distribution (i.e., the parameter c). The following code extends the MCMC simulation at two distinct settings, and graphs the trajectories. For comparison it also graphs the trajectory of the simulation algorithm explored in the previous set of numerical examples. The latter algorithm is the ideal case in which sample values are simulated independently across iterations. 

comp.parms$n.samp <- 300

par(mfrow = c(3, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right

comp.parms$theta.rad <- 0.95
model.parms$theta <- result$samp[length(result$samp)]
result <- run.mcmc.binom(model.parms, comp.parms)
plot(1:comp.parms$n.samp, result$samp, type="l", lty=1, lwd=1, col="black", xlab="", ylab="", main=paste("accept rate=", sprintf("%4.1f", 100*result$rate), "% (c=", comp.parms$theta.rad, ")", sep=""))

comp.parms$theta.rad <- 0.05
model.parms$theta <- result$samp[length(result$samp)]
result <- run.mcmc.binom(model.parms, comp.parms)
plot(1:comp.parms$n.samp, result$samp, type="l", lty=1, lwd=1, col="black", xlab="", ylab="", main=paste("accept rate=", sprintf("%4.1f", 100*result$rate), "% (c=", comp.parms$theta.rad, ")", sep=""))

plot(1:comp.parms$n.samp, rbeta(n=comp.parms$n.samp, shape1=alpha+x, shape2=beta+n-x), type="l", lty=1, lwd=1, col="black", xlab="", ylab="", main=paste("direct sampling", sep=""))

par(mfrow = c(1, 1))

# It is suggested in this plot that a smaller radius results in a higher acceptance rate, which is associated with a trajectory that moves more slowly through the parameter space. The radius that yields an acceptance rate of about 40% moves through the paramerer space more thoroughly; the corresponding trajectory more closely resembles the ideal case of direct sampling.

result$samp <- rbeta(n=comp.parms$n.samp, shape1=alpha+x, shape2=beta+n-x)

qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
n.qnt <- length(qprobs)
post.quant <- as.numeric(quantile(result$samp, probs=qprobs))
for (i.qnt in 1:n.qnt) {
	print(paste(sprintf("%4.1f", 100*qprobs[i.qnt]), "% quantile: ", sprintf("%6.4f", post.quant[i.qnt]), sep=""))
}
print(paste("mean: ", sprintf("%6.4f", mean(result$samp)), sep=""))
print(paste("sdev: ", sprintf("%6.4f", sd(result$samp)), sep=""))
yref <- 1400
hist(result$samp, xlim=c(0, 1), xlab="theta", main="simulated posterior sample")
for (i.qnt in 1:n.qnt) {
	lines(post.quant[i.qnt]*c(1,1), yref*c(0.05,0.5), lty=1, lwd=3, col="red")
}

# =================================================================================
# Analysis of AR(2) time-series data by MCMC
# =================================================================================

# To demonstrate the use of Bayesian methods for estimating the parameters of an autorgressive model, let use work with the fish population recruitment time series from previous examples. Recall that this time series among the data sets of the "astsa" R package, and is stored in the time-series object named "rec". A plot of the data is as follows.

data(data=rec, package="astsa")
ts.plot(rec, xlab="year", ylab="count", main="Fish population recruitment")

# These data are described as the "recruitment (index of the number of new fish) for a period of 453 months ranging over the years 1950-1987. Recruitment is loosely defined as an indicator of new members of a population to the first life stage at which natural mortality stabilizes near adult levels."

# In previous work with this time series, we applied an AR(2) model, as we will here.

# ---------------------------------------------------------------------------------
# Parameter estimation
# ---------------------------------------------------------------------------------

# A suitable MCMC algorithm for analysis of autoregressive time series is coded into the user-defined funtion below. As with the previous example, this function takes input arguemnts in the form of lists, and returns an output variable in that same form.

# The structure of the algorithm can be seen to follow Gibbs sampling, wherein each step in the cycle updates a distinct set of parameters. Updating itself is implemented using a Metropolis-Hastings step.

# The algorithm is set up to include initial unmeasured data values among the parameters of the model. It is thus suitable for analysis projects such that inference on "starting values" is part of the project's objectives. The MCMC algorithms used in the upcoming learning unit 5 do not have this feature.

# The algorithm not only simulates values of each parameter from the posterior distribution, but also simulates predicted values up to given number of steps ahead.

# Several subfunctions are defined to assist with the technical calculations. Among them, notice that two subfunctions are defined for proposing a candidate parameter value, one as a uniform shift and the other as a uniform shift with a reflecting lower boundary.

# The input arguments include variables for each set of parameter, and for tuning values associated with each set of parameters.

run.mcmc.ar <- function(model.parms, comp.parms) {
	# ----- ----- subfunctions ----- -----
	log.dgen.dens <- function(x, sigw2, mu, x0, phi) {
		n <- length(x)
		p <- length(phi)
		x.ext <- c(x0, x) - mu
		Usum <- 0
		for (t in (p+1):(n+p)) {
			Usum <- Usum + (x.ext[t] - sum(phi*x.ext[seq(from=t-1, to=t-p)]))^2
		}
		val <- -0.5*n*log(2*pi*sigw2) - 0.5*Usum/sigw2
		return(val)
	}
	log.prior.dens <- function(sigw2, mu, x0, phi) {
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
	# ----- ----- --- main --- ----- -----
	x <- model.parms$x
	n <- length(model.parms$x)
	p <- length(model.parms$phi)
	m <- model.parms$m
	curr.sigw2 <- model.parms$sigw2
	curr.mu <- model.parms$mu
	curr.x0 <- model.parms$x0
	curr.phi <- model.parms$phi
	n.samp <- comp.parms$n.samp
	sigw2.rad <- comp.parms$sigw2.rad
	mu.rad <- comp.parms$mu.rad
	x0.rad <- comp.parms$x0.rad
	phi.rad <- comp.parms$phi.rad
	cand.x0 <- curr.x0
	cand.phi <- curr.phi
	accept.sigw2.cnt <- 0
	accept.mu.cnt <- 0
	accept.x0.cnt <- rep(x=0, times=p)
	accept.phi.cnt <- rep(x=0, times=p)
	post.sigw2.samp <- numeric(length=n.samp)
	post.mu.samp <- numeric(length=n.samp)
	post.x0.samp <- matrix(data=NA, nrow=p, ncol=n.samp)
	post.phi.samp <- matrix(data=NA, nrow=p, ncol=n.samp)
	post.predx.samp <- matrix(data=NA, nrow=m, ncol=n.samp)
	curr.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.mu, curr.x0, curr.phi)
	curr.log.prior.dens <- log.prior.dens(curr.sigw2, curr.mu, curr.x0, curr.phi)
	curr.log.post.dens <- curr.log.dgen.dens + curr.log.prior.dens
	for (i.samp in 1:n.samp) {
		# Update sigw2
		cand.sigw2 <- propose.update.B(curr.sigw2, sigw2.rad)
		cand.log.dgen.dens <- log.dgen.dens(x, cand.sigw2, curr.mu, curr.x0, curr.phi)
		cand.log.prior.dens <- log.prior.dens(cand.sigw2, curr.mu, curr.x0, curr.phi)
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
		# Update mu
		cand.mu <- propose.update.A(curr.mu, mu.rad)
		cand.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, cand.mu, curr.x0, curr.phi)
		cand.log.prior.dens <- log.prior.dens(curr.sigw2, cand.mu, curr.x0, curr.phi)
		cand.log.post.dens <- cand.log.dgen.dens + cand.log.prior.dens
		log.accept.prob <- cand.log.post.dens - curr.log.post.dens
		unif <- runif(n=1, min=0, max=1)
		if (unif <= exp(log.accept.prob)) {
			accept.mu.cnt <- accept.mu.cnt + 1
			curr.mu <- cand.mu
			curr.log.dgen.dens <- cand.log.dgen.dens
			curr.log.prior.dens <- cand.log.prior.dens
			curr.log.post.dens <- cand.log.post.dens
		}
		# Update x0
		for (j in 1:p) {
			cand.x0[j] <- propose.update.A(curr.x0[j], x0.rad[j])
			cand.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.mu, cand.x0, curr.phi)
			cand.log.prior.dens <- log.prior.dens(curr.sigw2, curr.mu, cand.x0, curr.phi)
			cand.log.post.dens <- cand.log.dgen.dens + cand.log.prior.dens
			log.accept.prob <- cand.log.post.dens - curr.log.post.dens
			unif <- runif(n=1, min=0, max=1)
			if (unif <= exp(log.accept.prob)) {
				accept.x0.cnt[j] <- accept.x0.cnt[j] + 1
				curr.x0[j] <- cand.x0[j]
				curr.log.dgen.dens <- cand.log.dgen.dens
				curr.log.prior.dens <- cand.log.prior.dens
				curr.log.post.dens <- cand.log.post.dens
			}
		} #for (j in 1:p
		# Update phi
		for (j in 1:p) {
			cand.phi[j] <- propose.update.A(curr.phi[j], phi.rad[j])
			cand.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.mu, curr.x0, cand.phi)
			cand.log.prior.dens <- log.prior.dens(curr.sigw2, curr.mu, curr.x0, cand.phi)
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
		pastx <- x[seq(from=n, to=n-p+1)] - curr.mu
		for (t in 1:m) {
			predx.mean <- sum(curr.phi*pastx)
			predx.sd <- sqrt(curr.sigw2)
			predx[t] <- rnorm(n=1, mean=predx.mean, sd=predx.sd)
			pastx <- c(predx[t], pastx[1:p-1])
		}
		predx <- predx + curr.mu
		# Record sample values
		post.sigw2.samp[i.samp] <- curr.sigw2
		post.mu.samp[i.samp] <- curr.mu
		post.x0.samp[,i.samp] <- curr.x0
		post.phi.samp[,i.samp] <- curr.phi
		post.predx.samp[,i.samp] <- predx
	}
	samp <- list(sigw2=post.sigw2.samp, mu=post.mu.samp, x0=post.x0.samp, phi=post.phi.samp, predx=post.predx.samp)
	rate <- list(sigw2=accept.sigw2.cnt/n.samp, mu=accept.mu.cnt/n.samp, x0=accept.x0.cnt/n.samp, phi=accept.phi.cnt/n.samp)
	result <- list(samp=samp, rate=rate)
	return(result)
}

# To implement the algorithm, a first step is to calculate the Yule-Walker estimates and assign them as as initial values of the parameters. These estimates are calculated by the following code, which obtains sample autocovariance and autocorrelation values from R's built-in acf function, which subsequently defines the required vector-matrix objects.

acf.seq <- acf(rec, plot=FALSE)
acov.seq <- acf(rec, type="covariance", plot=FALSE)

p <- 2
gamma0 <- acov.seq$acf[1]
gamma.vect <- acov.seq$acf[1+1:p]
Gamma.mat <- matrix(data=0, nrow=p, ncol=p)
Gamma.mat[, 1] <- acov.seq$acf[1:p]
for (j in 2:p) {
	Gamma.mat[1:(j-1),j] <- Gamma.mat[seq(from=j, to=2),1]
	Gamma.mat[j:p,j] <- Gamma.mat[1:(p-j+1),1]
}
Gamma.inv <- solve(Gamma.mat)

phi.est <- Gamma.inv %*% gamma.vect
phi.est

sigw2.est <- as.numeric(gamma0 - t(phi.est) %*% gamma.vect)
sigw2.est

phi.sderr <- sqrt(sigw2.est*diag(Gamma.inv)[1:p] / length(rec))
phi.sderr

# The model parameters and initial values are specified below. The Yule-Walker estimates provide initial values for sigw2 and phi. Initial values for mu and the unmeasured data are set to the sample mean.

model.parms <- list()
model.parms$x <- as.numeric(rec)
model.parms$m <- 15
model.parms$sigw2 <- sigw2.est
model.parms$mu <- mean(model.parms$x)
model.parms$x0 <- rep(x=model.parms$mu, times=p)
model.parms$phi <- as.numeric(phi.est)

# Computational and tuning parameters are specified as follows. The Yule-Walker standard errors provide values of the proposal radius for phi. The proposal radius for sigw2 is specified at a guess. Proposal radius values for mu and the unmeasured data are set to the usual standard error of the sample mean. Note that the number of iterations is set to a small value for burining in the algorithm.

comp.parms <- list()
comp.parms$n.samp <- 100
comp.parms$sigw2.rad <- 50
comp.parms$mu.rad <- sd(model.parms$x)/length(model.parms$x)
comp.parms$x0.rad <- rep(x=comp.parms$mu.rad, times=p)
comp.parms$phi.rad <- phi.sderr

# The burn in phase is implemented as follows

result <- run.mcmc.ar(model.parms, comp.parms)

# The code below plots the trajectories of simulated values for each parameter, and checks the acceptance rates.

par(mfrow = c(2, 3))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(1:comp.parms$n.samp, result$samp$sigw2, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main="sigw2")
plot(1:comp.parms$n.samp, result$samp$mu, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main="mu")
for (j in 1:p) {
	plot(1:comp.parms$n.samp, result$samp$x0[j,], type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("x0[", j, "]", sep=""))
}
for (j in 1:p) {
	plot(1:comp.parms$n.samp, result$samp$phi[j,], type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi[", j, "]", sep=""))
}
par(mfrow = c(1, 1))

# The algorithm seems to have drifted into a region of high posterior probability, but the acceptance rates, obtained by the code below, suggest that further tuning may be helpful.

print(paste("Acceptance rate for sigw2: ", sprintf("%6.4f", result$rate$sigw2), sep=""))
print(paste("Acceptance rate for mu: ", sprintf("%6.4f", result$rate$mu), sep=""))
for (j in 1:p) {
	print(paste("Acceptance rate for x0[", j, "]: ", sprintf("%6.4f", result$rate$x0[j]), sep=""))
}
for (j in 1:p) {
	print(paste("Acceptance rate for phi[", j, "]: ", sprintf("%6.4f", result$rate$phi[j]), sep=""))
}

# The following code continues to run the algorithm from the current parameter values of the burn in step.

model.parms$sigw2 <- result$samp$sigw2[comp.parms$n.samp]
model.parms$mu <- result$samp$mu[comp.parms$n.samp]
model.parms$x0 <- result$samp$x0[,comp.parms$n.samp]
model.parms$phi <- result$samp$phi[,comp.parms$n.samp]
result <- run.mcmc.ar(model.parms, comp.parms)

# This next block of code adjusts the proposal radii to a multiple of the posterior standard deviations of the inital sample. This approach can often work well for tuning the algorithm, but it still requires attentive effort. 

comp.parms$sigw2.rad <- 3.0*sd(result$samp$sigw2)
comp.parms$mu.rad <- 3.5*sd(result$samp$mu)
for (j in 1:p) {
	comp.parms$x0.rad[j] <- 1.25*sd(result$samp$x0[j,])
}
for (j in 1:p) {
	comp.parms$phi.rad[j] <- 1.25*sd(result$samp$phi[j,])
}

# The algorithm is now run a little longer, starting from where the last run left off.

model.parms$sigw2 <- result$samp$sigw2[comp.parms$n.samp]
model.parms$mu <- result$samp$mu[comp.parms$n.samp]
model.parms$x0 <- result$samp$x0[,comp.parms$n.samp]
model.parms$phi <- result$samp$phi[,comp.parms$n.samp]
comp.parms$n.samp <- 300
result <- run.mcmc.ar(model.parms, comp.parms)

# The trajectories look a little better.

par(mfrow = c(2, 3))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(1:comp.parms$n.samp, result$samp$sigw2, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main="sigw2")
plot(1:comp.parms$n.samp, result$samp$mu, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main="mu")
for (j in 1:p) {
	plot(1:comp.parms$n.samp, result$samp$x0[j,], type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("x0[", j, "]", sep=""))
}
for (j in 1:p) {
	plot(1:comp.parms$n.samp, result$samp$phi[j,], type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi[", j, "]", sep=""))
}
par(mfrow = c(1, 1))

# ... and, if we loop through this adjustment step a few times, so do the acceptance rates

print(paste("Acceptance rate for sigw2: ", sprintf("%6.4f", result$rate$sigw2), sep=""))
print(paste("Acceptance rate for mu: ", sprintf("%6.4f", result$rate$mu), sep=""))
for (j in 1:p) {
	print(paste("Acceptance rate for x0[", j, "]: ", sprintf("%6.4f", result$rate$x0[j]), sep=""))
}
for (j in 1:p) {
	print(paste("Acceptance rate for phi[", j, "]: ", sprintf("%6.4f", result$rate$phi[j]), sep=""))
}

# The algorithm now appears to be burnt-in (i.e., it evolves within a region of high posterior probability) and well tuned, which will improve the accuracy of results. As a final step before summarizing our results, we run the algorithm for an extended time, starting from the current parameter values.

model.parms$sigw2 <- result$samp$sigw2[comp.parms$n.samp]
model.parms$mu <- result$samp$mu[comp.parms$n.samp]
model.parms$x0 <- result$samp$x0[,comp.parms$n.samp]
model.parms$phi <- result$samp$phi[,comp.parms$n.samp]
comp.parms$n.samp <- 10000
result <- run.mcmc.ar(model.parms, comp.parms)

# Given the effort we have put in, it is worthwhile at this point to save our results.

setwd(DataDirectory)
#save(model.parms, comp.parms, result, file="fish count MCMC.RData")
load(file="fish count MCMC.RData") #model.parms, comp.parms, result

# The following code produces histograms of the simulated posterior samples, one for each parameter, overlaid with several posterior quantiles. Posterior means and standard deviations are also calculated. A helper function is defined to manage the calculations.

plot.hist <- function(label, samp, yref=NA) {
	print(paste("*** ", label, " ***", sep=""))
	qprobs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
	n.qnt <- length(qprobs)
	post.quant <- as.numeric(quantile(samp, probs=qprobs))
	for (i.qnt in 1:n.qnt) {
		print(paste(sprintf("%4.1f", 100*qprobs[i.qnt]), "% quantile: ", sprintf("%6.4f", post.quant[i.qnt]), sep=""))
	}
	print(paste("          mean: ", sprintf("%6.4f", mean(samp)), sep=""))
	print(paste("          sdev: ", sprintf("%6.4f", sd(samp)), sep=""))
	hist(samp, xlab=label, main=paste(label, " posterior sample", sep=""))
	if (!is.na(yref)) {
		for (i.qnt in 1:n.qnt) {
			lines(post.quant[i.qnt]*c(1,1), yref*c(0.05,0.5), lty=1, lwd=3, col="red")	
		}
	}
}

par(mfrow = c(2, 3))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot.hist(label="sigw2", samp=result$samp$sigw2, yref=1200)
plot.hist(label="mu", samp=result$samp$mu, yref=1500)
yref.seq <- c(2500, 1500)
for (j in 1:p) {
	plot.hist(label=paste("x0[", j, "]", sep=""), samp=result$samp$x0[j,], yref=yref.seq[j])
}
yref.seq <- c(1500, 2000)
for (j in 1:p) {
	plot.hist(label=paste("phi[", j, "]", sep=""), samp=result$samp$phi[j,], yref=yref.seq[j])
}
par(mfrow = c(1, 1))

# The following code generates a plot of the last portion of the example time series, with simulated forecasts and quantiles

t1 <- 400
max.m <- dim(result$samp$predx)[1]
qprobs <- c(0.025, 0.5, 0.975)
n.qnt <- length(qprobs)
qmat <- matrix(data=NA, nrow=n.qnt, ncol=max.m)
for (m in 1:max.m) {
	qmat[,m] <- as.numeric(quantile(result$samp$predx[m,], probs=qprobs))
}
n <- length(model.parms$x)
ymin <- min(model.parms$x, qmat)
ymax <- max(model.parms$x, qmat)
yrng <- ymax-ymin
plot(t1:n, model.parms$x[t1:n], type="l", lty=1, lwd=1, col="black", xlim=c(t1-1,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="previous", ylab="current", main="")
points(t1:n, model.parms$x[t1:n], pch=1, cex=1, col="black")
lines(n:(n+max.m), c(model.parms$x[n],qmat[2,]), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), qmat[2,], pch=1, cex=1, col="red")
lines(n:(n+max.m), c(model.parms$x[n],qmat[1,]), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(model.parms$x[n],qmat[3,]), lty=2, lwd=1, col="blue")

