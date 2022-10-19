# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part A of learning unit 3
# =================================================================================
#
# =================================================================================
# Infinite-sum representations, truncated
# =================================================================================
#
# Infinite-sum representations of a time series are helpful mainly for theoretical purposes, e.g., to derive objects like the autocorrelation function. We will see, later in the course, a few practical uses for them. However, it can be worthwhile to explore them using simulation.

# We know, for example, that an AR(1) time series with autoregressive parameter phi1 had an infinite moving-average representation with coefficients psi_j = phi1^j. The following code simulates the AR(1) time series, and approximation to it using a truncated moving-average representation. Values for the white-noise standard deviation and autoregressive parameter are specfied as follows:

sigw <- 3.5
phi <- 0.9

# To produce the sample paths, we make use of the user-defined functions for simulating autoregressive and moving-average time series.  

ar.sim <- function(n, phi, sigw, x0=NA) {
	p <- length(phi)
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	if (is.na(x0)) {
		pastx <- rep(x=0, times=p)
	} else {
		pastx <- c(x0, rep(x=0, times=p-length(x0)))
	}
	ar.samp <- numeric(length=n)
	for (t in 1:n) {
		ar.samp[t] <- sum(phi*pastx) + wn.samp[t]
		pastx <- c(ar.samp[t], pastx[1:p-1])
	}
	ar.ts <- ts(ar.samp)
	return(ar.ts)
}

ma.sim <- function(n, theta, sigw, w0=NA) {
	q <- length(theta)
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	if (is.na(w0)) {
		pastw <- rep(x=0, times=q)
	} else {
		pastw <- c(w0, rep(x=0, times=q-length(w0)))
	}
	ma.samp <- numeric(length=n)
	for (t in 1:n) {
		ma.samp[t] <- sum(theta*pastw) + wn.samp[t]
		pastw <- c(wn.samp[t], pastw[1:q-1])
	}
	ma.ts <- ts(ma.samp)
	return(ma.ts)
}

# ---------------------------------------------------------------------------------
# Sample paths
# ---------------------------------------------------------------------------------

# A sample path from the AR(1) time series is simulated as follows

n <- 100
ar.ts <- ar.sim(n=n, phi=phi, sigw=sigw)
ts.plot(ar.ts, xlab="time", ylab="value", main=paste("AR(", length(phi), ")", sep=""))

# The moving-average version of the time series is simulated using just the first q=25 parameters of the infinite-sum representation. Shown next to each other, the two time series are seen to exhibit similar patterns.

q <- 25
theta <- phi^(1:q)
ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw)

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(ar.ts, xlab="time", ylab="value", main=paste("AR(", length(phi), ")", sep=""))
ts.plot(ma.ts, xlab="time", ylab="value", main=paste("MA(", length(theta), ")", sep=""))
par(mfrow = c(1, 1))

# Moreover, the sample ACFs and PACFs of the two time series are similar, though not exactly the same.

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(ar.ts)
pacf(ar.ts)
acf(ma.ts)
pacf(ma.ts)
par(mfrow = c(1, 1))

# The differences are partly due to simulation variability, and partly due to innacuracy of the approximation.

# In the infinite moving-average representation the moving-average parameters decrease exponentally to zero. We see in the truncated moving-average representation that the moving-average parameters come close to zero, but even with as many as q=25 values calculated, the last parameter is still noticeably different from zero.

theta

# If q is much smaller then distinctions between the two simulaed time-series become more apparent.

q <- 10
theta <- phi^(1:q)
ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw)

# The distinctions are not easy to see in the sample paths...

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(ar.ts, xlab="", ylab="value", main=paste("AR(", length(phi), ")", sep=""))
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))
par(mfrow = c(1, 1))

# ... but they are clear in the sample ACFs and PACFs.

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(ar.ts)
pacf(ar.ts)
acf(ma.ts)
pacf(ma.ts)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Non-zero mean time series
# ---------------------------------------------------------------------------------

# The following is an exploration of the AR(1) time series with a non-zero mean value specified as follows.

mu <- 50

# When approximating the AR(1) time series with a truncated moving-average representation, it is straightforward to incorporate the non-zero mean value into simulation. This is demonstrated in the following code (which extracts from the previous user-defined functions).

n <- 100
q <- 25
theta <- phi^(1:q)
wn.samp <- rnorm(n=n, mean=0, sd=sigw)
pastw <- rep(x=0, times=q)
ma.samp <- numeric(length=n)
for (t in 1:n) {
	ma.samp[t] <- mu + sum(theta*pastw) + wn.samp[t]
	pastw <- c(wn.samp[t], pastw[1:q-1])
}
ma.ts <- ts(ma.samp)
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))
abline(h=mu, lty=1, lwd=1, col="blue")

# Some nuance is requried to incorporate a non-zero mean into a direct simulation of the AR(1) time series. This is demonstrated as follows.

alpha <- (1-sum(phi))*mu
p <- length(phi)
wn.samp <- rnorm(n=n, mean=0, sd=sigw)
pastx <- rep(x=mu, times=p)
ar.samp <- numeric(length=n)
for (t in 1:n) {
	ar.samp[t] <- alpha + sum(phi*pastx) + wn.samp[t]
	pastx <- c(ar.samp[t], pastx[1:p-1])
}
ar.ts <- ts(ar.samp)

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(ar.ts, xlab="", ylab="value", main=paste("AR(", length(phi), ")", sep=""))
abline(h=mu, lty=1, lwd=1, col="blue")
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))
abline(h=mu, lty=1, lwd=1, col="blue")
par(mfrow = c(1, 1))

# =================================================================================
# ACF of a AR(p) time series
# =================================================================================

# In this example, we work with the difference equations that determine the ACF of an autoregressive time series.

# Supppse the time series is AR(p) with p=3 and autoregressive parameters phi1, phi2, and phi3 specified as follows.

p <- 3
phi1 <- 0.6
phi2 <- 0.25
phi3 <- -0.15

# Suppose further than the white-noise standard deviation, sigw, is specified as follows.

sigw <- 5

# ---------------------------------------------------------------------------------
# Initial ACF values
# ---------------------------------------------------------------------------------

# It is convenient to solve for the initial values by setting up the a matrix equation. 

# The difference equations to be solved are
#
# rho(0) = 1
# rho(1) = phi1*rho(0) + phi2*rho(1) + phi3*rho(2)
# rho(2) = phi1*rho(1) + phi2*rho(0) + phi3*rho(1)
#
# A suitable rearrangement that suggests a matrix formula is
#
# rho(0) = 1
# phi1*rho(0) - (1-phi2)*rho(1) + phi3*rho(2) = 0
# phi2*rho(0) + (phi1+phi3)*rho(1) - rho(2) = 0
#
# The matrix version of the equation is
#
# phi.mat %*% rho.vect = c.vect
#
# for the p x p matrix phi.mat and p x 1 vectors rho.vect and c.vect set up as follows.

phi.mat <- matrix(data=0, nrow=p, ncol=p)
phi.mat[1,1] <- 1
phi.mat[2,1] <- phi1
phi.mat[2,2] <- -(1-phi2)
phi.mat[2,3] <- phi3
phi.mat[3,1] <- phi2
phi.mat[3,2] <- phi1+phi3
phi.mat[3,3] <- -1
phi.mat

c.vect <- matrix(data=0, nrow=p, ncol=1)
c.vect[1,1] <- 1
c.vect

# The solution is 
#
# rho.vect = phi.inv %*% c.vect
#
# where phi.inv is the matrix inverse of phi.mat. This is calculated as follows.

phi.inv <- solve(phi.mat)
rho.vect <- phi.inv %*% c.vect
rho.vect

# Reading the entries of rho.vect, we see that
#
# rho(0) = 1
# rho(1) = 0.69
# rho(2) = 0.56

# ---------------------------------------------------------------------------------
# Additional ACF values
# ---------------------------------------------------------------------------------

# Once the inital ACF values are available, additional values are then calculated recursively from the difference equation
#
# rho(h) = phi1*rho(h-1) + phi2*rho(h-2) + phi3*rho(h-3)
#

max.lag <-9
phi <- c(phi1, phi2, phi3)
acf.seq <- c(as.numeric(rho.vect), rep(x=0, times=max.lag-p))
for (lag in p:max.lag) {
	acf.seq[lag+1] <- sum(phi*acf.seq[seq(from=lag, to=lag-p+1, by=-1)])
}
acf.seq

# We can check these values by simulating the time series and calculating the sample ACF.

ar.ts <- ar.sim(n=10000, phi=phi, sigw=sigw)
acf.chk <- acf(ar.ts, plot=FALSE)
cbind(acf.seq, acf.chk$acf[1:(max.lag+1)])

# ---------------------------------------------------------------------------------
# Variance
# ---------------------------------------------------------------------------------

# The variance of an xt value in the time series is given by the equation
#
# gamma(0) = sigw^2 / (1 - phi1*rho(1) - phi2*rho(2) - phi3*rho(3))
#

# This is calclated as follows

phi <- c(phi1, phi2, phi3)
rho <- as.numeric(acf.seq[2:(p+1)])
gamma0 <- sigw^2 / (1 - sum(phi*rho))
gamma0

# Compare the white-noise standard deviation with the standard deviation of an xt value:

sigw
sqrt(gamma0)

# We can check these values from the previously simulated sample path.

sd(ar.ts)

# =================================================================================
# Partial autocorrelation
# =================================================================================

# This next portion of the example explores some of the concepts and computational steps associated with partial autocorrelation.  

# ---------------------------------------------------------------------------------
# Regression framework
# ---------------------------------------------------------------------------------

# Partial autocorrelation concepts share some similarities with those of linear regression. Consider an AR(1) time series with parameters specified as follows.

sigw <- 5.2
phi <- 0.9

# The following code simualtes a sample path of this series.

n <- 100
ar.ts <- ar.sim(n=n, phi=phi, sigw=sigw)

# The following code generates a scatterplot of current time series values plotted against those at one time-point before (starting at time t=2).

plot(ar.ts[1:(n-1)], ar.ts[2:n], type="p", pch=16, cex=1, col="black", xlab="previous", ylab="current", main="")
abline(a=0, b=phi, lty=2, lwd=1, col="blue")
abline(a=0, b=1/phi, lty=2, lwd=1, col="blue")

# Two reference lines are drawn in the scatterplot, one with slope phi and the other with slope 1/phi to reflect that predictions in the partial autocorrelation computation go in both time directions.

# Next, let us pick a time-series value, x(t), and calculate the relevant fitted values for the lag-2 partial autocorrelation using the formulas
# 
# x-hat(t+2) = phi*x(t+1)
# x-hat(t) = phi*x(t+1)

t <- sample.int(n=n-2, size=1)
xhat2 <- phi*ar.ts[t+1]
xhat0 <- phi*ar.ts[t+1]

# Note that the two fitted value have exactly the same formula in this case. I have given them separate names to remind us that the two formulas would be different in more complex case.

# The following code produces a plot that illustrate the fitted value calculations. The plots zoom in to the region around the relevant values.

vmax <- max(ar.ts[t:(t+2)], xhat2, xhat0)
vmin <- min(ar.ts[t:(t+2)], xhat2, xhat0)
vrng <- vmax - vmin
vlim <- c(vmin-0.1*vrng, vmax+0.1*vrng)
plot(0, 0, type="l", lty=1, lwd=1, xlim=vlim, ylim=vlim, col="white", xlab="previous", ylab="current", main="")
points(ar.ts[1:(n-1)], ar.ts[2:n], pch=16, cex=1)
points(ar.ts[t+1], ar.ts[t+2], pch=16, cex=1, col="black")	
points(ar.ts[t+1], xhat2, pch=16, cex=1, col="blue")
lines(ar.ts[t+1]*c(1,1), c(ar.ts[t+2], xhat2), lty=1, lwd=1, col="blue")
abline(a=0, b=phi, lty=2, lwd=1, col="blue")
points(ar.ts[t], ar.ts[t+1], pch=16, cex=1, col="black")	
points(xhat0, ar.ts[t+1], pch=16, cex=1, col="blue")	
lines(c(ar.ts[t], xhat0), ar.ts[t+1]*c(1,1), lty=1, lwd=1, col="blue")
abline(a=0, b=1/phi, lty=2, lwd=1, col="blue")

# The vertical line connects the values x(t+2) and x-hat(t+2), relative to the prediction formula x-hat(t+2)=phi*x(t+1). The horizontal line connects the values x(t) and x-hat(t), relative to the prediction formula x-hat(t)=phi*x(t+1). To reflect that, in this plot, current values are on the vertical axis and previous values are on the horizontal axis, the latter prediction formula is depiceted as x(t+1)=(1/phi)*x-hat(t)

# ---------------------------------------------------------------------------------
# Calculation in AR(p) models
# ---------------------------------------------------------------------------------

# Partial autocorrelations may be efficiently calculated using the Durbin-Levinson algorithm, but to do so requires that we know the ACF of the time series. The sample PACF would be calculated from the sample ACF, which is a statistic that is easily calcuated from time-series data. To calculate the PACF of a time series model, we must first work out the model's theoretical ACF.

# We will need the user-defined function for implementing the Durbin-Levinson algorithm.

durbin.levinson.calc <- function(n, acf.seq) {
	phi.mat <- matrix(data=0, nrow=n, ncol=n)
	phi.mat[1,1] <- acf.seq[2]
	for (n0 in 2:n) {
		sum1 <- sum(phi.mat[n0-1, 1:(n0-1)]*acf.seq[seq(from=n0, to=2, by=-1)])
		sum2 <- sum(phi.mat[n0-1, 1:(n0-1)]*acf.seq[seq(from=2, to=n0, by=1)])
		phi.mat[n0, n0] <- (acf.seq[n0+1] - sum1) / (1- sum2)
		for (k in 1:(n0-1)) {
			phi.mat[n0, k] <- phi.mat[n0-1, k] - phi.mat[n0, n0]*phi.mat[n0-1, n0-k]
		}
	}
	result <- list(coeff=phi.mat[n,], pacf=diag(phi.mat), matrix=phi.mat)
	return(result)
}

# Let us consider the AR(3) time-series model that we worked with earlier. Its paramaters are

p <- 2
phi1 <- 1.41
phi2 <- -0.5
phi3 <- -0.15
sigw <- 8

# Borrowing code from before, the ACF of this time-series model is calculated as follows

phi.mat <- matrix(data=0, nrow=p, ncol=p)
phi.mat[1,1] <- 1
phi.mat[2,1] <- phi1
phi.mat[2,2] <- -(1-phi2)
phi.mat[2,3] <- phi3
phi.mat[3,1] <- phi2
phi.mat[3,2] <- phi1+phi3
phi.mat[3,3] <- -1

c.vect <- matrix(data=0, nrow=p, ncol=1)
c.vect[1,1] <- 1

phi.inv <- solve(phi.mat)
rho.vect <- phi.inv %*% c.vect

max.lag <-2
phi <- c(phi1, phi2)
acf.seq <- c(as.numeric(rho.vect), rep(x=0, times=max.lag-p))
for (lag in p:max.lag) {
	acf.seq[lag+1] <- sum(phi*acf.seq[seq(from=lag, to=lag-p+1, by=-1)])
}
acf.seq

# The PACF values are calculated by the code below

max.lag <- 3
result <- durbin.levinson.calc(n=max.lag, acf.seq)
round(result$pacf, digits=5)

# As we have done before, we can check these values by simulation. This time we calculate the sample PACF of a simulated time series.

ar.ts <- ar.sim(n=10000, phi=phi, sigw=sigw)
pacf.chk <- pacf(ar.ts, plot=FALSE)
round(cbind(result$pacf, pacf.chk$acf[1:max.lag,,]), digits=4)

# Observe that the PACF at lag h=3 is exactly phi3. This is reflects a general property that will be explored in other materials of this unit.

# =================================================================================
# Non-identifiability of moving-average time series
# =================================================================================

# We have learned that moving-average time-series models can be non-identifyable in the sense that two distinct sets of parameters can define time series with the same probabilistic properties. For example the MA(1) model with parameters

sigw.A <- 1
theta.A <- 5

# is probabilistically identical to the MA(1) model with parameters

sigw.B <- 5
theta.B <- 1/5

# We can explore this concept using simulation. The following code simulates time seriesfrom each of the distinct settings above, and generates plots of the sample paths.

n <- 100
ma.ts.A <- ma.sim(n=n, theta=theta.A, sigw=sigw.A)
ma.ts.B <- ma.sim(n=n, theta=theta.B, sigw=sigw.B)

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(ma.ts.A, xlab="time", ylab="value", main=paste("MA(", length(theta.A), ") series A", sep=""))
ts.plot(ma.ts.B, xlab="time", ylab="value", main=paste("MA(", length(theta.B), ") series B", sep=""))
par(mfrow = c(1, 1))

# The patterns in these two sample paths look roughly the same. The same is observed in the corresponding sample standard deviations...

sd(ma.ts.A)
sd(ma.ts.B)

# ... and (apart from sampling variability) in the ACFs and PACFs.

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(ma.ts.A)
pacf(ma.ts.A)
acf(ma.ts.B)
pacf(ma.ts.B)
par(mfrow = c(1, 1))

