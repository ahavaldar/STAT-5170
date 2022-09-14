# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part B of learning unit 2
# =================================================================================
#
# =================================================================================
# Sampling distribution of the SACF 
# =================================================================================
#
# Some of the properties of the sample autocorrelation function under hypothetical repeated sampling are covered in the notes. In the present example these properties are explored.

# In what follows we use R's built-in function "acf" to calculate autocorrelation values.

# ---------------------------------------------------------------------------------
# Central limit theorem under white noise 
# ---------------------------------------------------------------------------------

# When the SACF is calculated on a white-noise time series, the statistic is known to satisfy a central limit theorem, by which its sampling distribution at lag h is approximately normal with mean zero and standard deviation 1/sqrt(n).

# The following code checks this property by simulated repeated sampling. It examines the sampling distribution of the SACF for non-zero lag values up to 4. Notice that the lengh of each simulated white-noise time series is set initially to n = 5, which specifies a short time series, relative to what is typically encountered in practice. 

n <- 5
max.lag <- 4
n.rep <- 10000
sigw <- 1
ac.seq <- matrix(data=0, nrow=max.lag, ncol=n.rep)
for (i.rep in 1:n.rep) {
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	acf.stat <- acf(wn.samp, plot=FALSE)
	ac.seq[, i.rep] <- as.numeric(acf.stat$acf)[2:(max.lag+1)]
}

# Historgrams of the sumulated sampling distributions are generated as follows. Reference values at 0, -2/sqrt(n) and 2/sqrt(n) are indicated

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
for (lag in 1:4) {
	hist(ac.seq[lag,], xlim=c(-1,1), main=paste("lag=", lag, sep=""))
	abline(v=c(0, -2/sqrt(n), 2/sqrt(n)), lty=1, lwd=1, col="blue")
}
par(mfrow = c(1, 1))

# The sampling distributions appear quite different across the four lags, and do not look particularly bell shaped. However, this changes when the time-series length, n, is set to a larger value. The following code repeats the simulation with time-series length set to n=10.  

n <- 10
max.lag <- 4
n.rep <- 10000
sigw <- 1
ac.seq <- matrix(data=0, nrow=max.lag, ncol=n.rep)
for (i.rep in 1:n.rep) {
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	acf.stat <- acf(wn.samp, plot=FALSE)
	ac.seq[, i.rep] <- as.numeric(acf.stat$acf)[2:(max.lag+1)]
}

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
for (lag in 1:4) {
	hist(ac.seq[lag,], xlim=c(-1,1), main=paste("lag=", lag, sep=""))
	abline(v=c(0, -2/sqrt(n), 2/sqrt(n)), lty=1, lwd=1, col="blue")
}
par(mfrow = c(1, 1))

# Now all of the sampling distributons look about the same, and look bell-shaped, centered at the value zero.

# If the time series length is increased again, to n=25, then the sampling distributions become more conventrated around zero, indicating better precision of estimation.

n <- 25
max.lag <- 4
n.rep <- 10000
sigw <- 1
ac.seq <- matrix(data=0, nrow=max.lag, ncol=n.rep)
for (i.rep in 1:n.rep) {
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	acf.stat <- acf(wn.samp, plot=FALSE)
	ac.seq[, i.rep] <- as.numeric(acf.stat$acf)[2:(max.lag+1)]
}

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
for (lag in 1:4) {
	hist(ac.seq[lag,], xlim=c(-1,1), main=paste("lag=", lag, sep=""))
	abline(v=c(0, -2/sqrt(n), 2/sqrt(n)), lty=1, lwd=1, col="blue")
}
par(mfrow = c(1, 1))

# Try out this simulation on your own using other values of n.

# ---------------------------------------------------------------------------------
# Sampling distribution under an autocorrelated time series
# ---------------------------------------------------------------------------------

# Let us now modify the previous exploration of the sampling distribtuon of the SACF by replacing the white-noise time series with an MA(2) time series.

# For this we use the function defined previously for simulating a moving average time series.

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

# Let us specify the following MA(2) parameters

sigw <- 1
theta <- c(0.9, -0.5)

# The autocorrelations of the time series are known from a previous ACF formula. For the parameter values, above, the numbers are calculated as follows.

q <- length(theta)
q.ext <- q+1
theta.ext <- c(1, theta)
acov.form <- numeric(length=q.ext)
ac.form <- numeric(length=q)
for (lag in 0:q) {
	theta.lag <- c(theta.ext[(lag+1):q.ext], rep(x=0, times=lag))
	acov.form[lag+1] <- sigw^2*sum(theta.lag*theta.ext)
	if (lag > 0) {
		ac.form[lag] <- acov.form[lag+1]/acov.form[1]
	}
}
ac.form

# This gives the autocorrelations up to lag h=2; all higher-lag autocorrelations are zero.

# The simulated sampling distributions are obtained as follows, with time-series length set to n=15.  

n <- 15
max.lag <- 4
n.rep <- 10000
sigw <- 1
ac.seq <- matrix(data=0, nrow=max.lag, ncol=n.rep)
for (i.rep in 1:n.rep) {
	ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw)
	acf.stat <- acf(ma.ts, plot=FALSE)
	ac.seq[, i.rep] <- as.numeric(acf.stat$acf)[2:(max.lag+1)]
}

# Simulated means of each simulated sampling distribution are calculated as follows.

ac.sim <- numeric(length=max.lag)
for (lag in 1: max.lag) {
	ac.sim[lag] <- mean(ac.seq[lag,])
}
ac.sim

# The corresponding histograms are generated by the code below. Here, simulated and formula-derived autocorrelation value are displayed at the top of the relevant plot; the mean of each simulated sampling distribution is indicated by a vertical line.

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
for (lag in 1:4) {
	hist(ac.seq[lag,], xlim=c(-1,1), main=paste("lag=", lag, "; sim=", sprintf("%4.2f", ac.sim[lag]), ", form=", sprintf("%4.2f", ac.form[lag]), sep=""))
	abline(v=ac.sim[lag], lty=1, lwd=1, col="blue")
}
par(mfrow = c(1, 1))

# Here we see that the sampling distribtions are slightly asymmetric when the underlying autocorrelation is not zero, but still centered around that autocorrelation value. They are bell-shaped and centered at zero when the autocorrelation is zero. 

# Next, the simulation is repeated for an AR(1) time series. We do not yet have a formula for the autocorrelation function of an AR(1) time series, so our exploration relies on simulated values.

# Recall the function defined previously for simulating an autoregressive time series.

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

# Let us specify the following AR(1) parameters

sigw <- 1
phi <- 0.9

# The simulation and display of results are implemented as follows.

n <- 15
max.lag <- 4
n.rep <- 10000
sigw <- 1
ac.seq <- matrix(data=0, nrow=max.lag, ncol=n.rep)
for (i.rep in 1:n.rep) {
	ar.ts <- ar.sim(n=n, phi=phi, sigw=sigw)
	acf.stat <- acf(ar.ts, plot=FALSE)
	ac.seq[, i.rep] <- as.numeric(acf.stat$acf)[2:(max.lag+1)]
}

# The simulated means are

ac.sim <- numeric(length=max.lag)
for (lag in 1: max.lag) {
	ac.sim[lag] <- mean(ac.seq[lag,])
}
ac.sim

# The histograms are

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
for (lag in 1:4) {
	hist(ac.seq[lag,], xlim=c(-1,1), main=paste("lag=", lag, "; sim=", sprintf("%4.2f", ac.sim[lag]), sep=""))
	abline(v=ac.sim[lag], lty=1, lwd=1, col="blue")
}
par(mfrow = c(1, 1))

# Here, the asymmetries at non-zero autocorrelation values are seen more strongly.

# =================================================================================
# The SACF and SPACF on example data sets
# =================================================================================

# The next few sections of the numerical example takes use through simple displays of the sample ACF and PACF for a variety of time series. The SACF is calculated using R's built-in function "acf", and the SACF is calculated using its built-in function "pacf".

# ---------------------------------------------------------------------------------
# Sample ACF of simulated AR(1) time series
# ---------------------------------------------------------------------------------

# For a simulated AR(1) time series, with paramaters specified below, a plot of a single sample path is shown below.

n <- 500
sigw <- 1
phi <- 0.9
ar.ts <- ar.sim(n=n, phi=phi, sigw=sigw)
ts.plot(ar.ts, xlab="", ylab="value", main=paste("AR(", length(phi), ")", sep=""))

# Its sample ACF and PACF are displayed as follows.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(ar.ts)
pacf(ar.ts)
par(mfrow = c(1, 1))

# Note here the suggestion of a dampening ACF and PACF that becomes zero at lags higher than p=1.

# ---------------------------------------------------------------------------------
# Sample ACF of a simulated MA(2) time series
# ---------------------------------------------------------------------------------

# For a simulated MA(2) time series, with paramaters specified below, the plot of a single sample path is.

n <- 500
sigw <- 1
theta <- c(0.9, -0.5)
ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw)
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))

# Its sample ACF and PACF are

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(ma.ts)
pacf(ma.ts)
par(mfrow = c(1, 1))

# Here, there is the suggestion of a dampening PACF (with alternating values) and an ACF that becomes zero at lags higher than q=2.

# ---------------------------------------------------------------------------------
# Average monthly temperatures in Dubuque, IA
# ---------------------------------------------------------------------------------

# Recall the Dubuque, IA time series exhibits a strongly seasonal pattern, as can be seen in the plot below.

data(data=tempdub, package="TSA")
ts.plot(tempdub, xlab="year", ylab="temp (F)", main="Average monthly temperatures in Dubuque, IA")

# A similar periodic pattern appears in the sample ACF and PACF.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(tempdub)
pacf(tempdub)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Marriages in the Church of England
# ---------------------------------------------------------------------------------

# Yule's "Marriages in the Church of England" data not stationary, so we would not expect its sample ACF and PACF to have a straightforward interpretation. Recall the data have a decreasing trend:

setwd(DataDirectory)
marriages.raw <- scan(file="marriages.txt", sep="\n")
marriages.ts <- ts(marriages.raw, frequency=1, start=1866)
ts.plot(marriages.ts, xlab="year", ylab="percent", main="Marriages in the Church of England")

# Despite the time series's nonstationarity, the SACF and SPACF statistics can be calculated, as is done below. Non-stationarity is expressed through its sample ACF as a very slowly decaying pattern.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(marriages.ts)
pacf(marriages.ts)
par(mfrow = c(1, 1))

# The associated sample PAC suggests a single non-zero value at lag 1. Again, we do not have any helpful interpretation of these patterns.

# =================================================================================
# Linear prediction and partial autocorrelation
# =================================================================================

# In this section we explore the linear prediction concept and connect it to partial autocorrelation.

# ---------------------------------------------------------------------------------
# Fish population recruitment data
# ---------------------------------------------------------------------------------

# To start, let us introduce a data set with which we have not yet worked. These data are part of the "astsa" package, and are stored in the time-series object named "rec". A plot of the data is as follows.

data(data=rec, package="astsa")
ts.plot(rec, xlab="year", ylab="count", main="Fish population recruitment")

# In supporting materials to the "astsa" package, these data are described as the "recruitment (index of the number of new fish) for a period of 453 months ranging over the years 1950-1987. Recruitment is loosely defined as an indicator of new members of a population to the first life stage at which natural mortality stabilizes near adult levels."

# Plots of its sample ACF and PACF are

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(rec)
pacf(rec)
par(mfrow = c(1, 1))

# The values of the sample PACF at is first two lags are produced by the following code.

pacf.stat <- pacf(rec, plot=FALSE)$acf
pacf.stat[1:2]

# We will next use ideas from linear prediction to recalculate these values.

# ---------------------------------------------------------------------------------
# Linear prediction
# ---------------------------------------------------------------------------------

# Suppose we are interested in the process that gave rise to the fish recruitment data, and are in a scenario where we have collected n measurments, x1, ..., xn, and with these measurements we wish to predict the recruitment count at time n+m. For concreteness, let us say that n=10 and m=3, and the measurements are
# HW QUESTION 9,10
x0 <- c(84.7, 84.5, 85.9, 86.2, 87.8)
n <- length(x0)
m <- 2

# Because the basic formulas for linear prediction are developed for a mean-zero time series, we subtract the sample mean of time series from these values, then make our prediction, and then add the sample mean back at the end.

xbar <- 85
x0.shift <- x0 - xbar 

# We develop a linear prediction formula using the property that its coefficients are determined by the underlying ACF of the time series. We do not actually know the underlying ACF, but we can instead base our work on the sample ACF of the existing data set.
rec <- c(1,0.42, 0.43, 0.35, 0.41, 0.32, 0.29, 0.24, 0.20, 0.21)

acf.seq <- acf(rec, plot=FALSE)$acf

# One way to approach this problem is to set up the system of linear equations that determine the coefficients and then solve. Recall that this system is defined by the equations 
#
# gamma(m+k-1) = phi_n1*gamma(k-1) + ... + phi_nn*gamma(k-n)
#
# for k = 1, ..., n. It is set up for the present problem in matrix vector format as follows.

gam.vect <- matrix(data=NA, nrow=n, ncol=1)
gam.mat <- matrix(data=NA, nrow=n, ncol=n)
for (k in 1:n) {
	gam.vect[k,] <- acf.seq[m+k]
	for (h in 1:n) {
		gam.mat[k,h] <- acf.seq[abs(k-h)+1]
	}
}

# The solution is calculated by matrix inversion 

gam.inv <- solve(gam.mat)
phi.vect <- gam.inv %*% gam.vect
phi.vect

# Subsequently, the prediction is calculated as follows

x.vect <- matrix(data=x0.shift, nrow=n, ncol=1)
pred <- t(x.vect) %*% phi.vect
xbar + as.numeric(pred)

# The following code assembles this calculation into a user-defined function, and creates a plot that displays predictions at several values of m.

blpred <- function(m, x0, acf.seq, xbar) {
	n <- length(x0)
	gam.vect <- matrix(data=NA, nrow=n, ncol=1)
	gam.mat <- matrix(data=NA, nrow=n, ncol=n)
	for (k in 1:n) {
		gam.vect[k,] <- acf.seq[m+k]
		for (h in 1:n) {
			gam.mat[k,h] <- acf.seq[abs(k-h)+1]
		}
	}
	gam.inv <- solve(gam.mat)
	phi.vect <- gam.inv %*% gam.vect
	x.vect <- matrix(data=x0-xbar, nrow=n, ncol=1)
	pred <- xbar + as.numeric(t(x.vect) %*% phi.vect)
	result <- list(pred=pred, coeff=as.numeric(phi.vect))
	return(result)
}

max.m <- 5
plot(1:n, x0, type="l", lty=1, lwd=3, xlim=c(1,n+max.m), xlab="year", ylab="count", main="Predictions")
x.pred <- numeric(length=max.m)
for (m in 1:max.m) {
	result <- blpred(m, x0, acf.seq, xbar)
	x.pred[m] <- result$pred
}
points(n+1:max.m, x.pred, pch=1, cex=1, col="blue")	
lines(n+0:max.m, c(x0[n],x.pred), lty=1, lwd=3, col="blue")	

# Suppose instead that m=1. For this case, the Durbin-Levinson algorithm may be used to find the coeffcients of the prediciton formula. The following code creates a user-defined function that implements the Durbin-Levinson algorithm.

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

# We can compare the coefficients given by the Durbin-Levinson algorithm with those calculated directly by matrix inversion.

result1 <- durbin.levinson.calc(n=10, acf.seq)
result2 <- blpred(m=1, x0, acf.seq, xbar)
cbind(result1$coeff, result2$coeff)

# The double-scripted phi values in in the Durbin-Levinson algorithm are stored in an n x n matrix. The coefficients, above, are stored in the last row of that matrix.

cbind(result1$coeff, result1$matrix[n,])

# We can also compare the pacf values given by the Durbin-Levinson algorithm with those calculated by R's built in "pacf" function

cbind(result1$pacf, pacf(rec)$acf[1:n])

# These values are stored in the diagonal of the matrix of double-scripted phi values.

cbind(result1$pacf, diag(result1$matrix))
