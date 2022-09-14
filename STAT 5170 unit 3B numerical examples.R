# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part B of learning unit 3
# =================================================================================
#
# =================================================================================
# Simulated ARMA time series
# =================================================================================
#
# Adding to the previously defined functions for simulating time series, the following function may be used to simulate a sample path from an ARMA model.

arma.sim <- function(n, phi, theta, sigw, x0=NA, w0=NA) {
	p <- length(phi)
	q <- length(theta)
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	if (is.na(x0)) {
		pastx <- rep(x=0, times=p)
	} else {
		pastx <- c(x0, rep(x=0, times=p-length(x0)))
	}
	if (is.na(w0)) {
		pastw <- rep(x=0, times=q)
	} else {
		pastw <- c(w0, rep(x=0, times=q-length(w0)))
	}
	arma.samp <- numeric(length=n)
	for (t in 1:n) {
		arma.samp[t] <- sum(phi*pastx) + sum(theta*pastw) + wn.samp[t]
		pastx <- c(arma.samp[t], pastx[1:p-1])
		pastw <- c(wn.samp[t], pastw[1:q-1])
	}
	arma.ts <- ts(arma.samp)
	return(arma.ts)
}

# Consider an ARMA(1,1) model with parameters

sigw <- 3.25
phi <- 0.75
theta <- 0.5

# The following code simulates a time series from the is model and plots both its sample path

n <- 100
arma.ts <- arma.sim(n=n, phi=phi, theta=theta, sigw=sigw)
ts.plot(arma.ts, xlab="", ylab="value", main=paste("ARMA(", length(phi), ",", length(theta), ")", sep=""))

# Its ACF and PACF are displayed as follows.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(arma.ts)
pacf(arma.ts)
par(mfrow = c(1, 1))

# =================================================================================
# Infinite MA and AR representations of ARMA time series
# =================================================================================
#
# This next example illustrates the use of general recursive formulas for determining the coefficients of infinite MA or AR representations of an ARMA time series.
# 
# As a simple initial example, consider the ARMA(1,1) model with parameters 

phi <- 0.9
theta <- 0.5

# Suppose our interest is to determine the coefficients, psi(j), of its infinite MA representation.

# The coefficients up to index j=p=1 are determined by the formula
# 
# psi(1) = theta + phi
#
# The coefficients psi(p+k) for k at least 1 are determined by the formula
# 
# psi(k+1) = phi(1)*psi(k) 
# 
# The first several coefficeints are calculated as follows

max.k <- 5
psi <- numeric(length=max.k+1)
psi[1] <- theta + phi
for (k in 1:max.k) {
	psi[k+1] <- phi*psi[k] 
}
psi

# It may be helpful to know that this sequence is also given by the formula
# 
# psi(j) = (theta + phi)*phi^(j-1)
# 
# This is calculated as follows

(theta + phi)*phi^(1:(max.k+1)-1)

# Let us now practice on more a complicated model, an ARMA(3,2) model with parameters

phi1 <- 0.9
phi2 <- 0.22
phi3 <- -0.24

theta1 <- 0.3
theta2 <- -0.4

# The coefficients, psi(j), of its infinite MA representation are determined up to  index j=p=3 by the formula
# 
# psi(1) = theta1 + phi1
# psi(2) = theta2 + phi2 + phi1*psi(1)
# psi(3) = phi3 + phi1*psi(2) + phi2*psi(1)

# The coefficients psi(p+k) for k at least 1 are determined by the formula
# 
# psi(k+3) = phi1*psi(k+2) + phi2*psi(k+1) + phi3*psi(k)
# 
# The first several coefficients are calculated as follows

max.k <- 5
psi <- numeric(length=max.k+3)
psi[1] <- theta1 + phi1
psi[2] <- theta2 + phi2 + phi1*psi[1]
psi[3] <- phi3 + phi1*psi[2] + phi2*psi[1]
for (k in 1:max.k) {
	psi[k+3] <- phi1*psi[k+2] + phi2*psi[k+1] + phi3*psi[k]
}
psi

# As for its infinite AR representation, the coefficients, pi(j), up to index j=q=2 are determined by the formula
# 
# pi(1) = phi1 + theta1
# pi(2) = phi2 + theta2 - theta1*pi(1)

# The coefficient pi(q+k) for k=1 is
# 
# pi(3) = phi3 - theta1*pi(2) - theta2*pi(1)
# 
# The coefficients pi(q+k) for k at least 2 are determined by the formula
# 
# pi(k+2) = -theta1*pi(k+1) - theta2*pi(k)
# 
# The first several coefficeints are calcualted as follows

max.k <- 5
pi <- numeric(length=max.k+2)
pi[1] <- phi1 + theta1
pi[2] <- phi2 + theta2 - theta1*pi[1]
pi[3] <- phi3 - theta1*pi[2] - theta2*pi[1]
for (k in 2:max.k) {
	pi[k+2] <- -theta1*pi[k+1] - theta2*pi[k]
}
pi

# =================================================================================
# AR(p) forecasting
# =================================================================================
#
# The next several examples will demonstrate the use of m-step-ahead prediction formulas, and related formulas, on several time series models. We start with forecasting under an autoregressive model

# ---------------------------------------------------------------------------------
# Predictions
# ---------------------------------------------------------------------------------

# In a previous example, we worked with an AR(p) time-series model with p=3 and parameters

p <- 3
sigw <- 5
phi1 <- 0.6
phi2 <- 0.25
phi3 <- -0.15

# The following code simulates a sample path from this time-series model and generates a plot. We use a user-defined function defined previously to simulate the sample path.

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

n <- 100
phi <- c(phi1, phi2, phi3)
ar.ts <- ar.sim(n=n, phi=phi, sigw=sigw)
ts.plot(ar.ts, xlab="", ylab="value", main=paste("AR(", length(phi), ")", sep=""))

# Suppose our goal is to produce m-step-ahead forecasts, x.pred(n+m), starting with the last time point of the simulated time series. This is easily done using the recursive formula
#
# x.pred(n+m) = phi1*x.pred(n+m-1) + phi2*x.pred(n+m-2) + phi3*x.pred(n+m-3)
#
# where x.pred(t) = x(t) if t <= n.
#
# Several such m-step-ahead forecasts are calculated as follows

max.m <- 12
p <- length(phi)
predx <- numeric(length=max.m)
pastx <- ar.ts[seq(from=n, to=n-p+1)]
for (m in 1:max.m) {
	predx[m] = phi1*pastx[1] + phi2*pastx[2] + phi3*pastx[3]
	pastx <- c(predx[m], pastx[1:p-1])
}

# The simualated time series and predictions are plotted by the following code

plot(1:n, ar.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), xlab="previous", ylab="current", main="")
points(1:n, ar.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(ar.ts[n],predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")
abline(h=0, lty=2, lwd=1, col="blue")

# Observe how the forecasts decay to the mean of the time series.

# ---------------------------------------------------------------------------------
# Prediction bands
# ---------------------------------------------------------------------------------
#
# To add prediction bands to the plot, we need the first several coefficients of the model's infinite moving-average representation. These are calculated as follows.

max.k <- max.m-4
psi <- numeric(length=max.k+3)
psi[1] <- phi1
psi[2] <- phi2 + phi1*psi[1]
psi[3] <- phi3 + phi1*psi[2] + phi2*psi[1]
for (k in 1:max.k) {
	psi[k+3] <- phi1*psi[k+2] + phi2*psi[k+1] + phi3*psi[k]
}
psi

# With these values, the first several squared prediction errors are calculated below. 

pred.err <- numeric(length=max.m)
pred.err[1] <- sigw^2
for (m in 2:max.m) {
	pred.err[m] <- pred.err[m-1] + sigw^2*psi[m-1]^2
}
pred.err

# The following code converts the squared prediction errors to 95% prediction bands and adds them to the plot

alpha <- 0.05
z.cut <- qnorm(p=1-alpha/2, mean=0, sd=1)
lo.bnd <- -z.cut*sqrt(pred.err)
hi.bnd <- z.cut*sqrt(pred.err)

ymin <- min(ar.ts,predx+lo.bnd)
ymax <- max(ar.ts,predx+hi.bnd)
yrng <- ymax-ymin
plot(1:n, ar.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="previous", ylab="current", main="")
points(1:n, ar.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(ar.ts[n],predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")
lines(n:(n+max.m), c(ar.ts[n],predx+lo.bnd), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(ar.ts[n],predx+hi.bnd), lty=2, lwd=1, col="blue")

# ---------------------------------------------------------------------------------
# Innovations
# ---------------------------------------------------------------------------------

# The following code calculates the innovations on the simulated time series.

pred.samp <- numeric(length=n-p-1)
inno.samp <- numeric(length=n-p-1)
for (t in (p+1):n) {
	pred.val <- phi1*ar.ts[t-1] + phi2*ar.ts[t-2] + phi3*ar.ts[t-3]
	inno.val <- ar.ts[t] - pred.val
	pred.samp[t-p-1] = pred.val
	inno.samp[t-p-1] = inno.val
}

# The innvations are white noise.
inno.ts <- ts(inno.samp)
ts.plot(inno.ts, xlab="", ylab="value", main="innovations")

# The ACF and PACF provide an additional check:

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(inno.ts)
pacf(inno.ts)
par(mfrow = c(1, 1))

# Note also that the sample standard deviation of the innovations is close to the specified white-noise standard deviation:

c(sigw, sd(inno.samp))

# Suppose instead that the model is misspecified as an AR(1) with autoregressive parameter

phi1 <- 0.1

# The innovations are calculated under this model according to 

pred.samp <- numeric(length=n-2)
inno.samp <- numeric(length=n-2)
for (t in 2:n) {
	pred.val <- phi1*ar.ts[t-1]
	inno.val <- ar.ts[t] - pred.val
	pred.samp[t-2] = pred.val
	inno.samp[t-2] = inno.val
}

# These innovations exhibit patterns that deviate from what we would expect if they were white noise. It is hard to see those patterns just by plotting the values...

inno.ts <- ts(inno.samp)
ts.plot(inno.ts, xlab="", ylab="value", main="innovations")

# ...but the ACF and PACF are helpful diagnostic tools

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(inno.ts)
pacf(inno.ts)
par(mfrow = c(1, 1))

# =================================================================================
# ARMA(p,q) forecasting
# =================================================================================

# Let is consider forecasting under the ARMA(3,2) model we worked with earlier, with the following parameters.

phi1 <- 0.9
phi2 <- 0.22
phi3 <- -0.24

theta1 <- 0.3
theta2 <- -0.4

# Suppose further that the white-noise standard deviation is

sigw <- 1.75

# The following code simulates a sample path from this time-series model and generates a plot. 

n <- 100
phi <- c(phi1, phi2, phi3)
theta <- c(theta1, theta2)
arma.ts <- arma.sim(n=n, phi=phi, theta=theta, sigw=sigw)
ts.plot(arma.ts, xlab="", ylab="value", main=paste("ARMA(", length(phi), ",", length(theta), ")", sep=""))

# ---------------------------------------------------------------------------------
# Predictions
# ---------------------------------------------------------------------------------

# A prediction formula for this model is deduced from a truncated version of its infinite AR representation. As in a previous example before, the first several coefficients of this representation are calculated as follows.

max.k <- 20
pi <- numeric(length=max.k+2)
pi[1] <- phi1 + theta1
pi[2] <- phi2 + theta2 - theta1*pi[1]
pi[3] <- phi3 - theta1*pi[2] - theta2*pi[1]
for (k in 2:max.k) {
	pi[k+2] <- -theta1*pi[k+1] - theta2*pi[k]
}
pi

# The coefficients with subscripts above about ten are very near zero, and so just these first j=10 coefficients are used in the prediction formula.

# Several such m-step-ahead forecasts, starting at the end of the simulated time series, are calculated as follows

max.m <- 25
max.j <- 10
predx <- numeric(length=max.m)
pastx <- arma.ts[seq(from=n, to=n-max.j+1)]
for (m in 1:max.m) {
	predx[m] = sum(pi[1:max.j]* pastx)
	pastx <- c(predx[m], pastx[1:max.j-1])
}

# The simualated time series and predictions are plotted by the following code

plot(1:n, arma.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), xlab="previous", ylab="current", main="")
points(1:n, arma.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(arma.ts[n],predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")
abline(h=0, lty=2, lwd=1, col="blue")

# ---------------------------------------------------------------------------------
# Prediction bands
# ---------------------------------------------------------------------------------
#
# Prediction errors, for calcuating prediction bands, are calculated from the model's infinite moving-average representation. As we did before, the first several coefficients of this representation are calculated as follows. 

max.k <- max.m-4
psi <- numeric(length=max.k+3)
psi[1] <- theta1 + phi1
psi[2] <- theta2 + phi2 + phi1*psi[1]
psi[3] <- phi3 + phi1*psi[2] + phi2*psi[1]
for (k in 1:max.k) {
	psi[k+3] <- phi1*psi[k+2] + phi2*psi[k+1] + phi3*psi[k]
}

# Notice that these are calculated from the initial ARMA(3,2) model, not from the truncated infinite-AR representation.

# Having calculated the relevant coefficients, the first several squared prediction errors are calculated below. 

pred.err <- numeric(length=max.m)
pred.err[1] <- sigw^2
for (m in 2:max.m) {
	pred.err[m] <- pred.err[m-1] + sigw^2*psi[m-1]^2
}
pred.err

# 95% prediction bands are added to the plot by the following code.

alpha <- 0.05
z.cut <- qnorm(p=1-alpha/2, mean=0, sd=1)
lo.bnd <- -z.cut*sqrt(pred.err)
hi.bnd <- z.cut*sqrt(pred.err)

ymin <- min(arma.ts,predx+lo.bnd)
ymax <- max(arma.ts,predx+hi.bnd)
yrng <- ymax-ymin
plot(1:n, arma.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="previous", ylab="current", main="")
points(1:n, arma.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(arma.ts[n],predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")
lines(n:(n+max.m), c(arma.ts[n],predx+lo.bnd), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(arma.ts[n],predx+hi.bnd), lty=2, lwd=1, col="blue")

# The enormous size of the prediction bands partly reflects the relatively slow decay of the infinite-AR coefficients.

psi

# =================================================================================
# MA(q) forecasting
# =================================================================================

# In a previous example we worked with the MA(4) time-series model with the following parameters 

sigw <- 1
theta1 <- 0.1
theta2 <- 0.5
theta3 <- -0.8
theta4 <- -0.4

# A simulated single sample path and plot is generated as follows. We use a previously-defined function for simulation.

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

n <- 100
theta <- c(theta1, theta2, theta3, theta4)
ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw, w0=NA)
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))

# We used (a close variation of) the following code to calculate the autocovariance function of this time series, across h = 0, 1, ..., q. 

q <- length(theta)
q.ext <- q+1
theta.ext <- c(1, theta)
acov <- numeric(length=q.ext)
for (lag in 0:q) {
	theta.lag <- c(theta.ext[(lag+1):q.ext], rep(x=0, times=lag))
	acov[lag+1] <- sigw^2*sum(theta.lag*theta.ext)
}
acov

# In what follows, we use the Innovations Algorithm for forcasting. As preparation for coding that algorithm, we append zeros to the variable storing the autocovarance values in order to expand indexing

max.m <- 12
acov <- c(acov, rep(x=0, times=n+max.m))

# The Innovations Algorithm is coded below. The the first step of the algorithm is to calculate the innovations within the existing data. Notice that it simulataneously calculates predictions and squared prediction errors.

pred.seq <- numeric(length=n+max.m)
pred.err <- numeric(length=n+max.m)
theta.mat <- matrix(data=0, nrow=n+max.m-1, ncol=n+max.m-1)
pred.seq[1] <- 0
pred.err[1] <- acov[1]
t <- 1
theta.mat[1,1] <- acov[2]/pred.err[1]
pred.err[2] <- acov[1] - theta.mat[1,1]^2*pred.err[1]
pred.seq[2] <- theta.mat[1,1]*(ma.ts[1] - pred.seq[1])
for (t in 2:n) {
	tp1 <- t+1
	theta.mat[t,t] <- acov[tp1]/acov[1]
	for (j in 1:(t-1)) {
		theta.mat[t,t-j] <- acov[tp1-j]
		for (k in 0:(j-1)) {
			theta.mat[t,t-j] <- theta.mat[t,t-j] - theta.mat[j,j-k]*theta.mat[t,t-k]*pred.err[k+1]
		}
		theta.mat[t,t-j] <- theta.mat[t,t-j] / pred.err[j+1]
	}
	pred.err[tp1] <- acov[1]
	for (j in 0:(t-1)) {
		pred.err[tp1] <- pred.err[tp1] - theta.mat[t,t-j]^2*pred.err[j+1]
	}
	pred.seq[tp1] <- 0
	for (j in 1:t) {
		pred.seq[tp1] <- pred.seq[tp1] + theta.mat[t,j]*(ma.ts[t-j+1] - pred.seq[tp1-j+1])
	}
}

# The next step is to calculate the predictions beyond the the existing data. The only change in the code below, compared to that of the first step, is the indexing of main loop, and that of the predictions and prediction errors within the loop. 

for (t in (n+1):(n+max.m-1)) {
	tp1 <- t+1
	theta.mat[t,t] <- acov[tp1]/acov[1]
	for (j in 1:(t-1)) {
		theta.mat[t,t-j] <- acov[tp1-j]
		for (k in 0:(j-1)) {
			theta.mat[t,t-j] <- theta.mat[t,t-j] - theta.mat[j,j-k]*theta.mat[t,t-k]*pred.err[k+1]
		}
		theta.mat[t,t-j] <- theta.mat[t,t-j] / pred.err[j+1]
	}
	pred.err[tp1] <- acov[1]
	for (j in 1:(t-1)) {
		pred.err[tp1] <- pred.err[tp1] - theta.mat[t,j]^2*pred.err[t-j+1]
	}
	pred.seq[tp1] <- 0
	for (j in (t-n+1):t) {
		pred.seq[tp1] <- pred.seq[tp1] + theta.mat[t,j]*(ma.ts[t-j+1] - pred.seq[tp1-j+1])
	}
}

# The following code generates a plot of the simualated time series, predictions, and 95% prediction bands.

alpha <- 0.05
z.cut <- qnorm(p=1-alpha/2, mean=0, sd=1)
lo.bnd <- -z.cut*sqrt(pred.err[(n+1):(n+max.m)])
hi.bnd <- z.cut*sqrt(pred.err[(n+1):(n+max.m)])
predx <- pred.seq[(n+1):(n+max.m)]

ymin <- min(ma.ts,predx+lo.bnd)
ymax <- max(ma.ts,predx+hi.bnd)
yrng <- ymax-ymin
plot(1:n, ma.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="previous", ylab="current", main="")
points(1:n, ma.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(ma.ts[n],predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")
lines(n:(n+max.m), c(ma.ts[n],predx+lo.bnd), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(ma.ts[n],predx+hi.bnd), lty=2, lwd=1, col="blue")

# Several observations about this example are worth mentioning:
# 
# First, notice that the prediction at t=n+m is zero whenever m > q. This directly reflects the structure of the MA model, which specifies that m > q implies that x(t+m) would be uncorrelated with x(t), and so the prediction would be at the mean of the time series, which is zero.

# Upon closely examining the algorithm's calculations, we see that the vast majority of the of the entries in theta.mat are zero. We only find non-zero entries in the first q=4 columns.

theta.mat[1:10, 1:q]

theta.mat[1:10, (q+1):(q+10)]

# This, too, is a reflection of the MA model's autocovariance function, gamma(h), which is zero whenever h > q.

# Together, these properties suggest that the Innovations Algorithm may be customized to improve computational costs when forecasting under an MA mode. In other cases, it remains a suitable algorithm for general use. 
