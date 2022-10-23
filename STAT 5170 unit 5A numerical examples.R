# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part A of learning unit 5
# =================================================================================

# This set of examples will use several user-defined functions whose definitions are placed in a separate R file. The following code executes the contents of that file, thereby defining the functions.

#setwd(ScriptsDirectory)
source("STAT 5170 unit 5 functions.R")

# Functionality from a external package is also needed.

#install.packages("expm")
library(expm)

# =================================================================================
# Global temperature data
# =================================================================================

# Basic model-building ideas are illustrated on a time series of yearly global temperatures deviations from a baseline temperature across 1880-2009. The associated time sequence is standardize to prevent numerical errors.

data(gtemp, package = "astsa")
ts.plot(gtemp, xlab="year", ylab="temperature deviation", main="global temperature")
n <- length(gtemp)
std.time <- seq(from=0, to=1, length.out=n)

# ---------------------------------------------------------------------------------
# De-trend the data using regression: exploratory analysis
# ---------------------------------------------------------------------------------

# The following code illustrates the use of residual plots and autocorrelation statistics in an exploratory model-building analysis.

# We start by setting up the response vector and matrix of residuals for simple linear regression
x.vect <- as.matrix(gtemp)
Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time))
colnames(Z.mat) <- c("z0", "z1")

# The following code defines function that calculates basic regression summaries, which are repeatedly calculated below. We have defined a function like this before.

reg.summary <- function(x.vect, Z.mat) {
	n <- dim(Z.mat)[1]
	r <- dim(Z.mat)[2]
	ZpZ <- t(Z.mat) %*% Z.mat
	ZpX <- t(Z.mat) %*% x.vect
	ZpZ.inv <- solve(ZpZ)
	b.hat <- ZpZ.inv %*% ZpX
	res <- as.numeric(x.vect - Z.mat %*% b.hat)
	SSE <- sum(res^2)
	result <- list(b.hat=b.hat, res=res, SSE=SSE)
	return(result)
}

# The basic regression statistics are calculated by the following code; a plot of the fitted regression line is produced.

ls.stat <- reg.summary(x.vect, Z.mat)
b.hat <- ls.stat$b.hat

plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature")
lines(std.time, Z.mat %*% b.hat, lty=1, lwd=3, col="blue")

# This next set of code implements a residual analysis, in the form of a residual plot and plot of the estimated residuals' sample autocorrelation function.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="")
abline(h=0, lty=2, lwd=1, col="blue")
acf(ls.stat$res, main="")
par(mfrow = c(1, 1))

# The analysis is now repeated, but with a column that reflects a quadratic term in the mean function added to the matrix of regressors.

Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time, std.time ^2))
colnames(Z.mat) <- c("z0", "z1", "z2")
ls.stat <- reg.summary(x.vect, Z.mat)
b.hat <- ls.stat$b.hat

plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature")
lines(std.time, Z.mat %*% b.hat, lty=1, lwd=3, col="blue")

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="")
abline(h=0, lty=2, lwd=1, col="blue")
acf(ls.stat$res)
par(mfrow = c(1, 1))

# The analysis is repeated again, with the addition of a cubic term.

Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time, std.time ^2, std.time ^3))
colnames(Z.mat) <- c("z0", "z1", "z2", "z3")
ls.stat <- reg.summary(x.vect, Z.mat)
b.hat <- ls.stat$b.hat

plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature")
lines(std.time, Z.mat %*% b.hat, lty=1, lwd=3, col="blue")

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="")
abline(h=0, lty=2, lwd=1, col="blue")
acf(ls.stat$res)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Incorporate AR(1) autoregressive residual errors
# ---------------------------------------------------------------------------------

# The following example illustrates how the model may be extended by substituting the assumption of AR(1) residuals for the usual assumption that residuals are white noise.

# To see how this works, first obtain the estimated residuals from a standard regression analysis. Start with simple linear reggression.

x.vect <- as.matrix(gtemp)
Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time))
colnames(Z.mat) <- c("z0", "z1")
ls.stat <- reg.summary(x.vect, Z.mat)

# The following code calculates the sample autocorrelation function of the residuals, and from the lag-1 autocorrelation constructs an autocorrelation matrix that is reflective of AR(1) residuals.

p <- 1
acf.seq <- acf(ls.stat$res, plot=FALSE)
phi1 <- acf.seq$acf[2]
max.lag <- 25
acf.seq <- phi1^(0:max.lag)

# The autocorrelation matrix itself constructed using the "bandSparse" function, which is a specialized syntax for specifying banded matrices. The autocorrelations cut off after a maximum number of lags. Corresponding bands exceeding that lag are taken as zero. The output produced by the bandSparse function is a "sparse" matrix, in which no value is actually stored for the bands that are unassigned and taken to be zero. Instead the structure of the matrix is carefully recorded internally, in such a way that takes into account its unassigned entries in subsequent operations. This cuts down on the required amount of internal storage space and total computation time.

bands <- t(matrix(data=acf.seq[1:(max.lag+1)], nrow=max.lag+1, ncol=n))
R.mat <- bandSparse(n=n, k=0:max.lag, diag=bands, symm=TRUE)
R.mat[1:5, 1:5]

# The code that will be used for analysis is collected into user defined functions. The first one calculates an AR(1) autocorrelation function from a set of residuals. 

get.acf.ar1 <- function(res, max.lag) {
	acf.seq <- acf(res, plot=FALSE)
	phi1 <- acf.seq$acf[2]
	acf.seq <- phi1^(0:max.lag)
	return(acf.seq)
}

# The second one implements a modified least-squares analysis under the assumption of residual deviations with the specified autocorrelations.

reg.summary.acorr <- function(x.vect, Z.mat, acf.seq) {
	n <- dim(Z.mat)[1]
	r <- dim(Z.mat)[2]
	max.lagp1 <- length(acf.seq)
	max.lag <- max.lagp1-1
	diags <- t(matrix(data=acf.seq[1:max.lagp1], nrow= max.lagp1, ncol=n))
	R.mat <- bandSparse(n=n, k=0:max.lag, diag=diags, symm=TRUE) # The bandSparse() function is specialized syntax for banded matrices
	kappa.val <- kappa(R.mat) # Check for poor conditioning
	detR <- det(R.mat)
	if (detR <= 0) {
		kappa.val <- Inf
	}
	R.inv <- solve(R.mat)
	ZpRZ <- as.matrix(t(Z.mat) %*% R.inv %*% Z.mat)
	ZpRX <- as.matrix(t(Z.mat) %*% R.inv %*% x.vect)
	ZpRZ.inv <- solve(ZpRZ)
	b.hat <- ZpRZ.inv %*% ZpRX
	res <- x.vect - Z.mat %*% b.hat
	SSE <- as.numeric(t(res) %*% R.inv %*% res)
	result <- list(b.hat=b.hat, res=as.vector(res), SSE=SSE, kappa=kappa.val)
	return(result)
}

# Observe in this function that the usual least-squares formulas have been modified to incorporate the correlation matrix. While the usual least-squares estimator of the regression coefficients is
#
# b.hat = (Z^T*Z)^(-1)*Z^T*x,
#
# here it is 
#
# b.hat = (Z^T*R.inv*Z)^(-1) Z^T*R.inv*x,
#
# where R.inv is the inverse of the correlation matrix R.
#
# Similarly, whereas the usual sum-of-squared errors statistic is
#
# SSE = e1^2 + ... + en^2 = (x - Z*b.hat)^T*(x - Z*b.hat)
#
# having written ei for i'th entry of the vector of estimated residuals, x - Z*b.hat,
#
# here it is 
#
# SSE = (x - Z*b.hat)^T*R.inv*(x - Z*b.hat)
#
# The basic regression statistics and plots are now calcuated under the model formulation that incorporates AR(1) autoregressive residual errors. Plots are generated alongside parallel results generated under the assumption of uncorrelated residuals. Any distinctions are barely discernable in the plots. However, a comparison of AIC statistics, implemented below, indicates that the autoregressive model is substantially worthwhile for modeling relationships among the residual errors. In other words, AR(1) modeling of residuals may not have much of an effect on regression coefficient estimates, but may help to improve the precision of inferences such as confidence and prediction intervals.

ls.stat0 <- reg.summary(x.vect, Z.mat)
acf.seq <- get.acf.ar1(res=ls.stat0$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq)

par(mfrow = c(2, 1))
plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature: AR(1)")
lines(std.time, Z.mat %*% ls.stat$b.hat, lty=1, lwd=3, col="blue")
plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature: uncorr")
lines(std.time, Z.mat %*% ls.stat0$b.hat, lty=1, lwd=3, col="blue")
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="AR(1)")
abline(h=0, lty=2, lwd=1, col="blue")
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat0$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="uncorr")
abline(h=0, lty=2, lwd=1, col="blue")
acf(ls.stat$res)
acf(ls.stat0$res)
par(mfrow = c(1, 1))

# The remainer of the exploratory analysis is now repeated under the new model formulation.

# Add a quadratic term

Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time, std.time ^2))
colnames(Z.mat) <- c("z0", "z1", "z2")
ls.stat0 <- reg.summary(x.vect, Z.mat)
acf.seq <- get.acf.ar1(res=ls.stat0$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq)
b.hat <- ls.stat$b.hat

par(mfrow = c(2, 1))
plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature: AR(1)")
lines(std.time, Z.mat %*% ls.stat$b.hat, lty=1, lwd=3, col="blue")
plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature: uncorr")
lines(std.time, Z.mat %*% ls.stat0$b.hat, lty=1, lwd=3, col="blue")
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="AR(1)")
abline(h=0, lty=2, lwd=1, col="blue")
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat0$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="uncorr")
abline(h=0, lty=2, lwd=1, col="blue")
acf(ls.stat$res)
acf(ls.stat0$res)
par(mfrow = c(1, 1))

# Add a cubic term

Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time, std.time ^2, std.time ^3))
colnames(Z.mat) <- c("z0", "z1", "z2", "z3")
ls.stat0 <- reg.summary(x.vect, Z.mat)
acf.seq <- get.acf.ar1(res=ls.stat0$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq)
b.hat <- ls.stat$b.hat

par(mfrow = c(2, 1))
plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature: AR(1)")
lines(std.time, Z.mat %*% ls.stat$b.hat, lty=1, lwd=3, col="blue")
plot(std.time, gtemp, type="l", lty=1, lwd=1, col="black", xlab="standardized time", ylab="temperature deviation", main="global temperature: uncorr")
lines(std.time, Z.mat %*% ls.stat0$b.hat, lty=1, lwd=3, col="blue")
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="AR(1)")
abline(h=0, lty=2, lwd=1, col="blue")
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat0$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="uncorr")
abline(h=0, lty=2, lwd=1, col="blue")
acf(ls.stat$res)
acf(ls.stat0$res)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Predictive assessment: AIC
# ---------------------------------------------------------------------------------

# The functions defined below extend and replace the previous user-defined functions that calculate basic regression statistics so that they additionally calculate the AIC statistic.

reg.summary <- function(x.vect, Z.mat) {
	n <- dim(Z.mat)[1]
	r <- dim(Z.mat)[2]
	ZpZ <- t(Z.mat) %*% Z.mat
	ZpX <- t(Z.mat) %*% x.vect
	ZpZ.inv <- solve(ZpZ)
	b.hat <- ZpZ.inv %*% ZpX
	mu <- Z.mat %*% b.hat
	res <- as.numeric(x.vect - mu)
	SSE <- sum(res^2)
	AIC <- n*log(2*pi*SSE/n) + n + 2*(r+1)
	result <- list(b.hat=b.hat, mu=mu, res=res, SSE=SSE, AIC=AIC)
	return(result)
}

reg.summary.acorr <- function(x.vect, Z.mat, acf.seq, p.acf) {
	n <- dim(Z.mat)[1]
	r <- dim(Z.mat)[2]
	max.lagp1 <- length(acf.seq)
	max.lag <- max.lagp1-1
	diags <- t(matrix(data=acf.seq[1:max.lagp1], nrow= max.lagp1, ncol=n))
	R.mat <- bandSparse(n=n, k=0:max.lag, diag=diags, symm=TRUE) # The bandSparse() function is specialized syntax for banded matrices
	kappa.val <- kappa(R.mat) # Check for poor conditioning
	detR <- det(R.mat)
	if (detR <= 0) {
		kappa.val <- Inf
	}
	logdetR <- log(abs(detR))
	R.inv <- solve(R.mat)
	ZpRZ <- as.matrix(t(Z.mat) %*% R.inv %*% Z.mat)
	ZpRX <- as.matrix(t(Z.mat) %*% R.inv %*% x.vect)
	ZpRZ.inv <- solve(ZpRZ)
	b.hat <- ZpRZ.inv %*% ZpRX
	mu <- Z.mat %*% b.hat
	res <- x.vect - mu
	SSE <- as.numeric(t(res) %*% R.inv %*% res)
	AIC <- n*log(2*pi*SSE/n) + logdetR + n + 2*(r+p.acf)
	result <- list(b.hat=b.hat, mu=mu, res=as.numeric(res), SSE=SSE, AIC=AIC, kappa=kappa.val)
	return(result)
}

# The models we have been exploring are now cycled through again, this time to calcuate their associated AIC values.

Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time, std.time ^2, std.time ^3))
colnames(Z.mat) <- c("z0", "z1", "z2", "z3")

# Linear term only.
ls.stat <- reg.summary(x.vect, Z.mat[,1:2])
AIC <- ls.stat$AIC
print(paste("white-noise resids: AIC = ", sprintf("%8.3f", AIC), sep=""))
acf.seq <- get.acf.ar1(res=ls.stat$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq, p.acf=2)
AIC <- ls.stat$AIC
print(paste("AR(1) resids: AIC = ", sprintf("%8.3f", AIC), sep=""))

# Linear and quadratic
ls.stat <- reg.summary(x.vect, Z.mat[,1:3])
AIC <- ls.stat$AIC
print(paste("white-noise resids: AIC = ", sprintf("%8.3f", AIC), sep=""))
acf.seq <- get.acf.ar1(res=ls.stat$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq, p.acf=2)
AIC <- ls.stat$AIC
print(paste("AR(1) resids: AIC = ", sprintf("%8.3f", AIC), sep=""))

# Linear, quadratic, and cubic
ls.stat <- reg.summary(x.vect, Z.mat[,1:4])
AIC <- ls.stat$AIC
print(paste("white-noise resids: AIC = ", sprintf("%8.3f", AIC), sep=""))
acf.seq <- get.acf.ar1(res=ls.stat$res, max.lag)
ls.stat <- reg.summary.acorr(x.vect, Z.mat, acf.seq, p.acf=2)
AIC <- ls.stat$AIC
print(paste("AR(1) resids: AIC = ", sprintf("%8.3f", AIC), sep=""))

# ---------------------------------------------------------------------------------
# Implement a Bayesian analysis
# ---------------------------------------------------------------------------------

# The following analysis implements a Bayesian version of the previous analysis, using the user-defined functions that were defined at the start of these numerical examples. Recall that they we defined the user-defined functions by executing code that is stored in an external file.

# The first two functions defined in that file are described as follows.

# run.dsamp.reg.wn: This function is used for de-trending a time series by regression, under the assumption that the residuals are a white-noise time series. The core algorithm used in this function is essentially a copy of the algorithm we explored in the previous learning unit for regression analysis. The sample from the posterior distribution is simulated by direct sampling within a Gibbs strategy. The algorithm itself executes quite quickly.

# run.mcmc.reg.ar1: This function is used for de-trending a time series by regression, under the assumption that the residuals are an AR(1) time series. For a given correlation matrix R, the posterior distribution reflects that obtained under Bayesian regression as usual, except that the least squares statistics that define its mean vector and covariance matrix (of the regression-coefficient posterior distribution), and sum-of-squared error statisitcs (of the white-noise variance posterior distribution) incorporate the correlation matrix in the same manner as in the modified regression summary function defined above. The model does not assume a fixed starting value of the time-series of residual deviations, which distinguishes it from the AR model and computational algorithm we examined in the previous learning unit. The sample from the posterior distribution is simulated using a Gibbs strategy within an MCMC algorithm, using direct sampling to simulate the regression coefficients and white-noise variance, and an Metropolis-Hastings step to simulate the autoregressive parameter. The algorithm requires a burn in phase and monitoring in order to specify its tuning parameters.

# Refer to the external file where these functions are defined for explanations of the associated syntax, which is demonstrated below.

# The models we have been exploring are now cycled through again, wherein the relevant functions for Bayesian analysis are applied to calculated the associated DIC and WAIC values. Skip to line 450, below, to load the results of previously run analyses.

Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time, std.time ^2, std.time ^3))
colnames(Z.mat) <- c("z0", "z1", "z2", "z3")

model.parms <- list()
model.parms$x.vect <- x.vect
comp.parms <- list()
comp.parms$max.lag <- 25
comp.parms$max.kappa <- 1000

# Linear only
model.parms$Z.mat <- Z.mat[,1:2]

# Uncorrelated residuals
comp.parms$n.samp <- 1000
result <- run.dsamp.reg.wn(model.parms, comp.parms)
result.wn1 <- result
result$stat

# AR(1) residuals
# Initialize and burn in
model.parms$phi1 <- 0
comp.parms$n.samp <- 100
comp.parms$phi1.rad <- 0.20
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))

# Implement the main run
comp.parms$n.samp <- 300
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))
result.ar1 <- result
result$stat
result$rate

# Linear and quadratic
model.parms$Z.mat <- Z.mat[,1:3]

# Uncorrelated residuals
comp.parms$n.samp <- 1000
result <- run.dsamp.reg.wn(model.parms, comp.parms)
result.wn2 <- result
result$stat

# AR(1) residuals
# Initialize and burn in
model.parms$phi1 <- 0
comp.parms$n.samp <- 100
comp.parms$phi1.rad <- 0.20
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))

# Implement the main run
comp.parms$n.samp <- 300
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))
result.ar2 <- result
result$stat
result$rate

# Linear, quadratic, and cubic
model.parms$Z.mat <- Z.mat[,1:4]

# Uncorrelated residuals
comp.parms$n.samp <- 1000
result <- run.dsamp.reg.wn(model.parms, comp.parms)
result.wn3 <- result
result$stat

# AR(1) residuals
# Initialize and burn in
model.parms$phi1 <- 0
comp.parms$n.samp <- 100
comp.parms$phi1.rad <- 0.20
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))

# Implement the main run
comp.parms$n.samp <- 300
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))
result.ar3 <- result
result$stat
result$rate

setwd(DataDirectory)
#save(result.wn1, result.ar1, result.wn2, result.ar2, result.wn3, result.ar3, file="gtemp MCMC 5A.RData")
load(file="gtemp MCMC 5A.RData") #result.wn1, result.ar1, result.wn2, result.ar2, result.wn3, result.ar3

# Display results
result.wn <- list()
result.wn$label <- c("   linear", "quadratic", "    cubic", "  quartic")
result.wn$stat <- vector(mode="list", length=3)
result.wn$stat[[1]] <- result.wn1$stat
result.wn$stat[[2]] <- result.wn2$stat
result.wn$stat[[3]] <- result.wn3$stat

result.ar <- list()
result.ar$label <- c("   linear", "quadratic", "    cubic", "  quartic")
result.ar$stat <- vector(mode="list", length=3)
result.ar$stat[[1]] <- result.ar1$stat
result.ar$stat[[2]] <- result.ar2$stat
result.ar$stat[[3]] <- result.ar3$stat

disp.results <- function(result.wn, result.ar) {
	n.reg <- length(result.wn$stat)
	print("          |      DIC      |      WAIC     |")
	print("     mean |     wn  AR(1) |     wn  AR(1) |")
	print("----------+---------------+---------------+")
	for (i.reg in 1:n.reg) {
		disp.line <- paste(result.wn$label[i.reg], " | ", sprintf("%6.1f", result.wn$stat[[i.reg]]$DIC), " ", sprintf("%6.1f", result.ar$stat[[i.reg]]$DIC), " | ", sprintf("%6.1f", result.wn$stat[[i.reg]]$WAIC), " ", sprintf("%6.1f", result.ar$stat[[i.reg]]$WAIC), " |", sep="")
		print(disp.line)
	}
	print("----------+---------------+---------------+")
}
disp.results(result.wn, result.ar)

# =================================================================================
# Seasonal trends
# =================================================================================

# This next example illustates elementary concepts of seasonal variables.

# The following code simulates white-noise time series that is shifted by a single periodic component.

n <- 500
time.seq <- 1:500
omega <- 10/500
mu.seq <- 2*cos(2*pi*omega*time.seq + 0.6*pi)
x.seq <- mu.seq + rnorm(n=n, mean=0, sd=5)

plot(time.seq, x.seq, type="l", lty=2, lwd=1, xlab="time", ylab="x", main="")
lines(time.seq, mu.seq, lty=1, lwd=3, col="blue")

# A scaled periodogram is generated by the code below from the given formulas

scl.pdgram <- numeric(length=n/2)
coeff.seq <- matrix(data=0, nrow=2, ncol=n/2)
freq.seq <- seq(from=0, to=0.5*(n-1)/n, by=1/n)
for (omega in freq.seq) {
	b1.hat <- (2/n)*sum(x.seq*cos(2*pi*omega*time.seq))
	b2.hat <- (2/n)*sum(x.seq*sin(2*pi*omega*time.seq))
	i.col <- omega*n + 1
	coeff.seq[, i.col] <- c(b1.hat, b2.hat)
	scl.pdgram[i.col] <- b1.hat^2 + b2.hat^2
}

plot(freq.seq, scl.pdgram, type="l", lty=1, lwd=1, xlab="frequency", ylab="squared magnitude", main="")
abline(h=0, lty=2, lwd=1)

# The following code recreates portions of the scaled periodogram using regression. It illustrates that the regression coefficients at individual frequencies are unaffected by the presence of other seasonal terms (at different freqencies) in the regression model.

n.freq <- sample.int(n=10, size=1)
freq.samp <- sample(x=freq.seq, size=n.freq, replace=FALSE)
freq.samp[1] <- 10/500
x.vect <- as.matrix(x.seq)
Z.mat <- matrix(data=0, nrow=n, ncol=2*n.freq)
for (i.freq in 1:n.freq) {
	omega <- freq.samp[i.freq]
	Z.mat[, i.freq] <- cos(2*pi*omega*time.seq)
	Z.mat[, n.freq+i.freq] <- sin(2*pi*omega*time.seq)
}
result <- reg.summary(x.vect, Z.mat)
b.hat <- result$b.hat

plot(freq.seq, scl.pdgram, type="l", lty=1, lwd=1, xlab="frequency", ylab="squared magnitude", main="")
abline(h=0, lty=2, lwd=1)
for (i.freq in 1:n.freq) {
	omega <- freq.samp[i.freq]
	b1.hat <- b.hat[i.freq]
	b2.hat <- b.hat[n.freq+i.freq]
	sq.mag <- b1.hat^2 + b2.hat^2
	lines(omega*c(1,1), c(0,sq.mag), lty=1, lwd=3, col="blue")
}

# =================================================================================
# Modeling seasonality: Cardiovascular Mortality in LA
# =================================================================================

# The next example illustates the use of seasonal variable in a time-series analysis. The example data are of weekly measurements of cardiovascular mortality collected in Los Angeles over a ten year period.

data(cmort, package = "astsa")
ts.plot(cmort, xlab="time", ylab="mortality")

# ---------------------------------------------------------------------------------
# Cubic and periodic trends in cardiovascular mortality
# ---------------------------------------------------------------------------------

# The code below sets up a regression analysis that takes into account the possibility of the following trends:

# -- polynomial trends up to order 4
# -- yearly trends, with period 1/52 weeks

n <- length(cmort)
r <- 7
time.seq <- 0.01*(1:n - n/2) # Shift and scale time to avoid numerical overload
omega <- 100*(1/52)
x.vect <- as.matrix(cmort)
Z.mat <- matrix(data=0, nrow=n, ncol=r)
Z.mat[,1] <- rep(x=1, times=n)
Z.mat[,2] <- time.seq
Z.mat[,3] <- time.seq^2
Z.mat[,4] <- time.seq^3
Z.mat[,5] <- time.seq^4
Z.mat[,6] <- cos(2*pi*omega*time.seq)
Z.mat[,7] <- sin(2*pi*omega*time.seq)

# Plots of fitted values that are calculated a part of the basic regression summaries  are displayed below, each calculated according to a regression models that assumes white-noise residual deviations. The fitted values of two models are overlayed.

# Model 1 includes just the polynomial terms, up to quartic, with no seasonal variables

# Model 2 includes both the polynomial terms and the seasonal variables

# Values of AIC produced by these analyses are displayed in the plot.

plot(time.seq, as.numeric(x.vect), type="l", lty=1, lwd=1, xlab="shifted and scaled time", ylab="mortality")
result <- reg.summary(x.vect, Z.mat[,1:5])
lines(time.seq, result$mu, lty=1, lwd=3, col="blue")
text(x=0, y=125, labels=paste("AIC1 =", sprintf("%8.3f", result$AIC)), pos=4)
result <- reg.summary(x.vect, Z.mat)
lines(time.seq, result$mu, lty=1, lwd=3, col="blue")
text(x=0, y=120, labels=paste("AIC2 =", sprintf("%8.3f", result$AIC)), pos=4)

# ---------------------------------------------------------------------------------
# Predictive assessment
# ---------------------------------------------------------------------------------

# In what follows, the above analysis is refined to explore the impact of adding each polynomial terms, with and without seasonal variables. In addition to collecting the basic regression sumaries, a Bayesian analysis is implemented under the assumption of white-noise residual deviations to yield values of DIC and WAIC.

n.reg <- 8
reg.idx <- vector(mode="list", length=n.reg)
i.reg <- 1
reg.idx[[i.reg]] <- c(1, 2)
reg.idx[[i.reg+1]] <- c(reg.idx[[i.reg]], 6, 7)
i.reg <- 3
reg.idx[[i.reg]] <- c(1, 2, 3)
reg.idx[[i.reg+1]] <- c(reg.idx[[i.reg]], 6, 7)
i.reg <- 5
reg.idx[[i.reg]] <- c(1, 2, 3, 4)
reg.idx[[i.reg+1]] <- c(reg.idx[[i.reg]], 6, 7)
i.reg <- 7
reg.idx[[i.reg]] <- c(1, 2, 3, 4, 5)
reg.idx[[i.reg+1]] <- c(reg.idx[[i.reg]], 6, 7)

model.parms <- list()
model.parms$x.vect <- x.vect
comp.parms <- list()
comp.parms$n.samp <- 10000
comp.parms$max.lag <- 25
par(mfrow = c(3, 3))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
for (i.reg in 1:n.reg) {
	plot(time.seq, as.numeric(x.vect), type="l", lty=1, lwd=1, xlab="shifted and scaled time", ylab="mortality")
	result <- reg.summary(x.vect, as.matrix(Z.mat[,reg.idx[[i.reg]]]))
	lines(time.seq, result$mu, lty=1, lwd=3, col="blue")
	text(x=0, y=127, labels=paste("AIC = ", sprintf("%6.1f", result$AIC), sep=""), pos=4)
	model.parms$Z.mat <- as.matrix(Z.mat[,reg.idx[[i.reg]]])
	result <- run.dsamp.reg.wn(model.parms, comp.parms)
	text(x=0, y=119, labels=paste("DIC = ", sprintf("%6.1f", result$stat$DIC), sep=""), pos=4)
	text(x=0, y=111, labels=paste("WAIC = ", sprintf("%6.1f", result$stat$WAIC), sep=""), pos=4)
}
par(mfrow = c(1, 1))
