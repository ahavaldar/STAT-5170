# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part B of learning unit 5
# =================================================================================

# As with the part A of the learning unit, this set of examples will use several user-defined functions whose definitions are placed in a separate R file. The following code defines these functions.

setwd(ScriptsDirectory)
source("STAT 5170 unit 5 functions.R")

# We also need to define functions from an external package.

#install.packages("expm")
library(expm)

# =================================================================================
# Global temperature data
# =================================================================================

# We continue our exploration of model building on the global temperature data from before. Recall this is a time series of yearly global temperatures deviations from a baseline temperature across 1880-2009. The associated time sequence is standardize to prevent numerical errors.

data(gtemp, package = "astsa")
ts.plot(gtemp, xlab="year", ylab="temperature deviation", main="global temperature")
n <- length(gtemp)
std.time <- seq(from=0, to=1, length.out=n)

# ---------------------------------------------------------------------------------
# De-trending vs. differencing
# ---------------------------------------------------------------------------------

# In this part of our exploration we compare a de-trending verus a differencing approach. For this we redefine a function from before that calculate basic regression summaries.

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

# The comparison is made between the residuals in regression to the time series of differences. To start, let us de-trend the data using simple linear regression, and collect the residuals

x.vect <- as.matrix(gtemp)
Z.mat <- as.matrix(cbind(rep(x=1, times=n), std.time))
colnames(Z.mat) <- c("z0", "z1")
ls.stat <- reg.summary(x.vect, Z.mat)

# Next, we calculate the time series of differences.

gtemp.diff <- gtemp[2:n] - gtemp[1:(n-1)]

# Note how the differences are calculated. It can be convenient to obtain the same difference time series using the built-in function "diff".

gtemp.diff <- diff(gtemp)

# The values are the same

cbind(as.numeric(gtemp[2:n] - gtemp[1:(n-1)]), as.numeric(diff(gtemp)))

# Now that we have the residuals and difference time series, let us compare them graphically, in the manner of a residual plot. From the plot it is clear that it is the time series of differences looks more like a stationary time series.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(std.time, ls.stat$res, xlab="standardized time", type="p", pch=16, cex=1, ylab="residual", main="")
abline(h=0, lty=2, lwd=1, col="blue")
plot(std.time[2:n], gtemp.diff, xlab="standardized time", type="p", pch=16, cex=1, ylab="difference", main="")
abline(h=0, lty=2, lwd=1, col="blue")
par(mfrow = c(1, 1))

# It can be helpful to also compare the autocorrelation functions. Again, we will see that the patterns in the time series of differences look more like those of white noise

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
acf(ls.stat$res)
acf(gtemp.diff)
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Predictive assessment and model checking
# ---------------------------------------------------------------------------------

# The following code uses the user-defined functions defined at the top of this set of examples to calculate assessments of predictive performance and to check whether the model adequately accounts for autocorrelatiosn. The foud models to be compared are listed as follows. 

# Model 1: Regression model with uncorrelated residuals. Associated analysis results are calculated with the user-defined function "run.dsamp.reg.wn".

# Model 2: Regression model with AR(1) residuals. Associated analysis results are calculated with the user-defined function "run.mcmc.reg.ar1".

# Model 3: Time series of differences is white noise. Associated analysis results are calculated with the user-defined function "run.mcmc.diff.wn".

# Model 4: Time series of differences is AR(1). Associated analysis results are calculated with the user-defined function "run.mcmc.diff.ar1".

# Refer to the external file where these functions are defined for explanations of the associated syntax, which is demonstrated below.

# We start with analyses under the de-trending models. 

# ********** Analysis under Model 1: **********

# Specify the model and computational parameters.

model.parms <- list()
model.parms$x.vect <- x.vect
model.parms$Z.mat <- Z.mat
comp.parms <- list()
comp.parms$max.lag <- 25

# Use the "run.dsamp.reg.wn" function to implement the analysis.

# The code below displays the assessments of predictive performance and makes a plot of the posterior predictive p-values. In the plot, the solid black line marks the p-values associated with absolute predicted autocorrelation, the solid black line marks the p-values associated with maximum absolute predicted autocorrelation, and the solid green line marks the p-values associated with sum of absolute predicted autocorrelations.

comp.parms$n.samp <- 10000
result <- run.dsamp.reg.wn(model.parms, comp.parms)
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="posterior predictive p-values")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")
result$stat
result$pval

# We store the results in a new variable for reference later.

result.reg.wn <- result

# ********** Analysis under Model 2: **********

# Define additonal parameters and use the "run.mcmc.reg.ar1" function to implement the analysis. This function uses an MCMC algorithm, so we much burn in the simulation, and subsequently rerun it.

model.parms$phi1 <- 0
comp.parms$n.samp <- 100
comp.parms$phi1.rad <- 0.20
comp.parms$max.kappa <- 1000

# Initialize and burn in
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))

# Implement the main run
comp.parms$n.samp <- 300
result <- run.mcmc.reg.ar1(model.parms, comp.parms)
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="posterior predictive p-values")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")
par(mfrow = c(1, 1))
result$stat
result$pval
result$rate

# Store the results in a new variable for reference later.

result.reg.ar1 <- result

# Let us now move on the the analyses under the differencing models. 

# ********** Analysis under Model 3: **********

# This model involves an initial values of the time series, x0, which is treated as an additional parameter.

# Again, the simulation is carried out using an MCMC algorithm, which must be burnt in an then rerun.

# Initialize and burn in
model.parms <- list()
model.parms$x <- as.numeric(gtemp)
model.parms$x0 <- mean(diff(gtemp))
comp.parms <- list()
comp.parms$n.samp <- 100
comp.parms$m <- 15
comp.parms$max.kappa <- 1000
comp.parms$x0.rad <- 0.5
comp.parms$max.lag <- 25
result <- run.mcmc.diff.wn(model.parms, comp.parms)
model.parms$x0 <- result$samp$x0[comp.parms$n.samp]
plot(1:comp.parms$n.samp, result$samp$x0, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("x0 (acc.rate=", sprintf("%6.4f", result$rate$x0), ")", sep=""))

# Implement the main run
comp.parms$n.samp <- 300
result <- run.mcmc.diff.wn(model.parms, comp.parms)
model.parms$x0 <- result$samp$x0[comp.parms$n.samp]
par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(1:comp.parms$n.samp, result$samp$x0, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("x0 (acc.rate=", sprintf("%6.4f", result$rate$x0), ")", sep=""))
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="posterior predictive p-values")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")
par(mfrow = c(1, 1))
result$stat
result$pval
result$rate
result.diff.wn <- result

# ********** Analysis under Model 4: **********

# The final model is the most complicated one. It includes two parameters, the autogregressive parameter phi1 and the initial time series value x0, whose calculations must be monitored and tuned.

# Initialize and burn in
model.parms <- list()
model.parms$x <- as.numeric(gtemp)
model.parms$x0 <- mean(diff(gtemp))
model.parms$phi1 <- 0
comp.parms <- list()
comp.parms$n.samp <- 100
comp.parms$m <- 15
comp.parms$max.kappa <- 1000
comp.parms$x0.rad <- 0.4
comp.parms$phi1.rad <- 0.3
comp.parms$max.lag <- 25
result <- run.mcmc.diff.ar1(model.parms, comp.parms)
model.parms$x0 <- result$samp$x0[comp.parms$n.samp]
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(1:comp.parms$n.samp, result$samp$x0, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("x0 (acc.rate=", sprintf("%6.4f", result$rate$x0), ")", sep=""))
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))
par(mfrow = c(1, 1))

# Implement the main run
comp.parms$n.samp <- 300
result <- run.mcmc.diff.ar1(model.parms, comp.parms)
model.parms$x0 <- result$samp$x0[comp.parms$n.samp]
model.parms$phi1 <- result$samp$phi1[comp.parms$n.samp]
par(mfrow = c(3, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(1:comp.parms$n.samp, result$samp$x0, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("x0 (acc.rate=", sprintf("%6.4f", result$rate$x0), ")", sep=""))
plot(1:comp.parms$n.samp, result$samp$phi1, type="l", lty=1, lwd=1, xlab="iteration", ylab="value", main=paste("phi1 (acc.rate=", sprintf("%6.4f", result$rate$phi1), ")", sep=""))
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="posterior predictive p-values")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")
par(mfrow = c(1, 1))
result$stat
result$pval
result$rate

# It is striking how well the posterior predictive checks go through in this model.

# Store the results in a new variable for reference later.

result.diff.ar1 <- result

# Results based on longer simulation runs are stored in an external data file.

setwd(DataDirectory)
#save(result.reg.wn, result.reg.ar1, result.diff.wn, result.diff.ar1, file="gtemp MCMC 5B.RData")
load(file="gtemp MCMC 5B.RData") #result.reg.wn, result.reg.ar1, result.diff.wn, result.diff.ar1

# It is interesting to exampine the p-value plots of all four models in a single graphic. Here we see how distinctly well Model 4 accounts for autocorrelations, relative to the other three.

par(mfrow = c(2, 2))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
result <- result.reg.wn
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="de-trend, white noise")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")

result <- result.reg.ar1
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="de-trend, AR(1)")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")

result <- result.diff.wn
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="difference, white noise")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")

result <- result.diff.ar1
plot(c(0, comp.parms$max.lag+1), c(0.5, 0.5), type="l", lty=1, lwd=1, col="white", ylim=c(-0.05, 1.05), xlab="lag", ylab="p-value", main="difference, AR(1)")
abline(h=c(0.05, 0.95), lty=2, lwd=1, col="blue")
lines(1:comp.parms$max.lag, result$pval$indiv, lty=1, lwd=3, col="black")
lines(1:comp.parms$max.lag, result$pval$max, lty=1, lwd=3, col="red")
lines(1:comp.parms$max.lag, result$pval$sum, lty=1, lwd=3, col="green")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Forecasting
# ---------------------------------------------------------------------------------

# The user-defined functions associated with the differencing models additonally produce samples of predicted forecasts. Plots of these forecasts are generated by the following code.

t1 <- 90
qprobs <- c(0.025, 0.5, 0.975)
n.qnt <- length(qprobs)
n <- length(model.parms$x)
par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00))
result <- result.diff.wn
max.m <- dim(result$samp$fore)[1]
qmat <- matrix(nrow=n.qnt, ncol=max.m)
for (m in 1:max.m) {
	qmat[,m] <- as.numeric(quantile(result$samp$fore[m,], probs=qprobs))
}
ymin <- min(model.parms$x, qmat)
ymax <- max(model.parms$x, qmat)
yrng <- ymax-ymin
plot(t1:n, model.parms$x[t1:n], type="l", lty=1, lwd=1, col="black", xlim=c(t1-1,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="previous", ylab="current", main="difference, white noise")
points(t1:n, model.parms$x[t1:n], pch=1, cex=1, col="black")
lines(n:(n+max.m), c(model.parms$x[n],qmat[2,]), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), qmat[2,], pch=1, cex=1, col="red")
lines(n:(n+max.m), c(model.parms$x[n],qmat[1,]), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(model.parms$x[n],qmat[3,]), lty=2, lwd=1, col="blue")
result <- result.diff.ar1
max.m <- dim(result$samp$fore)[1]
qmat <- matrix(nrow=n.qnt, ncol=max.m)
for (m in 1:max.m) {
	qmat[,m] <- as.numeric(quantile(result$samp$fore[m,], probs=qprobs))
}
ymin <- min(model.parms$x, qmat)
ymax <- max(model.parms$x, qmat)
yrng <- ymax-ymin
plot(t1:n, model.parms$x[t1:n], type="l", lty=1, lwd=1, col="black", xlim=c(t1-1,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="previous", ylab="current", main="difference, AR(1)")
points(t1:n, model.parms$x[t1:n], pch=1, cex=1, col="black")
lines(n:(n+max.m), c(model.parms$x[n],qmat[2,]), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), qmat[2,], pch=1, cex=1, col="red")
lines(n:(n+max.m), c(model.parms$x[n],qmat[1,]), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(model.parms$x[n],qmat[3,]), lty=2, lwd=1, col="blue")
par(mfrow = c(1, 1))

# =================================================================================
# Exponentially-weighted moving-average models
# =================================================================================

# Our next set of numerical examples explores the properties of exponentially-weighted moving-average (EWMA) time series, and related calculations. To start, the following code creates a user-defined function that simulates a single sample path up of an EWMA time series.

ewma.sim <- function(n, sigw, lambda, delta=0, x0=0) {
	w0 <- rnorm(n=1, mean=0, sd=sigw)
	wn.samp <- rnorm(n=n, mean=0, sd=sigw)
	ewma.samp <- numeric(length=n)
	ewma.samp[1] <- delta + x0 + wn.samp[1] - lambda*w0
	for (t in 2:n) {
		ewma.samp[t] <- delta + ewma.samp[t-1] + wn.samp[t] - lambda*wn.samp[t-1]
	}
	ewma.ts <- ts(ewma.samp)
	return(ewma.ts)
}

# The following code executes the function and makes a plot of the simulated sample path.

n <- 25
x0 <- 0
sigw <- 1
lambda <- 0.75
delta <- 0.25
ewma.ts <- ewma.sim(n=n, sigw=sigw, lambda=lambda, delta=delta, x0=x0)
ts.plot(ewma.ts, xlab="", ylab="value", main=paste("EWMA: x0=", x0, ", delta=", delta, ", lambda=", lambda, sep=""))

# To get a sense of the properties of an EWMA time series, the following code independently simulates multiple such time series under the same parameters, and overlays them in a plot.

n.rep <- 5
ewma.seq <- matrix(data=0, nrow=n, ncol=n.rep)
for (i.rep in 1:n.rep) {
	ewma.seq[, i.rep] <- as.numeric(ewma.sim(n=n, sigw=sigw, lambda=lambda, delta=delta, x0=x0))
}

ymin <- min(ewma.seq); ymax <- max(ewma.seq); yrng <- ymax-ymin
plot(0, 0, type="l", lty=1, lwd=1, col="white", xlim=c(0, n+1), ylim=c(ymin-0.1*yrng, ymax+0.1*yrng), xlab="", ylab="value", main=paste("EWMA: x0=", x0, ", delta=", delta, ", lambda=", lambda, sep=""))
for (i.rep in 1:n.rep) {
	lines(1:n, ewma.seq[, i.rep], lty=1, lwd=1)
}

# ---------------------------------------------------------------------------------
# Means and variances
# ---------------------------------------------------------------------------------

# Formulas for the mean and variance function of an EWMA time series (xt), with a fixed starting value, are
#
# E[xt] = delta*t + x0
#
# and
#
# Var[xt] = (1 + (t-1)*(1-lambda)^2 + lambda^2)*sigw^2.
#
# where x0 is the initial value at time t=0, delta is a the drift parameter, and sigw is the standard deviation of the underlying white-noise sequence.
#
# The following uses simulation to numerically compute the expectations and variances of (randomly selected) individual measurements.

n.rep <- 10000
t <- sample.int(n=n, size=1, replace=FALSE)
sum1 <- 0
sum2 <- 0
for (i.rep in 1:n.rep) {
	x <- as.numeric(ewma.sim(n=n, sigw=sigw, lambda=lambda, delta=delta, x0=x0))
	sum1 <- sum1 + x[t]
	sum2 <- sum2 + x[t]^2
}
mean.sim <- sum1 / n.rep
var.sim <- sum2/n.rep - (sum1/n.rep)^2

# The means, variances, and autocorrelation given by the above formulas are calculated as follows.

mean.form <- delta*t+x0
var.form <- (1+(t-1)*(1-lambda)^2+lambda^2)*sigw^2

# Comparisons of the corresponding calculations are made by the code below

print(paste("Time value: t=", t, sep=""))
print(paste("E[xt]: simulated=", sprintf("%6.4f", mean.sim), ", formula=", sprintf("%6.4f", mean.form), sep=""))
print(paste("Var[xt]: simulated=", sprintf("%6.4f", var.sim), ", formula=", sprintf("%6.4f", var.form), sep=""))

# ---------------------------------------------------------------------------------
# Forecasting
# ---------------------------------------------------------------------------------

# To calculate an m-step-ahead forecast of an EWMA time series, first evaluate the one-step-ahead recursive formula
#
# x.pred(n+1;n) = (1-lambda)*x(n) + lambda*x.pred(n;n-1) + lambda^n*(delta+x0),
#
# and then
#
# x.pred(n+m;n) = (m-1)*delta + x.pred(n+1;n)
#
# Mean square prediction errors are given by the formula
#
# P(n+m;n) = (1 + (m-1)*(1-lambda)^2)*sigw^2
#

# These formulas are now demonstrated on a simulated EWMA time series.

n <- 25
x0 <- 0
sigw <- 1
lambda <- 0.75
delta <- 0.25
ewma.ts <- ewma.sim(n=n, sigw=sigw, lambda=lambda, delta=delta, x0=x0)
ts.plot(ewma.ts, xlab="time", ylab="value", main=paste("EWMA: x0=", x0, ", delta=", delta, ", lambda=", lambda, sep=""))

# The one-step-ahead forecast is calculated as follows

predx.measr <- numeric(length=n)
predx.measr[1] <- delta + x0
for (t in 2:n) {
	predx.measr[t] <- (1-lambda)*ewma.ts[t-1] + lambda*predx.measr[t-1] + lambda^(t-1)*(delta+x0)
}
onestep <- (1-lambda)*ewma.ts[n] + lambda*predx.measr[n] + lambda^n*(delta+x0)

# Forecasts further ahead are calcuated as

max.m <- 12
predx <- numeric(length=max.m)
predx[1] <- onestep
for (m in 2:max.m) {
	predx[m] = (m-1)*delta + predx[1]
}

# Plots of the simulated time series and predictions are generated by the following code

ymin <- min(ewma.ts,predx)
ymax <- max(ewma.ts,predx)
yrng <- ymax-ymin
plot(1:n, ewma.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="time", ylab="value", main="")
points(1:n, ewma.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(ewma.ts[n], predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")

# Prediction bands are calculated as follows

pred.err <- numeric(length=max.m)
for (m in 1:max.m) {
	pred.err[m] <- (1 + (m-1)*(1-lambda)^2)*sigw^2
}
pred.err

# The code below converts the squared prediction errors to 95% prediction bands and adds them to the plot

alpha <- 0.05
z.cut <- qnorm(p=1-alpha/2, mean=0, sd=1)
lo.bnd <- -z.cut*sqrt(pred.err)
hi.bnd <- z.cut*sqrt(pred.err)

ymin <- min(ewma.ts,predx+lo.bnd)
ymax <- max(ewma.ts,predx+hi.bnd)
yrng <- ymax-ymin
plot(1:n, ewma.ts, type="l", lty=1, lwd=1, col="black", xlim=c(0,n+max.m+1), ylim=c(ymin-0.1*yrng,ymax+0.1*yrng), xlab="time", ylab="value", main="")
points(1:n, ewma.ts, pch=1, cex=1, col="black")
lines(n:(n+max.m), c(ewma.ts[n],predx), lty=1, lwd=1, col="red")
points((n+1):(n+max.m), predx, pch=1, cex=1, col="red")
lines(n:(n+max.m), c(ewma.ts[n],predx+lo.bnd), lty=2, lwd=1, col="blue")
lines(n:(n+max.m), c(ewma.ts[n],predx+hi.bnd), lty=2, lwd=1, col="blue")

# =================================================================================
# Autoregressive integrated moving-average models
# =================================================================================

# Our exploration of ARIMA models will demonstrate concepts and calculations on a time series of GNP values from 1947 to 2002, which is loaded from an external file using the "read.table" function.

setwd(DataDirectory)
gnp96.raw <- read.table("gnp96.dat", header=FALSE, sep="")
year <- gnp96.raw[[1]]
gnp.ts <- gnp96.raw[[2]]

# The data are transformed into a growth-rate time series, which is defined as the differences of the log-transformed GNP time series.

gnp.ts <- ts(data=gnp96.raw[[2]], frequency=4, start=c(1947, 1))
dloggnp.ts <- ts(data=diff(log(gnp96.raw[[2]])), frequency=4, start=c(1947, 2))

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(gnp.ts, xlab="year", ylab="value", main="gross national product")
ts.plot(dloggnp.ts, xlab="year", ylab="growth rate", main="growth rate gross national product")
par(mfrow = c(1, 1))

# We seek an ARIMA model for the log-transformed GNP time series, but calculations are made as if we seek an ARMA model for the growth-rate time series. That is we apply the "arima" function directly to the growth-rate time series having set the difference parameter to d=0. We could instead calculate the log-transformed GNP time series and apply the "arima" function with difference parameter d=1 to that. However, implementing calculations the way that is done here more clearly illustrates concepts, and also simplifies the syntax of the "arima" function when a possibly non-zero mean is specified for the growth-rate time series. 

# ********** ARIMA(1, 1, 0) **********
# (Equivalent to ARMA(1, 0) applied to growth-rate)

p <- 1; d <- 0; q <- 0
dloggnp.tsa <- arima(dloggnp.ts, order=c(p,d,q))
tsdiag(dloggnp.tsa,gof.lag=25)
dloggnp.tsa$coef
sqrt(diag(dloggnp.tsa$var.coef))
sqrt(dloggnp.tsa$sigma2)
dloggnp.tsa$aic

# ********** ARIMA(0, 1, 2) **********
# (Equivalent to ARMA(0, 2) applied to growth-rate)

p <- 0; d <- 0; q <- 2
dloggnp.tsa <- arima(dloggnp.ts, order=c(p,d,q))
tsdiag(dloggnp.tsa,gof.lag=25)
dloggnp.tsa$coef
sqrt(diag(dloggnp.tsa$var.coef))
sqrt(dloggnp.tsa$sigma2)
dloggnp.tsa$aic

# ********** ARIMA(1, 1, 1) **********
# (Equivalent to ARMA(1, 1) applied to growth-rate)

p <- 1; d <- 0; q <- 1
dloggnp.tsa <- arima(dloggnp.ts, order=c(p,d,q))
tsdiag(dloggnp.tsa,gof.lag=25)
dloggnp.tsa$coef
sqrt(diag(dloggnp.tsa$var.coef))
sqrt(dloggnp.tsa$sigma2)
dloggnp.tsa$aic

# ********** ARIMA(1, 1, 2) **********
# (Equivalent to ARMA(1, 2) applied to growth-rate)

p <- 1; d <- 0; q <- 2
dloggnp.tsa <- arima(dloggnp.ts, order=c(p,d,q))
tsdiag(dloggnp.tsa,gof.lag=25)
dloggnp.tsa$coef
sqrt(diag(dloggnp.tsa$var.coef))
sqrt(dloggnp.tsa$sigma2)
dloggnp.tsa$aic

# ********** ARIMA(2, 1, 2) **********
# (Equivalent to ARMA(2, 2) applied to growth-rate)

p <- 2; d <- 0; q <- 2
dloggnp.tsa <- arima(dloggnp.ts, order=c(p,d,q))
tsdiag(dloggnp.tsa,gof.lag=25)
dloggnp.tsa$coef
sqrt(diag(dloggnp.tsa$var.coef))
sqrt(dloggnp.tsa$sigma2)
dloggnp.tsa$aic
