# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part A of learning unit 2
# =================================================================================
#
# =================================================================================
# Examples of time series
# =================================================================================
#
# This section of the numerical example demonstrates how to use R to read in time series data and generate plots, and demonstrates basic summary statistics as well. Several modes of reading in data are explored.

# ---------------------------------------------------------------------------------
# S&P 500 index
# ---------------------------------------------------------------------------------

# A version of the S&P 500 data set is included in an external R package called "astsa". External R packages add specialized functionality to R, are created by independent developers of data analysis (and other) software. They are kept in an online repository on the Comprehensive R Archive Network (CRAN) website, www.r-project.org.

# An an external R package may be installed using the function "install.packages", with the name of the package as its argument, in quotes. The following command installs the "astsa" package.

install.packages("astsa") 

# Once a package is installed, it remains on your machine for use in your current of any future R session. That is, it does not need to be installed again, even if you quit R. (You may need to reinstall the package if you remove R from your machine.)

# An external R package must still be loaded into any R session in which you want to use its functionality. This can be done using the "library" command, with the name of the package as its argument, WITHOUT quotes. The following command loads the "astsa" package into the current R session.

library(astsa)

# Even if a package is already installed on your machine, it must still be loaded into every new R session. That is, if you quit R, you must load the package again in the next session that you want access to the package.

# The "data" command may be used to generate a list of all data sets that are included with a package. Use the "package=" option to specify the relevant package, in quotes.

data(package="astsa")

# After submitting this command, a new widow pops up, displaying quite a large list of data sets, with short descriptions. Scrolling down, we can see a data set called "sp500w" with the description "Weekly Growth Rate of the Standard and Poor's 500". The following command returns the length of the time series:

n <- length(sp500w)
n

# From this output, we see there is about 10 years worth of weekly data. Moreover, documentation of the "astsa" package indicates that this time series is, specifically, the "Weekly closing returns of the SP 500 from 2003 to September, 2012". By reformating the data as a time series object, we can insert this information into the time series itself:

sp500wr.ts <- ts(data=sp500w, frequency=52, start=c(2003, 1))

# We can examine the attributes of the formatted time series, without filling up the screen by printing all of the data, using the "start", "end", and "frequency" commands:

start(sp500wr.ts)
end(sp500wr.ts)
frequency(sp500wr.ts)

# Thus, we see that the last measurement is for the 41st week of 2012.

# A plot of the time series is generated as follows:

ts.plot(sp500wr.ts, xlab="year", ylab="return", main="S&P 500 weekly returns")

# Observe the substantial distrubtion that appears in late 2008.

# The formula for a return is r(t) = [x(t) - x(t-1)]/x(t-1), which implies x(t) = [1+r(t)]*x(t-1). Thus, if we are able to specify a starting value, x(t0), for the last week of 2002, then we can reconstruct the weekly values of the S&P 500. For purposes of illustation, let us set the starting value to x(t0)=50 dollars. The time sequence of S&P 500 is then calculated as follows:


x0 <- 50
sp500wv <- numeric(length=n)
sp500wv[1] <- (1 + sp500w[1])*x0
for (t in 2:n) {
  sp500wv[t] <- (1 + sp500w[t])*sp500wv[t-1]
}
sp500wv.ts <- ts(data=sp500wv, frequency=52, start=c(2003, 1))

# A plot of the reconstructed time series is generated as follows

ts.plot(sp500wv.ts, xlab="year", ylab="value", main="S&P 500 weekly values")

# ---------------------------------------------------------------------------------
# Monthly beach-shop souvenir sales in Australia
# ---------------------------------------------------------------------------------

# A time series of souvenir-sales data in available in an 
#external file, formatted as comma-separated values (CSV). 

# The data may now be read in using the "read.csv" command:

souvenir.df <- read.csv(file="souvenir.csv", header=FALSE)
names(souvenir.df) <- "sales"
n <- length(souvenir.df$sales)

# The format of these data is that of a "data frame", 
#which is an object that is used in many statistical functions in R. A data frame typically consists of several columns of variables. The "souvenir.df" data frame just created consists of just one variable. The "names" command is used to name this variable.

# These are monthly data, whose first measurement is for 
#January, 1987. This information is incorporated into the data set 
#by converting to a time-series object:

souvenir.ts <- ts(data=souvenir.df$sales/1000, frequency=12, start=c(1987, 1))

# The units of data are also converted to $1000.

# A plot of the time series is generated as follows

ts.plot(souvenir.ts, xlab="year", ylab="sales ($1000)", main="Souvenir sales in Australia")

# An ordinary regression of these data against time may 
#be implemented using R's build in "lm" command, 
#which is used for the analysis of "linear models." 
#The command takes as input a data frame object that stores 
#the data to be analyzed. The following code adds a regressor 
#variable to "sales.df" data frame, while converting the units of 
#sales to $1000, and assigns names to the variables.

souvenir.df <- data.frame(souvenir.df$sales/1000, 1:n)
names(souvenir.df) <- c("sales", "time")

# The regression analysis is implemented as follows

souvenir.lm <- lm(sales ~ time, data=souvenir.df)
summary(souvenir.lm)

# A scatterplot with a fitted regression line, and residual plot 
# is generated as follows.

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(souvenir.df$time, souvenir.df$sales, pch=1, cex=1, xlab="", ylab="")
abline(souvenir.lm, lty=1, lwd=1, col="blue")
plot(souvenir.df$time, resid(souvenir.lm), pch=1, cex=1, xlab="", ylab="")
abline(h=0, lty=1, lwd=1, col="blue")
par(mfrow = c(1, 1))

# The patterns exhibited in these plots clearly indicate departures from the usual assumptions of regression.

# ---------------------------------------------------------------------------------
# Average summer temperatures in Munich
# ---------------------------------------------------------------------------------

# Data for average summmer temperatures in Munich from 1781 to 1988 are stored as an external R data file. They may be read in using the "load" command, as follows.

setwd(DataDirectory)
load(file="summer.RData")

# The data are stored in a data.frame named "summer.df". We can view the name of all variables in the data frame using the "names" command.

names(summer.df)

# Here we see that the data frame consists of just one variable, named "temp". 

# The dimensions of the data frame are returned by the "dim" command. The following code displays the dimentsions and stores the data frame's number rows in a variable.

dim(summer.df)
n <- dim(summer.df)[1]

# The following converts the data to a time-series object, while incorporating time information, and creates a plot.

summer.ts <- ts(data=summer.df$temp, frequency=1, start=1781)
ts.plot(summer.ts, xlab="year", ylab="temp (C)", main="Average summer temperatures in Munich")

# R includes a function for calculating a statistic called the sample autocorrelation function. The function automatically returns a plot of sample autocorrelation values by lag.

acf(summer.ts)

# The values calcuated by this statistic, and graphed in the plot, may be viewed by assigning the output of the function to a variable, and then viewing the stored values.

summer.acf <- acf(summer.ts)
summer.acf

# A close approximatuion to the acf values may be produced manually by calculating the sample correlation of the time series with its time-lagged values. This is done with the following code.

max.lag <- 23
acf.diy <- numeric(length=max.lag+1)
for (lag in 0:max.lag) {
  values <- summer.ts[(1+lag):n]
  lagged <- summer.ts[1:(n-lag)]
  acf.diy[lag+1] <- cor(values,lagged)
}

# The values do not exactly match up, but we at least have gained a sense of the sample acf statistic.

cbind(summer.acf$acf, acf.diy)

# We will carefully examine this statistic in part B of the learning unit.

# R data files can store much more then data: they may also store analysis results, functions definitions, and anything you create and store in a variable during an R session. To create an R data file, use the "save" command, and list the variables you want to store in before declaring the file name. The following stores the data and the output values of the autocorrelation function in a new external R data file file.

setwd(DataDirectory)
save(summer.df, summer.acf, file="summer_add.RData")

# ---------------------------------------------------------------------------------
# Average monthly temperatures in Dubuque, IA
# ---------------------------------------------------------------------------------

# Monthly temperature data for Dubuque, IA are provided in an external package called "TSA". Provided the package is installed on your machine, the data may be access without loading the entire package into the R session. This is accomplished using the following code'

install.packages("TSA") 
data(data=tempdub, package="TSA")
tempdub

# These data are stored as a time-series object, as can be checked using the "is.ts" command:

is.ts(tempdub)

# A plot of the data is generated as follows.

ts.plot(tempdub, xlab="year", ylab="temp (F)", main="Average monthly temperatures in Dubuque, IA")

# ---------------------------------------------------------------------------------
# Marriages in the Church of England
# ---------------------------------------------------------------------------------

# Yule's "Marriages in the Church of England" data are stored in an external data file in text format. The data in this file may be read using the "scan" command.

setwd(DataDirectory)
marriages.raw <- scan(file="marriages.txt", sep="\n")

# It is converted to a time-series object and plotted by the following code.

marriages.ts <- ts(marriages.raw, frequency=1, start=1866)
ts.plot(marriages.ts, xlab="year", ylab="percent", main="Marriages in the Church of England")

# ---------------------------------------------------------------------------------
# IBM stock
# ---------------------------------------------------------------------------------

# Daily proces for IBM stock are also stored in text format.

ibmv.raw <- scan(file="dailyibm.RData", sep="\n")
n <- length(ibmv.raw)

# The following code converts the stock-price data to a time-series object, creates a separate data set of stock returns, and generates a plot of both.

ibmv.ts <- ts(ibmv.raw, frequency=365, start=c(1980, 1))

ibmr.raw <- numeric(length=n)
for (t in 2:n) {
  ibmr.raw[t] <- (ibmv.raw[t] - ibmv.raw[t-1])/ibmv.raw[t-1]
}
ibmr.ts <- ts(ibmr.raw, frequency=365, start=c(1980, 2))

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(ibmv.ts, xlab="year", ylab="price", main="IBM stock prices")
ts.plot(ibmr.ts, xlab="year", ylab="return", main="IBM stock returns")
par(mfrow = c(1, 1))

# =================================================================================
# Elementary time series concepts
# =================================================================================
#
# The next few numerical examples illustrate several elementary time series concepts.

# ---------------------------------------------------------------------------------
# Random walk time series
# ---------------------------------------------------------------------------------

# Two key properties of a random walk is that the expectations, 
#variances, and autocovariances of individual measurements can 
#depend on time. The relevant formulas are 
#E[xt]=delta*t+x0, Var[xt]=t*sigw^2, and Cov[xs,xt]=min(s,t)*sigw^2, 
#where x0 is the initial value at time t=0, delta is a the drift 
#parameter, and sigw is the standard deviation of the underlying 
#white-noise sequence.

# In this example, we explore these properties by simulation. 
#To start, the following code creates a user-defined function 
#that simualtes a single sample path up to time t=25 of a random 
#walk.

rw.sim <- function(n, sigw, delta=0, x0=0) {
  wn.samp <- rnorm(n=n, mean=0, sd=sigw)
  rw.samp <- numeric(length=n)
  rw.samp[1] <- delta + x0 + wn.samp[1]
  for (t in 2:n) {
    rw.samp[t] <- delta + rw.samp[t-1] + wn.samp[t]
  }
  rw.ts <- ts(rw.samp)
  return(rw.ts)
}

# The following code executes the function and makes a 
#plot of the simulated sample path.

n <- 25
x0 <- 0
sigw <- 1
delta <- 0.25
rw.ts <- rw.sim(n=n, sigw=sigw, delta=delta, x0=x0)
ts.plot(rw.ts, xlab="", ylab="value", main=paste("Random walk: x0=", x0, ", delta=", delta, sep=""))

# It is difficult to get a sense of the properties of random walk 
#time series from just a single simulated series. The following 
#code independently simulates multiple random walks under the same 
#parameters.

n.rep <- 5
rw.seq <- matrix(data=0, nrow=n, ncol=n.rep)
for (i.rep in 1:n.rep) {
  rw.seq[, i.rep] <- as.numeric(rw.sim(n=n, sigw=sigw, delta=delta, x0=x0))
}

# The following code genarates a plot the overlaid simulated 
#time series.

ymin <- min(rw.seq); ymax <- max(rw.seq); yrng <- ymax-ymin
plot(0, 0, type="l", lty=1, lwd=1, col="white", xlim=c(0, n+1), ylim=c(ymin-0.1*yrng, ymax+0.1*yrng), xlab="", ylab="value", main=paste("Random walks: x0=", x0, ", delta=", delta, sep=""))
for (i.rep in 1:n.rep) {
  lines(1:n, rw.seq[, i.rep], lty=1, lwd=1)
}

# This gives a better sense of the random walk's behavior.

# The next portion of the example uses simulation to numerically 
#compute the expectations, variances, and autocovariance between 
#two (randomly selected) individual measurements.

n.rep <- 10000
rand.select <- sample.int(n=n, size=2, replace=FALSE)
s <- rand.select[1]
t <- rand.select[2]
sum1 <- c(0,0)
sum2 <- c(0,0)
sum12 <- 0
for (i.rep in 1:n.rep) {
  x <- as.numeric(rw.sim(n=n, sigw=sigw, delta=delta, x0=x0))
  sum1 <- sum1 + c(x[s], x[t])
  sum2 <- sum2 + c(x[s]^2, x[t]^2)
  sum12 <- sum12 + x[s]*x[t]
}
mean.sim <- sum1 / n.rep
var.sim <- sum2/n.rep - (sum1/n.rep)^2
acov.sim <- sum12/n.rep - (sum1[1]/n.rep)*(sum1[2]/n.rep)

# The formulas used in this simulation are motivates by the law of 
#large numbers, which establishes that the mean of a sequence of 
#independently simulated values is an accurate apprimation to the 
#expectation of the values being simulated. Simulating a larger 
#sequence improves the accuracy of the approximation.

# The means, variances, and autocorrelation given by the formulas 
#are calculated as follows.
# for probems 1 and 2
mean.form <- delta*c(s,t)+x0
var.form <- c(s,t)*sigw^2
acov.form <- min(s,t)*sigw^2

# The corresponding calculations are compared as follows.

print(paste("Time values: s=", s, ", t=",t, sep=""))
print(paste("E[xs]: simulated=", sprintf("%6.4f", mean.sim[1]), ", formula=", sprintf("%6.4f", mean.form[1]), sep=""))
print(paste("E[xt]: simulated=", sprintf("%6.4f", mean.sim[2]), ", formula=", sprintf("%6.4f", mean.form[2]), sep=""))
print(paste("Var[xs]: simulated=", sprintf("%6.4f", var.sim[1]), ", formula=", sprintf("%6.4f", var.form[1]), sep=""))
print(paste("Var[xt]: simulated=", sprintf("%6.4f", var.sim[2]), ", formula=", sprintf("%6.4f", var.form[2]), sep=""))
print(paste("Cov[xs, xt]: simulated=", sprintf("%6.4f", acov.sim), ", formula=", sprintf("%6.4f", acov.form), sep=""))

# ---------------------------------------------------------------------------------
# Autoregressive and moving average time series
# ---------------------------------------------------------------------------------

# The following two functions simulate a sample path of either a autoregressive or moving average time series.

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

# The following code executes the functions and makes separate plot of the simulated sample paths.

n <- 100
sigw <- 1
phi <- 0.7
x0 <- 0
theta <- c(0.9, -0.5)
w0 <- 0

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ar.ts <- ar.sim(n=n, phi=phi, sigw=sigw, x0=x0)
ts.plot(ar.ts, xlab="", ylab="value", main=paste("AR(", length(phi), ")", sep=""))
ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw, w0=w0)
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))
par(mfrow = c(1, 1))

# The code genarates plots of multiple simulated sample paths of the time series, each set overlaid within the same plot.

n.rep <- 5
ar.seq <- matrix(data=0, nrow=n, ncol=n.rep)
for (i.rep in 1:n.rep) {
  ar.seq[, i.rep] <- as.numeric(ar.sim(n=n, phi=phi, sigw=sigw, x0=x0))
}
ma.seq <- matrix(data=0, nrow=n, ncol=n.rep)
for (i.rep in 1:n.rep) {
  ma.seq[, i.rep] <- as.numeric(ma.sim(n=n, theta=theta, sigw=sigw, w0=w0))
}

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ymin <- min(ar.seq); ymax <- max(ar.seq); yrng <- ymax-ymin
plot(0, 0, type="l", lty=1, lwd=1, col="white", xlim=c(0, n+1), ylim=c(ymin-0.1*yrng, ymax+0.1*yrng), xlab="", ylab="value", main=paste("AR(", length(phi), ")", sep=""))
for (i.rep in 1:n.rep) {
  lines(1:n, ar.seq[, i.rep], lty=1, lwd=1)
}
ymin <- min(ma.seq); ymax <- max(ma.seq); yrng <- ymax-ymin
plot(0, 0, type="l", lty=1, lwd=1, col="white", xlim=c(0, n+1), ylim=c(ymin-0.1*yrng, ymax+0.1*yrng), xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))
for (i.rep in 1:n.rep) {
  lines(1:n, ma.seq[, i.rep], lty=1, lwd=1)
}
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Autocovariance function of a one-sided linear process
# ---------------------------------------------------------------------------------

# An MA(q) time series is a one-sided linear process such that all but the first q coefficients are zero, apart from the coefficient indexed at 0, which is 1. In this example, we use the autocovariance formula for linear processes to calculate autocovariances of the moving-average time series.

# We start by defining the parameters of the MA(q) time series

n <- 25
sigw <- 1
theta <- c(0.1, 0.5, -0.8, -0.4)

# A plot of a single sample path from this MA(q) time series is generated as follows.

ma.ts <- ma.sim(n=n, theta=theta, sigw=sigw, w0=NA)
ts.plot(ma.ts, xlab="", ylab="value", main=paste("MA(", length(theta), ")", sep=""))

# Examining the autocovariance function given in the notes, we see that the covariance between xs and xt depends only on the time distance h = |s-t| between the two points.

# Although the formula is given in terms of an infinite-sum, because only a finite number of terms are non-zero we can nevertheless evaluate the formula. This is done by the code below at time-distance values h = 0, 1, ..., q.

q <- length(theta)
q.ext <- q+1
theta.ext <- c(1, theta)
acov.form <- numeric(length=q.ext)
for (lag in 0:q) {
  theta.lag <- c(theta.ext[(lag+1):q.ext], rep(x=0, times=lag))
  acov.form[lag+1] <- sigw^2*sum(theta.lag*theta.ext)
}
acov.form

# The following code spot-checks these computed values, using simulation to numericaly calculate the autocorrelation Cov[xs, xt], wher s=t+h at random values of h and t.

lag <- sample.int(n=q+1, size=1) - 1
t <- sample.int(n=n-q-2, size=1)+2
n.rep <- 10000
s <- t + lag
sum1 <- c(0,0)
sum2 <- c(0,0)
sum12 <- 0
for (i.rep in 1:n.rep) {
  x <- as.numeric(ma.sim(n=n, theta=theta, sigw=sigw, w0=NA))
  sum1 <- sum1 + c(x[s], x[t])
  sum2 <- sum2 + c(x[s]^2, x[t]^2)
  sum12 <- sum12 + x[s]*x[t]
}
mean.sim <- sum1 / n.rep
var.sim <- sum2/n.rep - (sum1/n.rep)^2
acov.sim <- sum12/n.rep - (sum1[1]/n.rep)*(sum1[2]/n.rep)

# The corresponding calculations are compared as follows.

print(paste("lag=", lag, ", t=",t, sep=""))
print(paste("Cov[xs, xt]: simulated=", sprintf("%6.4f", acov.sim), ", formula=", sprintf("%6.4f", acov.form[lag+1]), sep=""))