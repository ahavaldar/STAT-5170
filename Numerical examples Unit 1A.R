# =================================================================================
# STAT 5170: Applied Time Series
# Numerical examples for part A of learning unit 1
# =================================================================================
#
# The aim of this set of examples is to familiarize the student with R and introduce a number of programming techniques that will be used in later examples.
#
# =================================================================================
# Time series objects
# =================================================================================
#
# R includes a specialized format for time series data. As an example consider the following data set of n=96 time-ordered measurements from a hypothetical process.

x.data <- c(85, 88, 94, 70, 58, 63, 59, 47, 49, 41, 46, 50, 64, 65, 79, 79, 93, 89, 103, 93, 88, 82, 74, 61, 63, 58, 53, 48, 45, 38, 38, 52, 45, 56, 67, 78, 66, 61, 72, 79, 94, 69, 88, 77, 73, 55, 64, 80, 71, 79, 99, 84, 89, 91, 80, 67, 60, 62, 73, 66, 72, 73, 81, 78, 77, 83, 73, 86, 73, 53, 58, 56, 62, 68, 76, 66, 70, 79, 82, 82, 86, 103, 105, 91, 83, 82, 85, 87, 80, 82, 81, 71, 81, 91, 96, 102)

# The time series format is implemented by converting this sequence of measurements to a time series object. Its advantage is that is stores information about the starting times, ending times, and intervals assoicated with the sequence.

# Let us suppose the data sequence above is of monthly measurements over an eight year period starting in June 2006. The conversion to a time series object is implemented by the following commands

x.ts <- ts(data=x.data, frequency=12, start=c(2008, 6))

# The data may now quickly be displayed in a table organized for monthly data

x.ts

# Similarly, a nice plot of the data is easily generated using the specialized "ts.plot" function

ts.plot(x.ts, ylab="measurement", main="hypothetical data")

# Quarterly data are also especially compatible with R's time series functionality, as is demonstrated with the following hypothtical data set. 

x.data <- c(29, 34, 34, 42, 42, 44, 44, 29, 27, 27, 26, 22, 28, 20, 15, 8, 16, 11, 12, 11, 15, 23, 27, 33, 32, 31, 22, 11, 21, 16, 11, 12, 16, 14, 16, 15, 11, 10, 8, 11, 17, 16, 20, 22, 30, 25, 33, 34)

x.ts <- ts(data=x.data, frequency=4, start=c(2009, 1))

x.ts

ts.plot(x.ts, ylab="measurement", main="hypothetical data")

# Other time intervals may also be specified, but may not be translated to familiar labels in printed output. For example, daily data organized within weeks is not automatically organized into a table reflecting that format in printed output, nor are the usual names of days of the week generated.

x.data <- c(671, 678, 673, 677, 691, 729, 722, 752, 766, 732, 729, 696, 642, 662, 600, 567, 591, 590, 555, 489, 518, 549, 562, 598, 593, 593, 575, 600, 593, 663, 663, 624, 639, 680, 656, 617, 607, 621, 639, 600, 555, 544, 524, 502, 504, 531, 503, 495, 491, 462, 491, 491, 496, 533, 603, 635, 675, 697, 668, 629, 663, 678, 675, 684, 652, 607, 580, 542, 531, 527, 517, 538, 560, 514, 527, 545, 572, 559, 539, 571, 571, 568, 574, 609, 619, 637, 689, 711, 706, 731, 715, 663, 631, 665, 621, 639, 612, 613, 651, 660, 619, 617, 612, 619, 615)

x.ts <- ts(data=x.data, frequency=7, start=c(4, 3))

x.ts

# However, note that the time-value of the final data point is calculated within the specified organizing system. Also, the week is displayed as the basic unit in a time series plot. 

ts.plot(x.ts, ylab="measurement", main="hypothetical data")

# =================================================================================
# Generating random numbers
# =================================================================================
#
# Often this course we will use simulation to explore the properties of time-series models, or to implement statistical inference procedures. For this, it is helful to refamiliarize ourselves with a few interesting probability distributions, and learn how to work with them in R.

# ---------------------------------------------------------------------------------
# Beta distributions
# ---------------------------------------------------------------------------------

# Random variables that have a beta distribution fall within the unit interval. The graph of a beta density with parameters alpha=1/2 and beta=1/2 is generated as follows.  

alpha <- 1/2
beta <- 1/2
x.grid <- seq(from=0, to=1, length.out=50)
p.grid <- dbeta(x=x.grid, shape1=alpha, shape2=beta)

plot(x.grid, p.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
abline(h=0, lty=1, lwd=1)
abline(v=c(0,1), lty=1, lwd=1)

# The following code uses the "rbeta" function to simulate an independent and identically distributed sample of beta random variables with common parameters alpha=1/2 and beta=1/2.

n.samp <- 1000
x.samp <- numeric(length=n.samp)
for (i.samp in 1:n.samp) {
  x.val <- rbeta(n=1, shape1=alpha, shape2=beta)
  x.samp[i.samp] <- x.val
}

# A histogram of the similated data, displayed with the correspnding density, is generated as follows:

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(x.grid, p.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
abline(v=c(0,1), lty=1, lwd=1)

hist(x.samp, breaks=20, xlab="", ylab="", main="")
par(mfrow = c(1, 1))

# Recall that a uniform distribution is just a beta distribution with parameters alpha=1 and beta=1. A uniform distribution in the unit square is understood as the distribution of pair of independent uniform random variables.

# The following code uses the "runif" function to simulate an independent and identically distributed sample of uniform bivariate random vectors on the square.

n.samp <- 1000
x.samp <- matrix(data=0, nrow=2, ncol=n.samp)
for (i.samp in 1:n.samp) {
  x.vect <- runif(n=2, min=0, max=1)
  x.samp[, i.samp] <- x.vect
}

# A plot of these data is generated as follows:

plot(x.samp[1,], x.samp[2,], pch=16, cex=1, xlab="", ylab="")

# ---------------------------------------------------------------------------------
# Normal distributions and t distributions
# ---------------------------------------------------------------------------------

# Recall that densities of the normal distributions and t distributions are bell shaped, but a t distribution has "thick tails." Overlaid example densities of each of these distributons are generated by the code below. Notice that the normal density is standard normal, and the t density has df=4 degrees of freedom.

x.min <- -3.5
x.max <- 3.5
x.grid <- seq(from=x.min, to=x.max, length.out=50)
z.grid <- dnorm(x=x.grid, mean=0, sd=1)
t.grid <- dt(x=x.grid, df=2)

plot(x.grid, z.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
lines(x.grid, t.grid, lty=2, lwd=3, col="red")
abline(h=0, lty=1, lwd=1)

# Simulated data from each distributions is generated by the following code, and displayed with the corresponding density plots as overlaid histograms.

n.samp <- 1000
z.samp <- numeric(length=n.samp)
t.samp <- numeric(length=n.samp)
for (i.samp in 1:n.samp) {
  z.val <- rnorm(n=1, mean=0, sd=1)
  t.val <- rt(n=1, df=2)
  z.samp[i.samp] <- z.val
  t.samp[i.samp] <- t.val
}

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(x.grid, z.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
lines(x.grid, t.grid, lty=2, lwd=3, col="red")
abline(h=0, lty=1, lwd=1)

hist(t.samp[which((t.samp > x.min) & (t.samp < x.max))], breaks=20, ylim=c(0,225), xlab="", ylab="", main="")
hist(z.samp, breaks=20, col="grey", add=TRUE)
par(mfrow = c(1, 1))

# Consider the distribution of a bivariate normal random vector such that the entries are independent and each has standard normal distribution. Simulated bivariate random vectors from this distribution, when plotted, form a circular shape.

n.samp <- 1000
z.samp <- matrix(data=0, nrow=2, ncol=n.samp)
for (i.samp in 1:n.samp) {
  z.vect <- rnorm(n=2, mean=0, sd=1)
  z.samp[, i.samp] <- z.vect
}

plot(z.samp[1,], z.samp[2,], pch=16, cex=1, xlab="", ylab="")

# Eliptical shapes arise by linear transformation of the random vectors. To simulate a sample of bivariate normal random vectors from a normal distributon with mean mu.vect and covariance matrix Sig.mat, start by simulating pairs of independent standard normal random vectors and apply the following transformation:
#
# x.vect = mu.vect + Sig.sqrt %*% z.vect
#
# where Sig.sqrt is the matrix square-root of Sig.mat, and %*% is R's syntax for matrix multiplication. The matrix square-root, Sig.sqrt, of Sig.mat is defined as the matrix such that
#
# Sig.mat = Sig.sqrt %*% Sig.sqrt
#
# It may be calculated using the "sqrtm" function, which is included in the package called "expm", which must be loaded during the R session.

#install.packages("expm")
library(expm)

# The following code simulates a sample of bivariate normal random vectors from a distributon with respective means mu1=100 and mu2=50, respective standard deviations sig1=30 and sig2=10, and correlation rho=0.75.

mu1 <- 100
mu2 <- 50
sig1 <- 30
sig2 <- 10
rho <- 0.75

mu.vect <- c(mu1, mu2)
Sig.mat <- matrix(data=0, nrow=2, ncol=2)
Sig.mat[1,1] <- sig1^2
Sig.mat[1,2] <- rho*sig1*sig2
Sig.mat[2,1] <- Sig.mat[1,2]
Sig.mat[2,2] <- sig2^2

Sig.sqrt <- sqrtm(Sig.mat)

n.samp <- 1000
z.samp <- matrix(data=0, nrow=2, ncol=n.samp)
x.samp <- matrix(data=0, nrow=2, ncol=n.samp)
for (i.samp in 1:n.samp) {
  z.vect <- rnorm(n=2, mean=0, sd=1)
  z.samp[, i.samp] <- z.vect
  x.samp[, i.samp] <- mu.vect + Sig.sqrt %*% z.vect
}

par(mfrow = c(2, 2))
plot(z.samp[1,], z.samp[2,], pch=16, cex=1, xlab="", ylab="")
plot(x.samp[1,], x.samp[2,], pch=16, cex=1, xlab="", ylab="")
par(mfrow = c(1, 1))

# Note that summary statistics for means, standard deviations, and correlations calculated in the simulated data produce values near to the values specified.

print(paste("mu1=", sprintf("%4.2f", mu1), "; ybar1=", sprintf("%4.2f", mean(x.samp[1,])), sep=""))
print(paste("mu2=", sprintf("%4.2f", mu2), "; ybar2=", sprintf("%4.2f", mean(x.samp[2,])), sep=""))
print(paste("sig1=", sprintf("%4.2f", sig1), "; s1=", sprintf("%4.2f", sd(x.samp[1,])), sep=""))
print(paste("sig2=", sprintf("%4.2f", sig2), "; s2=", sprintf("%4.2f", sd(x.samp[2,])), sep=""))
print(paste("rho=", sprintf("%6.4f", rho), "; r=", sprintf("%6.4f", cor(x.samp[1,], x.samp[2,])), sep=""))

# ---------------------------------------------------------------------------------
# Chi-square distributions
# ---------------------------------------------------------------------------------

# A chi-square distribution is asymmetric with a long tail to the right. The following code produces a plot of a chi-square density with df=3 degrees of freedom together with a histogram of a sample of values simulated from that distribution.

df <- 3
x.min <- 0
x.max <- 8
x.grid <- seq(from=x.min, to=x.max, length.out=50)
p.grid <- dchisq(x=x.grid, df=df)

n.samp <- 1000
x.samp <- numeric(length=n.samp)
for (i.samp in 1:n.samp) {
  x.val <- rchisq(n=1, df=df)
  x.samp[i.samp] <- x.val
}

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(x.grid, p.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
abline(v=0, lty=1, lwd=1)

hist(x.samp[which(x.samp < x.max)], breaks=20, xlab="", ylab="", main="")
par(mfrow = c(1, 1))

# A scaled inverse-chi-square distribution is defined as a transformation of the chi-square distributions. Specifically, if the random variable X has a chi-square distribution with df degrees of freedom, then, for a constant value lambda, the random variable
#
# Y = lambda / X
#
# has a scaled inverse-chi-square distribution with scale parameter lambda and df degrees of freedom.

# The density of a scaled inverse-chi-square distribution can be worked out using the transformation theorem. It is 
#
# p(y) = (0.5*lambda)^(0.5*df)*y^(-(0.5*df+1))*exp(-0.5*lambda/y)/gamma(0.5*df)
#
# However, this density formula is not needed to simulate a sample from the distribution. To simulate such a sample, just simulate values of the chi-square random variable and apply the transformation.

# The following code produces a plot of scaled inverse-chi-square density with scale parameter lambda=5 and df=3 degrees of freedom together with a sample histogram.

lambda <- 5
df <- 3
y.min <- 0
y.max <- 8
y.grid <- seq(from=y.min, to=y.max, length.out=50)
p.grid <- (0.5*lambda)^(0.5*df)*y.grid^(-(0.5*df+1))*exp(-0.5*lambda/y.grid)/gamma(0.5*df)

n.samp <- 1000
y.samp <- numeric(length=n.samp)
for (i.samp in 1:n.samp) {
  x.val <- rchisq(n=1, df=df)
  y.samp[i.samp] <- lambda / x.val
}

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(y.grid, p.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
abline(h=0, lty=1, lwd=1)
abline(v=0, lty=1, lwd=1)

hist(y.samp[which(y.samp < y.max)], breaks=20, xlab="", ylab="", main="")
par(mfrow = c(1, 1))

# ---------------------------------------------------------------------------------
# Triangular distributions
# ---------------------------------------------------------------------------------

# A triangular distribution has a density of the type displayed in the plot produced by the following code.

plot(c(0, 0.5, 2), c(0, 1, 0), type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
abline(h=0, lty=1, lwd=1)
abline(v=0.5, lty=2, lwd=1)

# Observe that this density has a peak of p(x) = 1 located at x=0.5, away from which the density decreases linearly to zero at values x=0 and x=2.

# It may not be immediately clear how one would go about simulating a sample from this distribution. One approach is to apply a general numerical strategy called inverse-CDF sampling. This strategy assumes that an explicit formula for the inverse-CDF of the distribution is available. Recall that the inverse-CDF of a random variable X is the function ICDF(u) such that ICDF(u)=x if, and only if, CDF(x)=u, where CDF(x) = P[x <= X] is the cumulative distribution function of X.

# One can check that the density of the particular triangular distribution, above, has the formula
#
# p(x) = 2*x if 0 <= x < 0.5; and (2/3)*(2-x) if 0.5 <= x < 2.
#
# The corresponding cumulative distribution function has the formula
#
# CDF(x) = x^2 if 0 <= x < 0.5; and (4/3)*x*(1-0.25*x) - 1/3 if 0.5 <= x < 2.
#
# The corresponding inverse-CDF has the formula
#
# ICDF(u) = sqrt(u) if 0 <= u < 0.25; and 2 - sqrt(3*(1-u)) if 0.25 <= x < 1

# These functions are coded as user-defined functions in R as follows

tri.dens <- function(x) {
  if ((0 <= x) & (x < 0.5)) {
    val <- 2*x
  } else if ((0.5 <= x) & (x <= 2)) {
    val <- (2/3)*(2-x)
  } else {
    val <- NA
  }
  return(val)
}

tri.cdf <- function(x) {
  if ((0 <= x) & (x < 0.5)) {
    val <- x^2
  } else if ((0.5 <= x) & (x <= 2)) {
    val <- (4/3)*x*(1-0.25*x) - 1/3
  } else {
    val <- NA
  }
  return(val)
}

tri.invcdf <- function(u) {
  if ((0 <= u) & (u < 0.25)) {
    val <- sqrt(u)
  } else if ((0.25 <= u) & (u <= 1)) {
    val <- 2 - sqrt(3*(1-u))
  } else {
    val <- NA
  }
  return(val)
}

# Plots of the CDF and inverse-CDF are generated by the following code

n.val <- 50
x.grid <- seq(from=0, to=2, length.out=n.val)
u.grid <- seq(from=0, to=1, length.out=n.val)
cdf.grid <- numeric(length=n.val)
invcdf.grid <- numeric(length=n.val)
for (i.val in 1:n.val) {
  x.val <- x.grid[i.val]
  cdf.grid[i.val] <- tri.cdf(x.val)
  u.val <- u.grid[i.val]
  invcdf.grid[i.val] <- tri.invcdf(u.val)
}

par(mfrow = c(2, 2))
plot(x.grid, cdf.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="", main="CDF")
abline(h=c(0,1), lty=1, lwd=1)

plot(u.grid, invcdf.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="", main="inverse-CDF")
abline(v=c(0,1), lty=1, lwd=1)
par(mfrow = c(1, 1))

# To simulate a single random variable from this distribtuon, the first step is to simulate a random variable U from a uniform distribution on the unit interval; the second step is to evaluate X=ICDF(U), the inverse-CDF of the triangular distribution at the value U. It can be shown that the random variable X has the desired triangular distribution.

# The following code produces a histogram of a sample of values simulated from the triangular distribution, alongside a plot of the density.

n.val <- 49
x.grid <- seq(from=0, to=2, length.out=n.val)
p.grid <- numeric(length=n.val)
for (i.val in 1:n.val) {
  x.val <- x.grid[i.val]
  p.grid[i.val] <- tri.dens(x.val)
}

n.samp <- 1000
x.samp <- numeric(length=n.samp)
for (i.samp in 1:n.samp) {
  u.val <- runif(n=1, min=0, max=1)
  x.samp[i.samp] <- tri.invcdf(u.val)
}

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
plot(x.grid, p.grid, type="l", lty=1, lwd=3, col="blue", xlab="", ylab="")
abline(h=0, lty=1, lwd=1)
abline(v=0.5, lty=2, lwd=1)

hist(x.samp, breaks=20, xlab="", ylab="", main="")
par(mfrow = c(1, 1))

# In theory, the inverse-CDF strategy applies to any distribution, however complex. In practice, however, it may not always be possible to implement the strategy within the practical limitations of the context. The difficuly is that the inverse-CDF must be known completely, which means that it must be possible to determine the constant factor that makes the density integrate to one. This can sometimes be accomplished numerically, but in complex problems the inverse-CDF sampling strategy is sometimes infeasible.

# =================================================================================
# Dependent and identically distributed samples
# =================================================================================
#
# Students of statistics quickly become familiar with the concept of an independent and identically distributed sample of values, but when first exposed to concepts in time series analysis it may take some time to grasp the notion of dependent sample measurements that are nevertheless identically distributed.

# The following is a simple, toy example of a simulation that produces a random time series of dependent values.

# The steps to produce the sample are as follows:

# Step 1: Specify an even-numbered sample size, n, and simulate an independent and identically distributed size-n sample from some distribuion.

# Step 2: Find the pair of values that are closest to each other, randomly place them next to each other in the time series. 

# Step 3: Repeat Step 2 for the values that have not yet been placed in the time series, until there are no more values to place.

# The simulation is coded, below, as a user-defined function, taking as input the initial i.i.d. sample.

sim.did.samp <- function(iid.samp) {
  # ----- ----- subfunctions ----- -----
  find.close.pair <- function(in.seq) {
    n.val <- length(in.seq)
    if (n.val <= 2) {
      out.seq <- in.seq
    } else {
      min.dist <- Inf
      for (i.val in 1:(n.val-1)) {
        for (j.val in (i.val+1):n.val) {
          dist <- abs(in.seq[i.val] - in.seq[j.val])
          if (dist < min.dist) {
            min.dist <- dist
            pair.idx <- c(i.val, j.val)
          }
        }
      }
      out.seq <- in.seq[c(pair.idx, setdiff(1:n.val, pair.idx))]
    }
    return(out.seq)
  }
  # ----- ----- --- main --- ----- -----
  n <- length(iid.samp)
  halfn <- n/2
  did.samp <- numeric(length=n)
  pair.ord <- seq(from=1, to=n, by=2)[sample.int(n=halfn, size=halfn, replace=FALSE)]
  rem.samp <- iid.samp
  for (i.pair in 1:(halfn-1)) {
    new.samp <- find.close.pair(rem.samp)
    rem.samp <- new.samp[3:length(rem.samp)]
    did.samp[c(pair.ord[i.pair], pair.ord[i.pair]+1)] <- new.samp[1:2]
  }
  i.pair <- halfn
  did.samp[c(pair.ord[i.pair], pair.ord[i.pair]+1)] <- rem.samp
  return(did.samp)
}

# The following code generates plots of the initial i.i.d. sample and dependent sample. 

n <- 36
iid.samp <- rnorm(n=n, mean=0, sd=1)
iid.ts <- ts(data=iid.samp, frequency=2, start=c(1,1))
did.samp <- sim.did.samp(iid.samp)
did.ts <- ts(data=did.samp, frequency=2, start=c(1,1))

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
ts.plot(iid.ts, ylab="measurement", main="i.i.d. sample")
ts.plot(did.ts, ylab="measurement", main="d.i.d. sample")
par(mfrow = c(1, 1))

# Observe in the plot of the dependent sample that adjacent values tend to be similar, provided the time index of the first value in the pair is odd.

# This pattern indicates there is dependency within the sample values, which is understood through the sample's time ordering. It is an example of a characteristic called autocorrelation, meaning that the time series is correlated with itself.

# A summary statistic measuring this type of autocorrelation is calculated as follows.

odd.ac <- cor(did.samp[seq(from=1, to=n-1, by=2)], did.samp[seq(from=2, to=n, by=2)])
odd.ac

# Notice autocorrelation is much smaller and is often negative within adjacent values where the time index of the first value in the pair is even. In other words, the values within those pairs often tend to be very different.

even.ac <- cor(did.samp[seq(from=2, to=n-2, by=2)], did.samp[seq(from=3, to=n-1, by=2)])
even.ac

# Here, we see that the autocorrelation within adjacent values depends on the starting index of the pair. This is a more complex pattern of autocorrelation then arises in most of the time-series models that we will study in the course. Specifically, much of our attention will be focused on what are called stationary time series, where autocorrelation between value is only determined by the time distance between the values, leaving the time location of the values irrelevant.

# The two calculations of autocorrelation made above can be explored more extensively using simulated repeated sampling. Here, the entire time series is repeatedly re-simulated, and each time the autocorrelation statistics are recalculated on the newly simulated data, and stored for later examination. The resulting sequences of summary-statistic values are simulations of the statistics' sampling distributions, which can be examined to get a sense of the probabilistic properties of the algorithm that simulates this type of dependent time series.

# Simulated repeated sampling is implemented by the following code.

n.rep <- 10000
n <- 36
odd.ac.seq <- numeric(length=n.rep)
even.ac.seq <- numeric(length=n.rep)
for (i.rep in 1:n.rep) {
  iid.samp <- rnorm(n=n, mean=0, sd=1)
  did.samp <- sim.did.samp(iid.samp)
  odd.ac <- cor(did.samp[seq(from=1, to=n-1, by=2)], did.samp[seq(from=2, to=n, by=2)])
  even.ac <- cor(did.samp[seq(from=2, to=n-2, by=2)], did.samp[seq(from=3, to=n-1, by=2)])
  odd.ac.seq[i.rep] <- odd.ac
  even.ac.seq[i.rep] <- even.ac
}

# The code below generates histograms of the sampling distributions for the two autocorrelation statistics. 

par(mfrow = c(2, 1))
par(mar=c(2.00, 2.00, 1.00, 1.00)) #bottom, left, top, and right
hist(odd.ac.seq, xlim=c(-1,1), main="odd")
hist(even.ac.seq, xlim=c(-1,1), main="even")
par(mfrow = c(1, 1))

# We see that adjacent pairs with an odd first-entry tend to be strongly correlated in a positive direction.  Autocorrelations are much weaker within adjacent pairs with an odd first-entry. The following calculation produces an approximation to the probability that the latter autocorrelation is negative:

rel.freq <- length(which(even.ac.seq < 0)) / n.rep
rel.freq

# From this value we see that the odds are slightly better than even that this autocorrelation statistic would be negative.