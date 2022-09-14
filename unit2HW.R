dat <- read.table('financial C.txt', header=FALSE)

# 7)
# SACF at lag h=1 for data
x <- acf(dat)
lag1 <- x$acf[2]
print(lag1)

#8) 
sd <- 1/sqrt(nrow(dat))
print(sd)

#9) 
corr <- c(1,0.42, 0.43, 0.35, 0.41, 0.32, 0.29, 0.24, 0.20, 0.21)

p = pacf(corr, 9, 0)
p


#10) 
# HW QUESTION 9,10
x0 <- c(64.7, 66.6, 61.9, 66.6, 68.3)
x0 <- c(84.7, 84.5, 85.9, 86.2, 87.8)
n <- length(x0)
m <- 2

# Because the basic formulas for linear prediction are developed for a mean-zero time series, we subtract the sample mean of time series from these values, then make our prediction, and then add the sample mean back at the end.

xbar <- 65
xbar <- 85
x0.shift <- x0 - xbar 

# We develop a linear prediction formula using the property that its coefficients are determined by the underlying ACF of the time series. We do not actually know the underlying ACF, but we can instead base our work on the sample ACF of the existing data set.
rec <- c(1.00, -0.50, 0.42, -0.45, 0.38, -0.42, 0.34, -0.22, 0.28, -0.27)
rec <- c(1.00, 0.42, 0.43, 0.35, 0.41, 0.32, 0.29, 0.24, 0.20, 0.21) 
acf.seq <- acf(rec, plot=FALSE)$acf

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

result <- blpred(m, x0, rec, xbar)
x.pred <- result$pred
result <- durbin.levinson.calc(n=n, acf.seq)
x.pacf <- result$pacf
print(x.pacf)

#9
result <- durbin.levinson.calc(n=9, rec)
x.pacf <- result$pacf
