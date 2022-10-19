mu = 44
sigw = 8
phi1 = -0.25
phi2 = 0.21
phi = c(phi1, phi2)
n <- 100
q <- 25
p <- length(phi)
wn.samp <- rnorm (n=n, mean=0, sd = sigw)
alpha <- (1-sum(phi))*mu
max.m <- 2
pastx <- c(33.2, 50.4, 32.2, 39.8, 30.7)
predx <- numeric(length=max.m)

for (m in 1:max.m){
  predx[m] = phi1*pastx[1] + phi2*pastx[2] + alpha + wn.samp[m]
  pastx <- c(predx[m], pastx[1:p-1])
}
pred.err[1] <- sigw^2
for (m in 2:max.m){
  pred.err[m] <- pred.err[m-1] + sigw^2*psi[m-1]^2
}
predx



## number 14
theta = 0.4
sigw = 6
n = 6
x = c(-1.9, 3.5, 14.6, 8.7, 3.5, 2.7)

cov0 = (1+theta^2) * sigw^2
cov1 = theta * sigw^2

p12 = cov0 - ((cov1/cov0)^2 * cov0)
p23 = cov0 - ((cov1/p12)^2 * p12)
p34 = cov0 - ((cov1/p23)^2 * p23)
p45 = cov0 - ((cov1/p34)^2 * p34)
p56 = cov0 - ((cov1/p45)^2 * p45)
p67 = cov0 - ((cov1/p56)^2 * p56)

### 13
x01 = 0
x12 = (cov1/cov0)*(x[1] - x01)
x23 = (cov1/p12)*(x[2] - x12)
x34 = (cov1/p23)*(x[3] - x23)
x45 = (cov1/p34)*(x[4] - x34)
x56 = (cov1/p45)*(x[5] - x45)
x67 = (cov1/p56) * (x[6]-x56)


