#1)
n <- 19
x <- 7 #0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0
alpha <- 1/2 
beta <- 1/2
post.mean <- (alpha+x) / (alpha+beta+n)
post.mean          # 0.375

#2) 
qprobs <- c(0.28)            # quantiles to input
post.quant <- qbeta(p=qprobs, shape1=alpha+x, shape2=beta+n-x)
post.quant       # 0.3088125

#3) 
res2 <- rbeta(100000, alpha+x, beta+n-x)
ans <- median(sqrt(res2 * (1-res2)))
ans           # 0.4823846

#4) 
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
model.parms <- list()
model.parms$x <- x
model.parms$n <- n
model.parms$alpha <- alpha
model.parms$beta <- beta
model.parms$theta <- 0.37
comp.parms <- list()
comp.parms$n.samp <- 100000
comp.parms$theta.rad <- 0.75        # radius
result <- run.mcmc.binom(model.parms, comp.parms)
hist(result$samp, xlim=c(0, 1), xlab="theta")
result$rate             # 0.3387

#5) 
test.stat <- function(x.seq) {
  n <- length(x.seq)
  stat <- sum(abs(x.seq[2:n] - x.seq[1:(n-1)]))
  return(stat)
}

x.obs <- c(0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0)
obs.tstat <- test.stat(x.obs)
samp.theta <- rbeta(n=100000, shape1=alpha+x, shape2=beta+n-x)
x.rep <- numeric()
samp.tstat <- numeric()
for(i.samp in 1:100000){
  x.rep <- rbinom(n=model.parms$n, size=1, prob=samp.theta[i.samp])
  samp.tstat[i.samp] = test.stat(x.rep)
}

pval <- sum(obs.tstat>samp.tstat) / comp.parms$n.samp
pval

#6
library(expm)
test.stat2 <- function(res) {
  n <- length(res)
  stat <- mean(abs(res[1:n-1] + res[2:n]))
  return(stat)
}

x.vect <-c(-7.2, -1.9, -5.7, 0.2,-1.3, 8.7, 15.4, 22.1,24.0)
Z.mat <- as.matrix(cbind(rep(x=1, times=9), 1:9))
reg.summary <- function(x.vect, Z.mat) {
  n <- dim(Z.mat)[1]
  k <- dim(Z.mat)[2]
  nu <- n-k
  ZpZ <- t(Z.mat) %*% Z.mat
  ZpX <- t(Z.mat) %*% x.vect
  ZpZ.inv <- solve(ZpZ)
  sqrt.ZpZ.inv <- sqrtm(ZpZ.inv) # The sqrtm() function is matrix square-root 
  b.hat <- ZpZ.inv %*% ZpX
  res <- as.vector(x.vect - Z.mat %*% b.hat)
  s2 <- sum(res^2) / nu
  out.list <- list(sqrt.ZpZ.inv=sqrt.ZpZ.inv, b.hat=b.hat, s2=s2, nu=nu, res=res)
  return(out.list)
}
result <- reg.summary(x.vect, Z.mat)
b.hat <- result$b.hat
b.hat          # 4.125

# or 
n.samp <- 100000
n <- dim(Z.mat)[1]
k <- dim(Z.mat)[2]
result <- reg.summary(x.vect, Z.mat)
sqrt.ZpZ.inv <- result$sqrt.ZpZ.inv
b.hat <- result$b.hat
s2 <- result$s2
nu <- result$nu
logsigw2.samp <- numeric(length=n.samp)
beta.samp <- matrix(nrow=k, ncol=n.samp)
mu.samp <- matrix(nrow=n, ncol=n.samp)
for (i.samp in 1:n.samp) {
  sigw2 <- nu*s2 / rchisq(n=1, df=nu)
  sigw <- sqrt(sigw2)
  beta <- b.hat + sqrt.ZpZ.inv %*% as.matrix(rnorm(n=k, mean=0, sd=sigw))
  mu.vect <- Z.mat %*% beta
  res <- as.numeric(x.vect - mu.vect) #?
  pred.vect <- mu.vect + as.matrix(rnorm(n=n, mean=0, sd=sigw))
  obs.tstat <- test.stat2(res) #?
  samp.tstat[i.samp] <- test.stat2(rnorm(n=n, mean=0, sd=sigw)) #?
  logsigw2.samp[i.samp] <- log(sigw2)
  beta.samp[, i.samp] <- as.vector(beta)
  mu.samp[, i.samp] <- as.vector(Z.mat %*% beta)
}

mean(beta.samp[2,])            # 4.124639

#7) 
mean(exp(0.5*logsigw2.samp))        # 4.63201

#8
qprobs <- c(0.36)
new <- beta1 + beta2 *9
mut <- quantile(new, probs=qprobs)
mut          # 21.70804 

#9

#10 
pval <- sum(obs.tstat > samp.tstat) / comp.parms$n.samp
pval      

#11
x.11 <-  c(3.4, -0.6, -8.9, -5.7, -2.1, 5.5, 10.5, 4.6)
run.mcmc.m0.ar1 <- function(model.parms, comp.parms) {
  # ----- ----- subfunctions ----- -----
  log.dgen.dens <- function(x, sigw2, x0, phi) {
    n <- length(x)
    p <- length(phi)
    x.ext <- c(x0, x)
    Usum <- 0
    for (t in (p+1):(n+p)) {
      Usum <- Usum + (x.ext[t] - sum(phi*x.ext[seq(from=t-1, to=t-p)]))^2
    }
    val <- -0.5*n*log(2*pi*sigw2) - 0.5*Usum/sigw2
    return(val)
  }
  log.prior.dens <- function(sigw2, x0, phi) {
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
  propose.update.C <- function(curr.parm, radius) {
    cand.parm <- curr.parm + runif(n=1, min=-radius, max=radius)
    if (cand.parm < -1) {
      cand.parm <- -1 - (cand.parm + 1)
    } else if (cand.parm > 1) {
      cand.parm <- 1 - (cand.parm - 1)
    }
    return(cand.parm)
  }
  # ----- ----- --- main --- ----- -----
  x <- model.parms$x
  n <- length(model.parms$x)
  p <- length(model.parms$phi)
  m <- model.parms$m
  curr.sigw2 <- model.parms$sigw2
  curr.x0 <- 0
  curr.phi <- model.parms$phi
  n.samp <- comp.parms$n.samp
  sigw2.rad <- comp.parms$sigw2.rad
  phi.rad <- comp.parms$phi.rad
  cand.phi <- curr.phi
  accept.sigw2.cnt <- 0
  accept.phi.cnt <- rep(x=0, times=p)
  post.sigw2.samp <- numeric(length=n.samp)
  post.phi.samp <- matrix(data=NA, nrow=p, ncol=n.samp)
  post.predx.samp <- matrix(data=NA, nrow=m, ncol=n.samp)
  curr.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.x0, curr.phi)
  curr.log.prior.dens <- log.prior.dens(curr.sigw2, curr.x0, curr.phi)
  curr.log.post.dens <- curr.log.dgen.dens + curr.log.prior.dens
  for (i.samp in 1:n.samp) {
    # Update sigw2
    cand.sigw2 <- propose.update.B(curr.sigw2, sigw2.rad)
    cand.log.dgen.dens <- log.dgen.dens(x, cand.sigw2, curr.x0, curr.phi)
    cand.log.prior.dens <- log.prior.dens(cand.sigw2, curr.x0, curr.phi)
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
    # Update phi
    for (j in 1:p) {
      cand.phi[j] <- propose.update.C(curr.phi[j], phi.rad[j])
      cand.log.dgen.dens <- log.dgen.dens(x, curr.sigw2, curr.x0, cand.phi)
      cand.log.prior.dens <- log.prior.dens(curr.sigw2, curr.x0, cand.phi)
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
    pastx <- x[seq(from=n, to=n-p+1)]
    for (t in 1:m) {
      predx.mean <- sum(curr.phi*pastx)
      predx.sd <- sqrt(curr.sigw2)
      predx[t] <- rnorm(n=1, mean=predx.mean, sd=predx.sd)
      pastx <- c(predx[t], pastx[1:p-1])
    }
    # Record sample values
    post.sigw2.samp[i.samp] <- curr.sigw2
    post.phi.samp[,i.samp] <- curr.phi
    post.predx.samp[,i.samp] <- predx
  }
  samp <- list(sigw2=post.sigw2.samp, phi=post.phi.samp, predx=post.predx.samp)
  rate <- list(sigw2=accept.sigw2.cnt/n.samp, phi=accept.phi.cnt/n.samp)
  result <- list(samp=samp, rate=rate)
  return(result)
}
acf.seq <- acf(x.11, plot=FALSE)
acov.seq <- acf(x.11, type="covariance", plot=FALSE)
p <- 1
gamma0 <- acov.seq$acf[1]
gamma.vect <- acov.seq$acf[1+1:p]
Gamma.mat <- matrix(data=0, nrow=p, ncol=p)
Gamma.mat[, 1] <- acov.seq$acf[1:p]
for (j in 1:p) {
  Gamma.mat[1:(j-1),j] <- Gamma.mat[seq(from=j, to=1),1]
  Gamma.mat[j:p,j] <- Gamma.mat[1:(p-j+1),1]
}
Gamma.inv <- solve(Gamma.mat)

#11
phi.est <- Gamma.inv %*% gamma.vect
phi.est              # 0.5663477

sigw2.est <- as.numeric(gamma0 - t(phi.est) %*% gamma.vect)
sigw2.est

phi.sderr <- sqrt(sigw2.est*diag(Gamma.inv)[1:p] / length(x.11))
phi.sderr

#12
model.parms <- list()
model.parms$x <- x.11
model.parms$m <- 4
model.parms$sigw2 <- sigw2.est
model.parms$phi <- phi.est
comp.parms <- list()
comp.parms$n.samp <- 100000
comp.parms$sigw2.rad <- 1.5*sigw2.est
comp.parms$phi.rad <- phi.sderr
result <- run.mcmc.m0.ar1(model.parms, comp.parms)
qprobs <- c(0.28,0.5, 0.64)
quant.sigw2 <- quantile(result$samp$sigw2, probs=qprobs)
quant.phi1 <- quantile(quantile(result$samp$phi[1,], probs=qprobs))
quant.predx <- quantile(result$samp$predx[model.parms$m,], probs=qprobs)
      
# 12
quant.phi1 <- quantile(result$samp$phi[1,], probs=qprobs)
quant.phi1          # 0.6535969 

#13 
quant.sigw2 <- quantile(result$samp$sigw2, probs=qprobs)
quant.sigw2         # 29.15047 

#14
quant.predx <- quantile(result$samp$predx[model.parms$m,], probs=qprobs)
quant.predx     # -3.2227792  

