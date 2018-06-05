# helper functions

expit <- function(x){ exp(x)/(1+exp(x)) }; logit <- function(x) { log(x/(1-x))}

b0 <- 1 #b0 = 0 means CRT

pi.x <- function(alpha, x) { 
  
  pz <- expit(alpha + b0*x) 
  
  pz[which(pz > 0.999 | is.nan(pz))] <- 0.999
  pz[which(pz < 0.001)] <- 0.001
  
  return(pz)
  
}

lambda1.x <- function(beta, x) { 
  
  pa1 <- expit(beta + 2*x) 
  pa1[which(pa1 > 0.999 | is.nan(pa1))] <- 0.999
  pa1[which(pa1 < 0.001)] <- 0.001
  
  return(pa1)
  
}

lambda0.x <- function(gamma, x) { 
  
  pa0 <- expit(gamma - x) 
  
  pa0[which(pa0 > 0.999 | is.nan(pa0))] <- 0.999
  pa0[which(pa0 < 0.001)] <- 0.001
  
  return(pa0)
  
}

pa.x <- function(alpha, beta, gamma, x) {
  
  lambda1.x(beta,x)*pi.x(alpha,x) + lambda0.x(gamma,x)*(1-pi.x(alpha,x))
  
}

b1 <- b2 <- b3 <- 1 # they are all set to 1 in the poster

# E{Y|X, A=a, Z=z} = E{Y|X, A = a}
mean.y <- function(x, a, delta) { b1 + b2*x + b3*a + delta*a*x }

# calculate E{Y|X,Z=1} = E{Y|X,A=1,Z=1}lambda1(X) + E{Y|X,A=0,Z=1}(1-lambda1(X))
mu1.x <- function(beta, delta, x) {
  
  mean.y(x, 1, delta)*lambda1.x(beta, x) + 
    mean.y(x, 0, delta)*(1 - lambda1.x(beta, x))
  
}

mu0.x <- function(gamma, delta, x) {
  
  mean.y(x, 1, delta)*lambda0.x(gamma, x) + 
    mean.y(x, 0, delta)*(1 - lambda0.x(gamma, x))
  
}

mu.x <- function(beta, gamma, delta, x) {
  
  mu1.x(beta, delta, x) - mu0.x(gamma, delta, x)
  
}

# compute other DR estimator (from plug-in estimator in Frolich & Melly (2013))
get_psi1 <- function(y, a, z, x, piz, mu0) {
  
  est.x <- (y - (1-z)/(1-piz)*(y - mu0) - mu0)/mean(a)
  var.est <- var(est.x); se <- sqrt(var.est/length(y)); est <- mean(est.x)
  
  out <- list(est = est, se = se, var = var.est, ci = est + c(-1, 1)*1.96*se)
  
  return(out)
  
}

get_calpha <- function(n, lo, hi) {
  
  delta.seq <- seq(-5, 5, length.out = 1e4)
  fn <- function(calpha) { pnorm(calpha + sqrt(n)*(hi - lo)/max(c(hi,lo))) -
      pnorm(-calpha) }
  vals <- sapply(delta.seq, fn)
  delta <- delta.seq[which.min(abs(0.95 - vals))]
  
}

# compute proposed estimator
# mu0 <- mu0hat; mu1 <- mu1hat; lambda1 <- lambda1hat; lambda0 <- lambda0hat
# beta10 <- beta10hat
get_psi2 <- function(y, a, z, x, piz, mu0, mu1, lambda1, lambda0, beta10) {
  
  est.x.num.up <-  piz*(z/piz*(y - mu1) - (1-z)/(1-piz)*(y - mu0) + mu1 - mu0) +
    (1-z)/(1-piz)*(y*a - lambda0*beta10) + lambda0*beta10
  est.x.den.up <- a - (z-piz)*(lambda1 - lambda0)
  
  est.x.num.lo <- est.x.num.up - (a - lambda0)*(1-z)/(1-piz) - lambda0
  est.x.den.lo <- est.x.den.up
  
  est.up <- mean(est.x.num.up)/mean(est.x.den.up)
  est.lo <- mean(est.x.num.lo)/mean(est.x.den.lo)
  
  var.est.up <- var((est.x.num.up-est.up*est.x.den.up)/mean(est.x.den.up))
  se.up <- sqrt(var.est.up/length(y))
  var.est.lo <- var((est.x.num.lo-est.lo*est.x.den.lo)/mean(est.x.den.lo))
  se.lo <- sqrt(var.est.lo/length(y))
  
  calpha <- get_calpha(length(y), est.lo, est.up)
  
  ci <- c(est.lo - calpha*se.lo, est.up + calpha*se.up)
  
  se <- (ci[2] - ci[1])/(2*calpha)*sqrt(n)
  
  out <- list(est = mean(c(est.lo, est.up)), se = se, 
              var = se^2, var.lo.bd  = var.est.lo, var.up.bd  = var.est.up, 
              ci = ci)
  
  return(out)
  
}

# delta(x) = E{Y|A = 1, Z = 1, X} - E{Y|A = 0, Z = 1, X}
beta.x <- function(delta, x) {
  
  mean.y(x, 1, delta) - mean.y(x, 0, delta)
  
}

eff.infl.curve <- function(y, a, z, x, alpha, beta, gamma, delta, p, att, lo) {
  
  pz <- pi.x(alpha, x); pa0 <- lambda0.x(gamma, x); pa1 <- lambda1.x(beta, x)
  
  mu1 <- mu1.x(beta, delta, x); mu0 <- mu0.x(gamma, delta, x)
  beta10 <- mean.y(x, 1, delta)
  mu <- mu1 - mu0; lambda <- pa1 - pa0
  
  infl.curve.up <- z*(y - mu1) - (1-z)*pz/(1-pz)*(y - mu0) + pz*mu +
    
    (1-z)/(1-pz)*( y*a - pa0*beta10 ) + pa0*beta10 -
    
    att*( pz*(z/pz*(a - pa1) + pa1) + (1-pz)*((1-z)/(1-pz)*(a - pa0) + pa0) )
  
  infl.curve.lo <- infl.curve.up - (a - pa0)*(1-z)/(1-pz) - pa0
  
  out <- ifelse(lo, infl.curve.lo/p, infl.curve.up/p)
  
  return(out)
  
}

get_eff_bound <- function(alpha, beta, gamma, delta, p, att) {
  
  inn.fn1 <- Vectorize(function(a, z, x, lo) {
    
    distrExIntegrate(function(y) { 
      
      eff.infl.curve(y, a, z, x, alpha, beta, gamma, delta, p, att, 
                     lo)^2*dnorm(y, mean.y(x,a,delta)) }, -200, 200) } 
    )
  
  inn.fn2 <- function(x,lo) {
    
    inn.fn1(1, 1, x,lo)*lambda1.x(beta, x)*pi.x(alpha, x) +
      inn.fn1(0, 1, x,lo)*(1 - lambda1.x(beta, x))*pi.x(alpha, x) +
      inn.fn1(1, 0, x,lo)*lambda0.x(gamma, x)*(1 - pi.x(alpha, x)) +
      inn.fn1(0, 0, x,lo)*(1 - lambda0.x(gamma, x))*(1 - pi.x(alpha, x))
    
  }
  
  eff.bound.up <- distrExIntegrate(function(x) { inn.fn2(x, FALSE)*dnorm(x) }, 
                                   -200, 200)
  eff.bound.lo <- distrExIntegrate(function(x) { inn.fn2(x, TRUE)*dnorm(x) }, 
                                   -200, 200)
  
  out <- c(eff.bound.lo, eff.bound.up)
  
  if(any(out) <= 0) stop("Some variances are negative!")
  
  return(out)
  
}

get_par <- function(fn, lo, hi, par) {
  
  delta.seq <- seq(lo, hi, length.out = 1e4)
  vals <- sapply(delta.seq, fn)
  delta <- delta.seq[which.min(abs(par - vals))]
  
  if(abs(fn(delta) - par) >= 0.1)  {
    stop(paste0("Not able to find value to achieve desired parameter.",
                "\n Desired: ", par,
                "\n Achieved: ", 
                round(fn(delta), 2)))
  }
  
  return(delta)
  
}
