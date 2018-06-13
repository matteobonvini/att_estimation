# helper functions

expit <- function(x){ exp(x)/(1+exp(x)) }; logit <- function(x) { log(x/(1-x))}

b0 <- 1 #b0 = 0 means CRT


# to avoid problems with positivity assumption, we truncate the prop. scores
truncate_ps <- function(pscores) {
  
  eps <- 1e-8
  ps <- pscores
  ps[which(ps > 1-eps | is.nan(ps))] <- 1-eps; ps[which(ps < eps)] <- eps
  return(ps)
  
}

pi.x <- function(alpha, x, trunc) { 
  
  if(trunc) return(truncate_ps(expit(alpha + b0*x)))
  else return(expit(alpha + b0*x))
  
}

lambda1.x <- function(beta, x, trunc) { 
  
  if(trunc) return(truncate_ps(expit(beta - x))) 
  else return(expit(beta - x))
  
}

lambda0.x <- function(gamma, x, trunc) { 
  
  if(trunc) return(truncate_ps(expit(gamma + x/2)))
  else return(expit(gamma + x/2))
  
}

pa.x <- function(alpha, beta, gamma, x, trunc) {
  
  lambda1.x(beta,x,trunc)*pi.x(alpha,x,trunc) + 
    lambda0.x(gamma,x,trunc)*(1-pi.x(alpha,x,trunc))
  
}

b1 <- b2 <- b3 <- 1 # they are all set to 1 in the poster

# E{Y|X, A=a, Z=z} = E{Y|X, A = a}
mean.y <- function(x, a, delta) { b1 + b2*x + b3*a + delta*a*x }

# calculate E{Y|X,Z=1} = E{Y|X,A=1,Z=1}lambda1(X) + E{Y|X,A=0,Z=1}(1-lambda1(X))
mu1.x <- function(beta, delta, x, trunc) {
  
  mean.y(x, 1, delta)*lambda1.x(beta, x, trunc) + 
    mean.y(x, 0, delta)*(1 - lambda1.x(beta, x, trunc))
  
}
# calculate E{Y|X,Z=0} = E{Y|X,A=1,Z=0}lambda0(X) + E{Y|X,A=0,Z=0}(1-lambda0(X))
mu0.x <- function(gamma, delta, x, trunc) {
  
  mean.y(x, 1, delta)*lambda0.x(gamma, x, trunc) + 
    mean.y(x, 0, delta)*(1 - lambda0.x(gamma, x, trunc))
  
}

mu.x <- function(beta, gamma, delta, x, trunc) {
  
  mu1.x(beta, delta, x, trunc) - mu0.x(gamma, delta, x, trunc)
  
}

# compute other DR estimator (from plug-in estimator in Frolich & Melly (2013))
get_psi1 <- function(y, a, z, x, piz, mu0) {
  
  est.x <- (y - (1-z)/(1-piz)*(y - mu0) - mu0)/mean(a)
  var.est <- var(est.x); se <- sqrt(var.est/length(y)); est <- mean(est.x)
  
  out <- list(est = est, se = se, var = var.est, ci = est + c(-1, 1)*1.96*se)
  
  return(out)
  
}

# Using formula from Imbens & Mansky (2004)
get_calpha <- function(n, lo, hi) {
  
  delta.seq <- seq(-5, 5, length.out = 1e4)
  
  fn <- function(calpha) { pnorm(calpha + sqrt(n)*(hi - lo)/max(c(hi,lo))) -
      pnorm(-calpha) }; vals <- sapply(delta.seq, fn)
  
  delta <- delta.seq[which.min(abs(0.95 - vals))]
  
}

# compute proposed estimator
get_psi2 <- function(y, a, z, x, piz, mu0, mu1, lambda1, lambda0, beta10) {
  
  est.x.num.up <-  piz*(z/piz*(y - mu1) - (1-z)/(1-piz)*(y - mu0) + mu1 - mu0) +
    (1-z)/(1-piz)*(y*a - beta10) + beta10
  
  est.x.den.up <- piz*(z/piz*(a - lambda1) + lambda1) +
    (1-piz)*((1-z)/(1-piz)*(a - lambda0) + lambda0)
  
  est.x.num.lo <- est.x.num.up - (1-z)/(1-piz)*(a - lambda0) - lambda0
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
              var = se^2, var.lo.bd = var.est.lo, var.up.bd = var.est.up, 
              ci = ci, calpha = calpha)
  
  return(out)
  
}

eff.infl.curve <- function(y, a, z, x, alpha, beta, gamma, delta, p, att, lo, 
                           trunc) {
  
  pz <- pi.x(alpha, x, trunc)
  
  pa0 <- lambda0.x(gamma, x, trunc); pa1 <- lambda1.x(beta, x, trunc)
  
  mu1 <- mu1.x(beta, delta, x, trunc)
  mu0 <- mu0.x(gamma, delta, x, trunc)
  
  beta10 <- mean.y(x, 1, delta)
  
  infl.curve.up <- 1/p*( pz*(z/pz*(y - mu1) - (1-z)/(1-pz)*(y - mu0) + mu1 - mu0) +
    
    (1-z)/(1-pz)*( y*a - pa0*beta10 ) + pa0*beta10 -
    
    att*( pz*(z/pz*(a - pa1) + pa1) + (1-pz)*((1-z)/(1-pz)*(a - pa0) + pa0) ) )
  
  infl.curve.lo <- infl.curve.up - 1/p*( (a - pa0)*(1-z)/(1-pz) - pa0 )
  
  if(lo) out <- infl.curve.lo else out <- infl.curve.up
  
  return(out)
  
}

get_eff_bound <- function(alpha, beta, gamma, delta, p, att, trunc) {
  
  inn.fn <- function(y, x, lo) {
    
    dy <- function(y, a, x) { return(dnorm(y, mean.y(x, a, delta))) } 
    da <- function(a, z, x) { 
      
      if(a==1 & z==1) return(lambda1.x(beta,x,trunc))
      if(a==1 & z==0) return(lambda0.x(gamma,x,trunc))
      if(a==0 & z==1) return(1-lambda1.x(beta,x,trunc))
      if(a==0 & z==0) return(1-lambda0.x(gamma,x,trunc))
      
    }
    dz <- function(z, x) {
      
      if(z==1) return(pi.x(alpha,x,trunc))
      if(z==0) return(1-pi.x(alpha,x,trunc))
      
    }
    dx <- function(x) { dnorm(x) }

    dens11 <- dy(y,1,x)*da(1,1,x)*dz(1,x)*dx(x)
    dens10 <- dy(y,1,x)*da(1,0,x)*dz(0,x)*dx(x)
    dens01 <- dy(y,0,x)*da(0,1,x)*dz(1,x)*dx(x)
    dens00 <- dy(y,0,x)*da(0,0,x)*dz(0,x)*dx(x)
    
    out <- {
      
      eff.infl.curve(y,1,1,x,alpha,beta,gamma,delta,p,att,lo,trunc)^2*dens11 +
      eff.infl.curve(y,1,0,x,alpha,beta,gamma,delta,p,att,lo,trunc)^2*dens10 +
      eff.infl.curve(y,0,1,x,alpha,beta,gamma,delta,p,att,lo,trunc)^2*dens01 +
      eff.infl.curve(y,0,0,x,alpha,beta,gamma,delta,p,att,lo,trunc)^2*dens00
      
    }

    return(out)

  }

  inn.fn1 <- Vectorize(function(x, lo) {
    distrExIntegrate(function(y) { inn.fn(y, x, lo) }, -500, 500)
  })

  eff.bound.lo <- distrExIntegrate(function(x) { inn.fn1(x, TRUE) }, -500, 500)
  eff.bound.up <- distrExIntegrate(function(x) { inn.fn1(x, FALSE) }, -500, 500)
  
  out <- c(eff.bound.lo, eff.bound.up)
  
  if(any(out <= 0)) stop("Efficiency bounds are negative!")
  
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
beta.x <- function(alpha, gamma, beta, x) {
  
  mean.y(x, 1, beta) - mean.y(x, 0, beta)
  
}

