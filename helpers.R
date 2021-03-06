# helper functions

expit <- function(x){ exp(x)/(1+exp(x)) }; logit <- function(x) { log(x/(1-x))}

# to avoid problems with positivity assumption, we truncate the prop. scores
truncate_ps <- function(pscores, eps = 1e-8) {
  
  ps <- pscores
  ps[which(ps > 1-eps | is.nan(ps))] <- 1-eps; ps[which(ps < eps)] <- eps
  return(ps)
  
}

pi.x <- function(alpha, x, trunc) { 
  
  probs = expit(alpha + x + (x/10)^2 - 1*(x > 0.5))
  if(trunc) return(truncate_ps(probs))
  else return(expit(probs))
  
}

lambda1.x <- function(beta, x, trunc) { 
  
  probs = expit(beta + (x/10)^2 + 1*(x >= -0.5))
  
  if(trunc) return(truncate_ps(probs)) 
  else return(probs)
  
}

lambda0.x <- function(gamma, x, trunc) { 
  
  probs <- expit(gamma - (x/10)^3 + 1*(x >= 1.5))
  
  if(trunc) return(truncate_ps(probs))
  else return(expit(probs))
  
}

pa.x <- function(alpha, beta, gamma, x, trunc) {
  
  lambda1.x(beta,x,trunc)*pi.x(alpha,x,trunc) + 
    lambda0.x(gamma,x,trunc)*(1-pi.x(alpha,x,trunc))
  
}

b1 <- b2 <- b3 <- 1 # they are all set to 1 in the poster

# E{Y|X, A = a, Z = z} = E{Y|X, A = a}
mean.y <- function(x, a, delta, trunc) { 
  
  probs <- a*expit(-delta + x^2 + 1*(x <= 0.2)) + 
    (1-a)*(delta - x^2)
  
  if(trunc) return(truncate_ps(probs))
  if(!trunc) return(probs)
    
}

# calculate E{Y|X,Z=1} = E{Y|X,A=1,Z=1}lambda1(X) + E{Y|X,A=0,Z=1}(1-lambda1(X))
mu1.x <- function(beta, delta, x, trunc) {
  
  mean.y(x, 1, delta, trunc)*lambda1.x(beta, x, trunc) + 
    mean.y(x, 0, delta, trunc)*(1 - lambda1.x(beta, x, trunc))
  
}
# calculate E{Y|X,Z=0} = E{Y|X,A=1,Z=0}lambda0(X) + E{Y|X,A=0,Z=0}(1-lambda0(X))
mu0.x <- function(gamma, delta, x, trunc) {
  
  mean.y(x, 1, delta, trunc)*lambda0.x(gamma, x, trunc) + 
    mean.y(x, 0, delta, trunc)*(1 - lambda0.x(gamma, x, trunc))
  
}

mu.x <- function(beta, gamma, delta, x, trunc) {
  
  mu1.x(beta, delta, x, trunc) - mu0.x(gamma, delta, x, trunc)
  
}

lambda.x <- function(beta,gamma,x,trunc) {
  
  lambda1.x(beta,x,trunc) - lambda0.x(gamma,x,trunc)
  
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
get_phi <- function(t, z, pix, mu) {
  
  return(z/pix*(t - mu) + mu)
  
}

get_psi2 <- function(y, a, z, x, piz, mu0, mu1, lambda1, lambda0, beta10, 
                     experiment) {
  
  est.x.num.up <-  piz*(get_phi(y,z,piz,mu1) - get_phi(y,1-z,1-piz,mu0)) +
    get_phi(y*a,1-z,1-piz,beta10)
  
  est.x.den.up <- piz*(get_phi(a,z,piz,lambda1)) + 
    (1-piz)*get_phi(a,1-z,1-piz,lambda0)
  
  est.x.num.lo <- est.x.num.up - get_phi(a,1-z,1-piz,lambda0)
  est.x.den.lo <- est.x.den.up
  
  if (!experiment) {
    
    est.x.num.up <- est.x.num.up + (mu1 - mu0)*(z - piz)
    est.x.num.lo <- est.x.num.lo + (mu1 - mu0)*(z - piz)
    
    est.x.den.up <- est.x.den.up + (lambda1 - lambda0)*(z - piz)
    est.x.den.lo <- est.x.den.lo + (lambda1 - lambda0)*(z - piz)
    
  }
  
  est.up <- mean(est.x.num.up)/mean(est.x.den.up)
  est.lo <- mean(est.x.num.lo)/mean(est.x.den.lo)
  
  var.est.up <- var((est.x.num.up-est.up*est.x.den.up)/mean(est.x.den.up))
  se.up <- sqrt(var.est.up/length(y))
  var.est.lo <- var((est.x.num.lo-est.lo*est.x.den.lo)/mean(est.x.den.lo))
  se.lo <- sqrt(var.est.lo/length(y))
  
  calpha <- get_calpha(length(y), est.lo, est.up)
  
  ci <- c(est.lo - calpha*se.lo, est.up + calpha*se.up)

  out <- list(est = mean(c(est.lo, est.up)), 
              var.lo.bd = var.est.lo, var.up.bd = var.est.up, 
              ci = ci, calpha = calpha, est.lo = est.lo, est.up = est.up)
  
  return(out)
  
}

eff.infl.curve <- function(y, a, z, x, alpha, beta, gamma, delta, p, att, lo, 
                           trunc, experiment) {
  
  pz <- pi.x(alpha, x, trunc)
  
  pa0 <- lambda0.x(gamma, x, trunc); pa1 <- lambda1.x(beta, x, trunc)
  
  mu1 <- mu1.x(beta, delta, x, trunc)
  mu0 <- mu0.x(gamma, delta, x, trunc)
  
  beta10 <- mean.y(x, 1, delta, trunc)
  
  infl.curve.up <- 1/p*( pz*(get_phi(y,z,pz,mu1) - get_phi(y,1-z,1-pz,mu0) ) + 
                           get_phi(y*a,1-z,1-pz,pa0*beta10) - 
                           att*(pz*get_phi(a,z,pz,pa1) + 
                                  (1-pz)*get_phi(a,1-z,1-pz,pa0) ) )
  
  infl.curve.lo <- infl.curve.up - 1/p*( get_phi(a,1-z,1-pz,pa0) )
  
  if(!experiment) {
    
    infl.curve.up <- infl.curve.up + (mu1 - mu0 - att*(pa1 - pa0))*(z - pz)/p
    infl.curve.lo <- infl.curve.lo + (mu1 - mu0 - att*(pa1 - pa0))*(z - pz)/p
    
  }
  
  if(lo) out <- infl.curve.lo else out <- infl.curve.up
  
  return(out)
  
}

get_eff_bound <- function(alpha, beta, gamma, delta, p, att, trunc, experiment){
  
  inn.fn <- function(y, x, lo) {
    
    dy <- function(y, a, x) { 
      
      return( mean.y(x,a,delta,trunc)^y*( 1-mean.y(x,a,delta,trunc) )^(1-y) )

    } 
    
    da <- function(a, z, x) { 
      
      if(a==1 & z==1) return(lambda1.x(beta,x,trunc))
      if(a==1 & z==0) return(lambda0.x(gamma,x,trunc))
      if(a==0 & z==1) return(1-lambda1.x(beta,x,trunc))
      if(a==0 & z==0) return(1-lambda0.x(gamma,x,trunc))
      
    }
    
    dz <- function(z, x) { 
      
      return( pi.x(alpha,x,trunc)^z*(1-pi.x(alpha,x,trunc))^(1-z) )
      
    }
    
    dx <- function(x) { dnorm(x) }

    dens11 <- dy(y,1,x)*da(1,1,x)*dz(1,x)*dx(x)
    dens10 <- dy(y,1,x)*da(1,0,x)*dz(0,x)*dx(x)
    dens01 <- dy(y,0,x)*da(0,1,x)*dz(1,x)*dx(x)
    dens00 <- dy(y,0,x)*da(0,0,x)*dz(0,x)*dx(x)
    
    out <- {
      
      eff.infl.curve(y,1,1,x,alpha,beta,gamma,delta,p,att,lo,trunc,
                     experiment)^2*dens11 +
      eff.infl.curve(y,1,0,x,alpha,beta,gamma,delta,p,att,lo,trunc,
                     experiment)^2*dens10 +
      eff.infl.curve(y,0,1,x,alpha,beta,gamma,delta,p,att,lo,trunc,
                     experiment)^2*dens01 +
      eff.infl.curve(y,0,0,x,alpha,beta,gamma,delta,p,att,lo,trunc,
                     experiment)^2*dens00
      
    }

    return(out)

  }
  
  inn.fn1 <- Vectorize(function(x, lo) { inn.fn(0, x, lo) + inn.fn(1, x, lo) } )

  eff.bound.lo <- distrExIntegrate(function(x) { inn.fn1(x, TRUE) }, 
                                   -Inf, Inf)
  eff.bound.up <- distrExIntegrate(function(x) { inn.fn1(x, FALSE) }, 
                                   -Inf, Inf)
  
  out <- c(eff.bound.lo, eff.bound.up)
  
  if(any(out <= 0)) stop("Efficiency bounds are negative!")
  
  return(out)
  
}

get_difference <- function(alpha, beta, gamma, delta, p, att, trunc=TRUE) {
  
  inn.fn <- function(x) { 
    
    varz <- pi.x(alpha,x,TRUE)*(1-pi.x(alpha,x,TRUE))
    
    out <- ((mu.x(beta,gamma,delta,x,TRUE) - 
        att*lambda.x(beta,gamma,x,TRUE))^2*varz)/p^2
    
    return(out)
  }
  
  diff <- distrExIntegrate(function(x) { inn.fn(x)*dnorm(x) }, -Inf, Inf)
  
  return(diff)
  
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