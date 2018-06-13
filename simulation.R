setwd("C:/Users/matte/Desktop/crt onesided noncompliance/att_estimation")

library(ranger); library(dplyr); library(distrEx); source("./helpers.R")

# Set parameters for simulation
set.seed(2436); nsim <- 500; counter <- 1; n <- 500; nsplits <- 2

# zvals = P(Z = 1); attvals = ATT; lambdaZvals <- E{lambdaZ(X)}
zvals <- c(0.2, 0.4, 0.5, 0.6, 0.8); attvals <- 3
lambda1vals <- 0.7; lambda0vals <- 0.3

# unobs = E{Y^0|A = 1, Z = 0, X} in [0,1], which cannot be observed
unobsvals <- 0.5

cols <- c("att", "mean.z", "mean.a", "lambda0", "lambda1", "unobs",
          "psi1","psi1.ci1", "psi1.ci2", 
          "psi2","psi2.ci1","psi2.ci2", 
          "psi1.var", "psi2.var", 
          "psi2.var.lo", "psi2.var.up",
          "eff.bound.lo", "eff.bound.up")

nrows1 <- nsim*length(attvals)*length(zvals)*length(lambda1vals)
nrows <- nrows1*length(lambda0vals)*length(unobsvals)

res <- as.data.frame(matrix(nrow=nrows, ncol=length(cols)))
colnames(res) <- cols

for(att in attvals) {
  
  for (zval in zvals) {
    
    for(unobs in unobsvals) {
      
      for(lambda1 in lambda1vals) {
        
        for(lambda0 in lambda0vals) {
          
          # calculate the value for alpha in P(Z = 1|X) = expit(alpha + x)
          # to achieve the desired value for P(Z = 1) = zval
          pz <- function(alpha) { 
            
            distrExIntegrate(function(x) { pi.x(alpha,x,FALSE)*dnorm(x) }, 
                             -200, 200) 
            
          }; alpha <- get_par(pz, -20, 20, zval)
          
          pa1 <- function(beta) { 
            
            distrExIntegrate(function(x) { lambda1.x(beta,x,TRUE)*dnorm(x) }, 
                             -200, 200) 
            
          }; beta <- get_par(pa1, -20, 20, lambda1)
          
          pa0 <- function(gamma) { 
            
            distrExIntegrate(function(x) { lambda0.x(gamma,x,TRUE)*dnorm(x) }, 
                             -200, 200) 
            
          }; gamma <- get_par(pa0, -200, 200, lambda0)
          
          # calculate p = P(A = 1)
          p <- distrExIntegrate(function(x) { 
            pa.x(alpha,beta,gamma,x,FALSE)*dnorm(x) }, 
                                -200, 200)
          
          num <- function(delta) { distrExIntegrate(function(x) { 
            
            (pi.x(alpha,x,TRUE)*mu.x(beta,gamma,delta,x,TRUE) + 
               lambda0.x(gamma,x,TRUE)*(mean.y(x,1,delta) - unobs))*dnorm(x) }, 
            
            -200, 200) }; att.fn <- function(delta) { num(delta)/p }
          
          delta <- get_par(att.fn, -200, 200, att)

          res[counter:(counter+nsim-1), c("eff.bound.lo")] <- 
            get_eff_bound(alpha, beta, gamma, delta, p, att, TRUE)[1]
          
          res[counter:(counter+nsim-1), c("eff.bound.up")] <- 
            get_eff_bound(alpha, beta, gamma, delta, p, att, TRUE)[2]
          
          for(i in 1:nsim) { 
            
            mess <- paste0("ATT = ", att, "; P(Z = 1) = ", zval,
                           "; P(A = 1) = ", round(p, 2), 
                           "; E{Y^0|A = 1, Z = 0, X} = ", unobs,
                           "; sim = ", i)
            if(i%%1 == 0) print(mess); flush.console() 
            
            # Simulate covariate (x), treat. assignment (z), treat. received (a)
            x <- rnorm(n); piz <- pi.x(alpha, x,TRUE); z <- rbinom(n, 1, piz)
            a <- rbinom(n, 1, z*lambda1.x(beta,x,TRUE) + 
                          (1-z)*lambda0.x(gamma,x,TRUE))
            y <- rnorm(n, mean.y(x,a,delta), 1)
            
            mu0hat <- mu1hat <- lambda1hat <- lambda0hat <- beta10hat <- 
              rep(NA,n)
            
            df <- data.frame(y = y, a = a, z = z, x = x)
            
            s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
            
            for (vfold in 1:nsplits) {
              
              train.row <- s != vfold; test.row <- s == vfold
              
              if (nsplits == 1) { train.row <- test.row }
              
              train <- df[train.row,]; test <- df[test.row,]
              train1 <- train[train$z==1,]; train0 <- train[train$z==0,]
              
              # Estimation of the nuisance regression functions
              
              if(sum(train0$a) > 0) {
                
                lambda0fit <- glm(a~x, data=train0, family=binomial)
                lambda0hat[test.row] <- predict(lambda0fit,newdata=test,
                                                type="response")
                
                beta10fit <- ranger(y*a~x, data=train0)
                beta10hat[test.row] <- predict(beta10fit, data=test)$predictions
                
              } else { lambda0hat[test.row] <- beta10hat[test.row] <- 0 }
              
              mu0fit <- ranger(y~x, data=train0)
              mu0hat[test.row] <- predict(mu0fit, data=test)$predictions
              
              mu1fit <- ranger(y~x, data=train1)
              mu1hat[test.row] <- predict(mu1fit, data=test)$predictions
              
              lambda1fit <- glm(a~x, data=train1, family=binomial)
              lambda1hat[test.row] <- predict(lambda1fit,newdata=test,
                                              type="response")
              
            }
            
            psi1 <- get_psi1(y,a,z,x,piz,mu0hat)
            
            psi2 <- get_psi2(y,a,z,x,piz,mu0hat,mu1hat,lambda1hat,lambda0hat,
                             beta10hat)
            
            res[counter, c("att", "mean.z", "mean.a", 
                           "lambda0", "lambda1", "unobs")] <- c(att, zval, p, 
                                                                lambda0, 
                                                                lambda1, unobs)
            
            res[counter, c("psi1","psi1.ci1", "psi1.ci2")] <- c(psi1$est, 
                                                                psi1$ci)
            res[counter, c("psi2","psi2.ci1", "psi2.ci2")] <- c(psi2$est, 
                                                                psi2$ci)
            
            res[counter, c("psi1.var", "psi2.var")] <- c(psi1$var, psi2$var)
            
            res[counter, c("psi2.var.lo", "psi2.var.up")] <- c(psi2$var.lo.bd,
                                                               psi2$var.up.bd)
            
            
            counter <- counter + 1
            
          }
        }
      }
    }
  }
}

agg.res <- res %>% group_by(att, mean.z, mean.a, lambda0, lambda1, unobs) %>%
  
  summarize(
    bias.psi1 = mean(psi1 - att),
    bias.psi2 = mean(psi2 - att),
    
    cvg.psi1 = mean(psi1.ci1 < att & att < psi1.ci2),
    cvg.psi2 = mean(psi2.ci1 < att & att < psi2.ci2),
    
    length.ci.psi1 = mean(psi1.ci2 - psi1.ci1),
    length.ci.psi2 = mean(psi2.ci2 - psi2.ci1),
    
    var.psi1 = mean(psi1.var),
    var.psi2 = mean(psi2.var),
    var.psi2.lo = mean(psi2.var.lo),
    var.psi2.up = mean(psi2.var.up),
    
    eff.ratio.1to2 = mean(psi1.var/psi2.var),
    
    eff.bound.lo = mean(eff.bound.lo),
    eff.bound.up = mean(eff.bound.up),
    
    psi2.to.eff.lo = mean(var.psi2.lo/eff.bound.lo),
    psi2.to.eff.up = mean(var.psi2.up/eff.bound.up))

rounded.agg.res <- round(agg.res, 2)

save(res, file = "./results_twosided.RData")
