glm_negb_zero <- function(data,
                          formula.negb,
                          formula.log){

  #MLE estimation of theta
  ##########################
  theta.MLE <- function(par, curr.mu, y, pis, weights, neg=F){
    ll <- sum(weights*log(dnbinom(y, mu=curr.mu, size=par)))
    return(ifelse(neg, -ll, ll))
  }

  #Parameter initliazations
  ##########################
  par1 <- model.frame(formula.negb, data=data)
  par2 <- model.frame(formula.log, data=data)
  y <- model.response(par1)
  X.negb <- model.matrix(formula.negb, data=par1)
  X.logit <- model.matrix(formula.log, data=par2)
  offset <- as.vector(model.offset(par1))
  betas <- matrix(0, nrow=ncol(X.negb), ncol=1)
  alphas <- matrix(0, nrow=ncol(X.logit), ncol=1)
  pred.means <- 1
  theta <- 1
  i <- 1
  maxrep=1000

  if(is.null(offset)){
    offset <- 0
  } else{
    offset <- log(offset)
  }

  repeat{
    #Alpha Part
    ############
    eta.logit <- offset + X.logit %*% alphas
    pred.logs <- exp(eta.logit)
    pred.zeros <- pred.logs/(1+pred.logs)
    delta <- (y == 0) * pred.zeros/(pred.zeros + (1-pred.zeros)*((theta/(pred.means+theta))**theta))
    tXW.logit <- t(X.logit * as.vector(pred.zeros*(1-pred.zeros)))
    tXWX.logit <- tXW.logit %*% X.logit
    z.logit <- (eta.logit) + (delta-pred.zeros)/(pred.zeros*(1-pred.zeros))
    tXWz.logit <- tXW.logit %*% z.logit
    new.alphas <- solve(tXWX.logit, tXWz.logit)
    ss1 <- sum((new.alphas-alphas)**2)
    alphas <- new.alphas
    #Negative Binomial Part
    ########################
    eta <- offset + X.negb %*% betas
    pred.means <- exp(eta)
    theta <- optim(par=theta,
                          fn=theta.MLE,
                          curr.mu=pred.means,
                          weights=1-delta,
                          method="Brent",
                          neg=T,
                          y=y,
                          lower=1e-4,
                          upper=1000)$par
    tXW.negb <- t(X.negb * as.vector((1-delta)*(pred.means**2/((pred.means**2/theta+pred.means)))))
    tXWX.negb <- tXW.negb %*% X.negb
    z <- (eta)+(y-pred.means)*(1/pred.means)
    tXWz.negb <- tXW.negb %*% z
    betas.new <- solve(tXWX.negb, tXWz.negb)
    ss2 <- sum((betas.new-betas)**2)
    betas <- betas.new
    if(ss1 < 1e-6 && ss2 < 1e-6){
      std.error.negb <- sqrt(diag(solve(tXWX.negb)))
      std.error.log <- sqrt(diag(solve(tXWX.logit)))
      break
    }else{
      if(i == maxrep){
        std.error.negb <- sqrt(diag(solve(tXWX.negb)))
        std.error.log <- sqrt(diag(solve(tXWX.logit)))
        warning("Maximum number of iterations reached.")
        break
      }else{
        i <- i+1
      }
    }
  }

  fit.dat <- list(coefficients=list(count=as.vector(betas),
                                    zero=as.vector(alphas),
                                    std.error.count=std.error.negb,
                                    std.error.zero=std.error.log))

  # Coefficient Pvals and Confidence Intervals
  #############################################
  test.stat <- rep(NA, times = length(fit.dat$coefficients$count))
  p.vals <- rep(NA, times = length(fit.dat$coefficients$count))
  asymp.CI.lower <- rep(NA, times = length(fit.dat$coefficients$count))
  asymp.CI.higher <- rep(NA, times = length(fit.dat$coefficients$count))
  for(i in 1:length(fit.dat$coefficients$count)){
    SE <- std.error.negb[i]
    test.stat[i] <- fit.dat$coefficients$count[i]/SE
    p.vals[i] <- 2*(1-pnorm(abs(test.stat[i])))
    crit <- qnorm(0.975)
    asymp.CI.lower[i] <- fit.dat$coefficients$count[i]-crit*SE
    asymp.CI.higher[i] <- fit.dat$coefficients$count[i]+crit*SE
  }

  test.stat.alpha <- rep(NA, times = length(fit.dat$coefficients$zero))
  p.vals.alpha <- rep(NA, times = length(fit.dat$coefficients$zero))
  asymp.CI.lower.alpha <- rep(NA, times = length(fit.dat$coefficients$zero))
  asymp.CI.higher.alpha <- rep(NA, times = length(fit.dat$coefficients$zero))
  for(i in 1:length(fit.dat$coefficients$zero)){
    SE <- std.error.log[i]
    test.stat.alpha[i] <- fit.dat$coefficients$zero[i]/SE
    p.vals.alpha[i] <- 2*(1-pnorm(abs(test.stat[i])))
    crit <- qnorm(0.975)
    asymp.CI.lower.alpha[i] <- fit.dat$coefficients$zero[i]-crit*SE
    asymp.CI.higher.alpha[i] <- fit.dat$coefficients$zero[i]+crit*SE
  }

  # Model summary/data
  #####################
  fit.dat <- c(fit.dat,
               list(residuals = as.list(y - pred.means),
                    summary = list(count=list(Z=test.stat,
                                              asymp.CI.lower=asymp.CI.lower,
                                              asymp.CI.higher=asymp.CI.higher,
                                              p.vals=p.vals),
                                   zero=list(Z=test.stat.alpha,
                                             asymp.CI.lower=asymp.CI.lower.alpha,
                                             asymp.CI.higher=asymp.CI.higher.alpha,
                                             p.vals=p.vals.alpha)),
                    fitted.values = as.list((1-pred.zeros)*pred.means),
                    df.residuals = nrow(data)-length(betas)-length(alphas),
                    pred.zero = as.list(pred.zeros),
                    call = match.call(),
                    terms = list(count=terms(formula.negb),
                                 zero=terms(formula.log)),
                    model = list(y=y, x=X.negb),
                    theta = theta))
  names(fit.dat$coefficients$count) <- colnames(X.negb)
  names(fit.dat$coefficients$zero) <- colnames(X.logit)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)
  names(fit.dat$pred.zero) <- 1:length(fit.dat$pred.zero)

  if (theta>=999){
    warning("Theta diverges, perhaps use a Poisson regression model?")
  }

  class(fit.dat) <- "glm_negb_zero"
  return(fit.dat)
}
