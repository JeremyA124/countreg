glm_pois_GP2 <- function(data,
                         formula){

  dgpois <- function(y, lambda, alpha){
    probs <- numeric(length(y))
    comp1 <- lambda/(1+alpha*lambda)
    comp2 <- (lambda+alpha*y)^(y-1)/factorial(y)
    comp3 <- exp((-lambda*(1+alpha*y))/(1+alpha*lambda))
    probs <- comp1*comp2*comp3
    return(probs)
  }

  alpha.MLE <- function(par,
                        y,
                        fits,
                        neg=F){
    ll <- sum(log(dgpois(y, fits, par)))
    return(ifelse(neg, -ll, ll))
  }

  #Parameter initliazations
  ##########################
  par <- model.frame(formula, data=data)
  y <- model.response(par)
  X <- model.matrix(formula, data=par)
  offset <- model.offset(par)
  betas <- matrix(0, nrow=ncol(X),ncol=1)
  alpha <- 1

  if(is.null(offset)){
    offset <- 0
  }

  #IWLS algorithm model fit
  ##########################
  repeat{
    eta <- offset + X %*% betas
    pred.means <- exp(eta)
    alpha <- optim(par = alpha,
                   fn=alpha.MLE,
                   fits=pred.means,
                   method = "Brent",
                   neg=T,
                   y=y, lower=1e-4, upper=1000)$par
    tXW <- t(X * as.vector(pred.means/(1+alpha*pred.means)^2))
    tXWX <- tXW %*% X
    z <- eta+(y-pred.means)/(pred.means)
    tXWz <- tXW %*% z
    betas.new <- solve(tXWX, tXWz)
    ss <- sum((betas.new-betas)**2)
    betas <- betas.new
    if(ss < 1e-6){
        std.error <- sqrt(diag(solve(tXWX)))
        break
    }
  }

  fit.dat <- list(coefficients=list(betas=as.vector(betas),
                                    std.error=std.error))

  # Coefficient Pvals and Confidence Intervals
  #############################################
  test.stat <- rep(NA, times = length(fit.dat$coefficients$betas))
  p.vals <- rep(NA, times = length(fit.dat$coefficients$betas))
  asymp.CI.lower <- rep(NA, times = length(fit.dat$coefficients$betas))
  asymp.CI.higher <- rep(NA, times = length(fit.dat$coefficients$betas))
  for(i in 1:length(fit.dat$coefficients$betas)){
    SE <- std.error[i]
    test.stat[i] <- fit.dat$coefficients$betas[i]/SE
    p.vals[i] <- 2*(1-pnorm(abs(test.stat[i])))
    crit <- qnorm(0.975)
    asymp.CI.lower[i] <- fit.dat$coefficients$betas[i]-crit*SE
    asymp.CI.higher[i] <- fit.dat$coefficients$betas[i]+crit*SE
  }

  # Model summary/data
  #####################
  fit.dat <- c(fit.dat,
               list(residuals = as.list(y - pred.means),
                    summary = list(Z=test.stat,
                                   asymp.CI.lower=asymp.CI.lower,
                                   asymp.CI.higher=asymp.CI.higher,
                                   p.vals=p.vals),
                    fitted.values = as.list(pred.means),
                    df.residuals = nrow(data)-length(betas)-1,
                    call = match.call(),
                    terms = terms(formula),
                    model = list(y=y, x=X),
                    alpha=alpha))
  names(fit.dat$coefficients$betas) <- colnames(X)
  names(fit.dat$coefficients$std.error) <- colnames(X)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)

  class(fit.dat) <- "glm_pois_GP2"
  return(fit.dat)
}
