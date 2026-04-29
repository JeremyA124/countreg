#' @importFrom stats model.frame model.response model.matrix pnorm qnorm dnbinom optim

glm_negb <- function(data,
                     formula,
                     offset = log(1)){
  #Parameter initliazations
  ##########################
  par <- model.frame(formula, data=data)
  y <- model.response(par)
  X <- model.matrix(formula, data=par)
  offset <- as.vector(model.offset(par))
  betas <- matrix(0, nrow=ncol(X),ncol=1)
  theta <- 1

  if(is.null(offset)){
    offset <- 0
  } else{
    offset <- log(offset)
  }

  #MLE estimation of theta
  ##########################
  theta.MLE <- function(par, curr.mu, y, neg=F){
    ll <- sum(log(dnbinom(y, size=par, mu=curr.mu)))
    return(ifelse(neg, -ll, ll))
  }

  #IWLS algorithm model fit
  ##########################
  repeat{
    eta <- offset + X %*% betas
    pred.means <- exp(eta)
    theta <- optim(par = theta,
                   fn=theta.MLE,
                   curr.mu=pred.means,
                   method = "Brent",
                   neg=T,
                   y=y, lower=0, upper=1000)$par
    tXW <- t(X * as.vector(pred.means^2/(pred.means^2/theta+pred.means)))
    tXWX <- tXW %*% X
    z <- (eta)+(y-pred.means)/pred.means
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
                    df.residuals = nrow(data)-length(betas),
                    call = match.call(),
                    terms = terms(formula),
                    model = list(y=y, x=X),
                    theta = theta))
  names(fit.dat$coefficients$betas) <- colnames(X)
  names(fit.dat$coefficients$std.error) <- colnames(X)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)

  if (theta>=999){
    warning("Theta diverges, perhaps use a Poisson regression model?")
  }

  class(fit.dat) <- "glm_negb"
  return(fit.dat)

}
