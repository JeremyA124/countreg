glm_pois_zero <- function(data,
                          formula.pois,
                          formula.log,
                          offset = log(1)){

  #Parameter initliazations
  ##########################
  par1 <- stats::model.frame(formula.pois, data=data)
  par2 <- stats::model.frame(formula.log, data=data)
  y <- stats::model.response(par1)
  X.pois <- stats::model.matrix(formula.pois, data=par1)
  X.logit <- stats::model.matrix(formula.log, data=par2)
  betas <- matrix(0, nrow=ncol(X.pois), ncol=1)
  alphas <- matrix(0, nrow=ncol(X.logit), ncol=1)
  pred.means <- exp(offset)

  #IWLS algorithm model fit
  ##########################
  repeat{
    #Alpha Part
    ############
    eta.logit <- offset + X.logit %*% alphas
    pred.logs <- exp(eta.logit)
    pred.zeros <- pred.logs/(1+pred.logs)
    pred.zeros <- pmax(pmin(pred.zeros, 1 - 1e-100), 1e-100)
    delta <- (y == 0) * pred.zeros/(pred.zeros + (1-pred.zeros)*exp(-pred.means))
    tXW.logit <- t(X.logit * as.vector(pred.zeros*(1-pred.zeros)))
    tXWX.logit <- tXW.logit %*% X.logit
    z.logit <- eta.logit + (delta-pred.zeros)/(pred.zeros*(1-pred.zeros))
    tXWz.logit <- tXW.logit %*% z.logit
    new.alphas <- solve(tXWX.logit, tXWz.logit)
    ss1 <- sum((new.alphas-alphas)**2)
    alphas <- new.alphas
    #Poisson Part
    ##############
    eta <- offset + X.pois %*% betas
    pred.means <- exp(eta)
    tXW.pois <- t(X.pois * as.vector((1-delta)*pred.means))
    tXWX.pois <- tXW.pois %*% X.pois
    z <- eta+(y-pred.means)/pred.means
    tXWz.pois <- tXW.pois %*% z
    betas.new <- solve(tXWX.pois, tXWz.pois)
    ss2 <- sum((betas.new-betas)**2)
    betas <- betas.new
    if(ss1 < 1e-6 && ss2 < 1e-6){
      std.error.pois <- sqrt(diag(solve(tXWX.pois)))
      std.error.log <- sqrt(diag(solve(tXWX.logit)))
      break
    }
  }

  fit.dat <- list(coefficients=list(count=as.vector(betas),
                                    zero=as.vector(alphas),
                                    std.error.count=std.error.pois,
                                    std.error.zero=std.error.log))

  # Coefficient Pvals and Confidence Intervals
  #############################################
  test.stat <- rep(NA, times = length(fit.dat$coefficients$count))
  p.vals <- rep(NA, times = length(fit.dat$coefficients$count))
  asymp.CI.lower <- rep(NA, times = length(fit.dat$coefficients$count))
  asymp.CI.higher <- rep(NA, times = length(fit.dat$coefficients$count))
  for(i in 1:length(fit.dat$coefficients$count)){
    SE <- std.error.pois[i]
    test.stat[i] <- fit.dat$coefficients$count[i]/SE
    p.vals[i] <- 2*(1-stats::pnorm(abs(test.stat[i])))
    crit <- stats::qnorm(0.975)
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
    p.vals.alpha[i] <- 2*(1-stats::pnorm(abs(test.stat[i])))
    crit <- stats::qnorm(0.975)
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
                    terms = list(count=terms(formula.pois),
                                 zero=terms(formula.log)),
                    model = list(y=y, x=X.pois)))
  names(fit.dat$coefficients$count) <- colnames(X.pois)
  names(fit.dat$coefficients$zero) <- colnames(X.logit)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)
  names(fit.dat$pred.zero) <- 1:length(fit.dat$pred.zero)

  class(fit.dat) <- "glm_pois_zero"
  return(fit.dat)
}
