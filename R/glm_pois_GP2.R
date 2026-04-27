glm_pois_GP2 <- function(data,
                         formula,
                         offset = log(1)){

  alpha.estimate <- function(y,
                             fits,
                             num_obs,
                             num_parm){
    df <- num_obs-num_parm
    num <- sum((abs(y-fits)/sqrt(fits)-1)*fits^(-1))

    return(num/df)
  }

  #Parameter initliazations
  ##########################
  par <- stats::model.frame(formula, data=data)
  y <- stats::model.response(par)
  X <- stats::model.matrix(formula, data=par)
  betas <- matrix(0, nrow=ncol(X),ncol=1)

  #IWLS algorithm model fit
  ##########################
  repeat{
    eta <- offset + X %*% betas
    pred.means <- exp(eta)
    alpha <- alpha.estimate(y,
                            pred.means,
                            nrow(data),
                            ncol(X))
    tXW <- t(X * as.vector(1/(1+pred.means*alpha)^2))
    tXWX <- tXW %*% X
    z <- eta+(y-pred.means)/(pred.means*(1+pred.means*alpha)^2)
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
  test.stat <- rep(NA, times = length(fit.dat$coefficients))
  p.vals <- rep(NA, times = length(fit.dat$coefficients))
  asymp.CI.lower <- rep(NA, times = length(fit.dat$coefficients))
  asymp.CI.higher <- rep(NA, times = length(fit.dat$coefficients))
  for(i in 1:length(fit.dat$coefficients)){
    SE <- std.error[i]
    test.stat[i] <- fit.dat$coefficients$betas[i]/SE
    p.vals[i] <- 2*(1-stats::pnorm(abs(test.stat[i])))
    crit <- stats::qnorm(0.975)
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
                    alpha=alpha))
  names(fit.dat$coefficients$betas) <- colnames(X)
  names(fit.dat$coefficients$std.error) <- colnames(X)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)

  class(fit.dat) <- "glm_pois_GP2"
  return(fit.dat)
}
