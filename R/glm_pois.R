#' @importFrom stats model.frame model.response model.matrix pnorm qnorm qt pt

glm_pois <- function(data,
                     formula,
                     quasi = F){

  #Parameter initialization
  ##########################
  par <- model.frame(formula, data=data)
  y <- model.response(par)
  X <- model.matrix(formula, data=par)
  offset <- as.vector(model.offset(par))
  betas <- matrix(0, nrow=ncol(X),ncol=1)

  if(is.null(offset)){
    offset <- 0
  } else{
    offset <- log(offset)
  }

  #IWLS algorithm model fit
  ##########################
  repeat{
    eta <- offset + X %*% betas
    pred.means <- exp(eta)
    tXW <- t(X * as.vector(pred.means))
    tXWX <- tXW %*% X
    z <- (eta)+(y-pred.means)/pred.means
    tXWz <- tXW %*% z
    betas.new <- solve(tXWX, tXWz)
    ss <- sum((betas.new-betas)**2)
    betas <- betas.new
    if(ss < 1e-6){
      disp.ratio <- sum(((y-pred.means)**2)/pred.means)/(nrow(data)-length(betas))
      if (quasi) {
        std.error <- sqrt(disp.ratio * diag(solve(tXWX)))
      } else {
        std.error <- sqrt(diag(solve(tXWX)))
      }
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
    if(quasi){
      p.vals[i] <- 2*(1-pt(abs(test.stat[i]), df=nrow(data)-length(betas)))
      crit <- qt(0.975, df=nrow(data)-length(betas))
      asymp.CI.lower[i] <- fit.dat$coefficients$betas[i]-crit*SE
      asymp.CI.higher[i] <- fit.dat$coefficients$betas[i]+crit*SE
    } else{
      p.vals[i] <- 2*(1-pnorm(abs(test.stat[i])))
      crit <- qnorm(0.975)
      asymp.CI.lower[i] <- fit.dat$coefficients$betas[i]-crit*SE
      asymp.CI.higher[i] <- fit.dat$coefficients$betas[i]+crit*SE
    }
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
                    model = list(y=y, x=X)))
  names(fit.dat$coefficients$betas) <- colnames(X)
  names(fit.dat$coefficients$std.error) <- colnames(X)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)

  # Quasi-Poisson add-ons
  ########################
  if(quasi){
    names(fit.dat$summary)[1] <- "t"
    fit.dat <- c(fit.dat, list(dispersion=disp.ratio))
    if(disp.ratio>0.9 & disp.ratio<1.1){
      warning("Dispersion ratio near 1, quasi unecessitated")
    }
  }

  class(fit.dat) <- "glm_pois"
  return(fit.dat)

}
