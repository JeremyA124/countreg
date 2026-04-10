#' @importFrom stats model.frame model.response model.matrix pnorm qnorm qt pt

glm_pois <- function(data,
                     formula,
                     offset = log(1),
                     quasi = F){

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
    W <- diag(as.vector(pred.means))
    z <- eta+(y-pred.means)/pred.means
    inv.Hess <- solve(t(X)%*%W%*%X)
    betas.new <- inv.Hess%*%t(X)%*%W%*%z
    ss <- sum((betas.new-betas)**2)
    betas <- betas.new
    std.error <- sqrt(diag(inv.Hess))
    if(ss < 1e-6){
      break
    }
  }

  if(quasi){
    disp.ratio <- sum(((y-pred.means)**2)/pred.means)/(nrow(data)-length(betas))
    inv.Hess <- disp.ratio*inv.Hess
  }

  fit.dat <- data.frame(betas=betas,std.error=std.error)

  # Coefficient Pvals and Confidence Intervals
  #############################################
  test.stat <- rep(NA, times = nrow(fit.dat))
  p.vals <- rep(NA, times = nrow(fit.dat))
  asymp.CI.lower <- rep(NA, times = nrow(fit.dat))
  asymp.CI.higher <- rep(NA, times = nrow(fit.dat))
  for(i in 1:nrow(fit.dat)){
    SE <- sqrt(inv.Hess[i,i])
    test.stat[i] <- fit.dat$betas[i]/SE
    if(quasi){
      p.vals[i] <- 2*(1-stats::pt(abs(test.stat[i]), df=nrow(data)-length(betas)))
      crit <- stats::qt(0.975, df=nrow(data)-length(betas))
      asymp.CI.lower[i] <- fit.dat$betas[i]-crit*SE
      asymp.CI.higher[i] <- fit.dat$betas[i]+crit*SE
    } else{
      p.vals[i] <- 2*(1-stats::pnorm(abs(test.stat[i])))
      crit <- stats::qnorm(0.975)
      asymp.CI.lower[i] <- fit.dat$betas[i]-crit*SE
      asymp.CI.higher[i] <- fit.dat$betas[i]+crit*SE
    }
  }

  fit.dat <- data.frame(fit.dat,
                        Z=test.stat,
                        asymp.CI.lower=asymp.CI.lower,
                        asymp.CI.higher=asymp.CI.higher,
                        p.vals=p.vals)
  if(quasi){
    names(fit.dat)[2] <- "t"
    if(disp.ratio>0.9 & disp.ratio<1.1){
      warning("Dispersion ratio near 1, quasi unecessitated")
    }
  }

  return(fit.dat)

}
