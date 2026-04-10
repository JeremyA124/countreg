glm_negb <- function(data,
                     formula,
                     offset = log(1)){
  #Parameter initliazations
  ##########################
  par <- stats::model.frame(formula, data=data)
  y <- stats::model.response(par)
  X <- stats::model.matrix(formula, data=par)
  betas <- matrix(0, nrow=ncol(X),ncol=1)
  theta <- 1

  #MLE estimation of theta
  ##########################
  theta.MLE <- function(par, curr.mu, y, neg=F){
    ll <- sum(log(stats::dnbinom(y, size=par, mu=curr.mu)))
    return(ifelse(neg, -ll, ll))
  }

  #IWLS algorithm model fit
  ##########################
  repeat{
    eta <- offset + X %*% betas
    pred.means <- exp(eta)
    theta <- stats::optim(par = theta,
                          fn=theta.MLE,
                          curr.mu=pred.means,
                          method = "Brent",
                          neg=T,
                          y=y, lower=0, upper=1000)$par
    W <- diag(as.vector(pred.means**2/((1/theta)*pred.means**2+pred.means)))
    z <- eta+(y-pred.means)/((1/theta)*pred.means**2+pred.means)
    inv.Hess <- solve(t(X)%*%W%*%X)
    betas.new <- inv.Hess%*%t(X)%*%W%*%z
    ss <- sum((betas.new-betas)**2)
    betas <- betas.new
    std.error <- sqrt(diag(inv.Hess))
    if(ss < 1e-6){
      break
    }
  }

  fit.dat <- data.frame(betas=betas,std.error=std.error)
  print(theta)

  # Coefficient Pvals and Confidence Intervals
  #############################################
  Zs <- rep(NA, times = nrow(fit.dat))
  p.vals <- rep(NA, times = nrow(fit.dat))
  asymp.CI.lower <- rep(NA, times = nrow(fit.dat))
  asymp.CI.higher <- rep(NA, times = nrow(fit.dat))
  for(i in 1:nrow(fit.dat)){
    SE <- sqrt(inv.Hess[i,i])
    Zs[i] <- fit.dat$betas[i]/SE
    p.vals[i] <- 2*(1-stats::pnorm(abs(Zs[i])))
    crit <- stats::qnorm(0.975)
    asymp.CI.lower[i] <- fit.dat$betas[i]-crit*SE
    asymp.CI.higher[i] <- fit.dat$betas[i]+crit*SE
  }

  fit.dat <- data.frame(fit.dat,
                        Z=Zs,
                        asymp.CI.lower=asymp.CI.lower,
                        asymp.CI.higher=asymp.CI.higher,
                        p.vals=p.vals)

  return(fit.dat)

}
