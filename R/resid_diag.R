resid_diag <- function(mod){
  require(patchwork)
  require(ggplot2)

  if (!(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero", "glm_pois_GP2"))){
    stop("Model class not supported")
  }

  pznbinom <- function(q, mu, theta, pi) {
    probs <- numeric(length(q))
    if(length(pi) == 1) pi <- rep(pi, length(q))
    if(length(mu) == 1) mu <- rep(mu, length(q))
    probs[q < 0] <- 0
    idx <- q >= 0
    probs[idx] <- pi[idx] + (1 - pi[idx]) * pnbinom(q[idx], size = theta, mu = mu[idx])

    return(probs)
  }

  pzpois <- function(q, lambda, pi) {
    probs <- numeric(length(q))
    probs[q < 0] <- 0
    idx <- q >= 0
    probs[idx] <- pi[idx] + (1 - pi[idx]) * ppois(q[idx], lambda[idx])

    return(probs)
  }

  ppoisgp2 <- function(q, lambda, alpha){
    probs <- numeric(length(q))
    for(j in 1:length(q)){
      for(i in 0:q[j]){
        comp1 <- lambda/(1+alpha*lambda)
        comp2 <- (lambda+alpha*i)^(i-1)/factorial(i)
        comp3 <- exp((-lambda*(1+alpha*i))/(1+alpha*lambda))
        probs[j] <- probs[j] + comp1*comp2*comp3
      }
    }
    return(probs)
  }

  ris.compute <- function(mod.type, y, fits, pis=0, alpha=0){
    if(mod.type=="glm_pois"){
      ris <- ((y-fit.vals)/sqrt(fit.vals))
    } else if(mod.type=="glm_pois_zero"){
      mu <- (1-pis)*fits
      var <- sqrt(fits*(1-pis)*(1+pis*fits))
      ris <- y-mu/var
    } else{
      ris <- y-fits/(sqrt(fits)*(1+fits*alpha))
    }

    return(ris^2)
  }

  disp.compute <- function(mod.type, y, fits, df, theta, pis, alpha){
    if(mod.type %in% c("glm_pois", "glm_negb")){
      disp <- sum(((y-fits)/sqrt(fits^2/theta+fits))^2)/df
    } else if(mod.type %in% c("glm_pois_zero", "glm_negb_zero")){
      sigma <- sqrt((1-pis)*fits*(1+fits/theta+pis*fits))
      disp <- sum(((y-fits)/sigma)^2)/df
    } else{
      disp <- sum((y-fits)^2/(sqrt(fits)*(1+fits*alpha)))/df
    }

    return(disp)
  }

  # Pearson Residuals Calculations and Plot
  ##########################################
  counts <- mod$model$y
  if(class(mod) %in% c("glm_pois", "glm_negb")){pis <- 0}
     else{pis <- unlist(mod$pred.zero)}
  fit.vals <- as.vector(unlist(mod$fitted.values))
  theta <- ifelse(class(mod)%in%c("glm_negb", "glm_negb_zero"), mod$theta, 1e100)
  alpha <- ifelse(class(mod)=="glm_pois_GP2", mod$alpha, 1e100)
  ris <- ris.compute(class(mod), counts, fit.vals, pis, alpha)

  p1 <- ggplot(data = data.frame(fit.val=fit.vals,
                                 ri=ris)) +
    geom_smooth(aes(x=fit.val, y=ri),
              color = 'red',
              linetype='dotted') +
    geom_point(aes(x=fit.val, y=ri)) +
    geom_hline(yintercept = 1, linetype = 'dotted') +
    theme_bw() +
    labs(title = expression(Pearson~Residuals~over~mu*"'s"),
         y = expression(r[i]^2),
         x = expression(mu[i]))

  # Randomized Quantized Residuals Calculation and Plot
  ######################################################
  rqr <- rep(NA, times=length(fit.vals))
  for (i in 1:length(fit.vals)){
    if (class(mod) =="glm_pois"){
      ai <- ppois(counts[i]-1, lambda = fit.vals[i])
      bi <- ppois(counts[i], lambda = fit.vals[i])
    } else if(class(mod)== "glm_pois_zero"){
      ai <- pzpois(counts[i]-1, lambda = fit.vals[i], pi=pis[i])
      bi <- pzpois(counts[i], lambda = fit.vals[i], pi=pis[i])
    } else if(class(mod) == "glm_pois_GP2"){
      ai <- ppoisgp2(counts[i]-1, lambda = fit.vals[i], alpha=alpha)
      bi <- ppoisgp2(counts[i], lambda = fit.vals[i], alpha=alpha)
    } else if(class(mod) == "glm_negb_zero"){
      ai <- pznbinom(counts[i]-1, mu = fit.vals[i], theta=mod$theta, pi=pis[i])
      bi <- pznbinom(counts[i], mu = fit.vals[i], theta=mod$theta, pi=pis[i])
    } else if(class(mod) == "glm_negb"){
      ai <- pnbinom(counts[i]-1, size=theta, mu = fit.vals[i])
      bi <- pnbinom(counts[i], size=theta, mu = fit.vals[i])
    }
    ui <- ai * runif(1) + (bi-ai)
    ui <- max(min(ui, 1-10^(-6)), 10^(-6))
    rqr[i] <- qnorm(ui)
  }

  dispersion <- disp.compute(class(mod),
                             counts,
                             fit.vals,
                             mod$df.residuals,
                             theta=theta,
                             pis=pis,
                             alpha=alpha)
  p2 <- ggplot(data=data.frame(fit.val=fit.vals,
                               rqrs=rqr)) +
    geom_hline(yintercept=0, linetype="dotted")+
    geom_point(aes(x=fit.val, y=rqrs)) +
    theme_bw() +
    labs(title = "Randomized Quantized Residuals",
          x = bquote(mu),
          y = "RQR")
  p3 <- ggplot(data=data.frame(rqrs=rqr)) +
    stat_qq(aes(sample=rqrs)) +
    stat_qq_line(aes(sample=rqrs)) +
    theme_bw() +
    ggtitle(paste("Dispersion Ratio =", round(dispersion, 4)))

  if(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_pois_GP2")){
    print(p1+p2/p3)
  } else {
    print(p2+p3)
  }

  diag.message <- ifelse(dispersion>=2, "Overdispersed, try quasi or negative-binomial?",
                     ifelse(dispersion<=0.85, "Underdispersed, try quasi or proceed with caution", "Passed"))
  message(paste("Model fit check:",diag.message))
  message("NOTE: Choose different model structure if wedge shape appears in Pearson Residuals")
  message("(Poisson models only)")
}
