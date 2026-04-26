resid_diag <- function(mod){
  require(patchwork)
  require(ggplot2)

  if (!(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero"))){
    stop("Model class not supported")
  }

  pznbinom <- function(q, mu, theta, pi) {
    prob <- pi + (1 - pi) * pnbinom(q, size = theta, mu = mu)
    return(prob)
  }

  pzpois <- function(q, lambda, pi) {
    prob <- pi + (1 - pi) * ppois(q, lambda)
    return(prob)
  }

  ris.compute <- function(mod.type, y, fits, pis){
    if(mod.type=="glm_pois"){
      ris <- ((counts-fit.vals)/sqrt(fit.vals))
    } else {
      mu <- (1-pis)*fits
      sigma <- sqrt(fits*(1-pis)*(1+pis*fits))
      ris <- (y-mu)/sigma
    }

    return(ris^2)
  }

  disp.compute <- function(mod.type, y, fits, df, theta, pis){
    if(mod.type %in% c("glm_pois", "glm_negb")){
      disp <- sum(((y-fits)/sqrt(fits^2/theta+fits))^2)/df
    } else{
      sigma <- sqrt((1-pis)*fits*(1+fits/theta+pis*fits))
      disp <- sum(((y-fits)/sigma)^2)/df
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
  ris <- ris.compute(class(mod), counts, fit.vals, pis)

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
    } else if(class(mod) == "glm_negb_zero"){
      ai <- pznbinom(counts[i]-1, mu = fit.vals[i], theta=mod$theta, pi=pis[i])
      bi <- pznbinom(counts[i], mu = fit.vals[i], theta=mod$theta, pi=pis[i])
    } else {
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
                             pis=pis)
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

  if(class(mod) %in% c("glm_pois", "glm_pois_zero")){
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
