#' Residual Diagnostics for Count Regression models
#'
#' \code{resdi_diag} is used to generate standardized pearson residuals and randomized
#' quantile residuals for evaluation and diagnostics of models in \code{\link{countmods}}.
#'
#' @param mod Models of classes built by \code{\link{countmods}}
#'
#' @details
#' \code{resid_diag} utilizes randomized quantile residuals to test for
#' goodness-of-fit and dispersion of residuals. The following algorithm is used to
#' generate randomized quantile residuals.
#' \deqn{a_i=P(Y_i=y_i-1)}
#' \deqn{b_i=P(Y_i=y_i)}
#' \deqn{u_i=\mathrm{random uniform}(a_i,b_i)}
#' \deqn{e_i=\Phi(u_i)}
#' Where \eqn{\Phi} denotes the CDF of the normal distribution.
#'
#' In addition, standardized pearson residuals are computed for models in the
#' poisson family (\code{\link{glm_pois}}, \code{\link{glm_pois_zero}}, \code{\link{glm_pois_GP2}}).
#' Standardized pearson residuals provide another metric for determining dispersion and validation
#' of model choice.
#'
#' @author
#' Implementation of \code{resid_diag} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' Long, J. S. (1997).
#' Regression Models for Categorical and Limited Dependent Variables.
#' Sage Publications.
#'
#' Deb, P., & Trivedi, P. K. (1997).
#' Demand for medical care by the elderly: A finite mixture approach.
#' Journal of Applied Econometrics, 12(3), 313--336.
#'
#' @examples
#' ## Example using biochemists dataset from pscl (Campbell & Mahon, 1974)
#' utils::data(bioChemists, package = "pscl")
#' mod <- glm_pois(data=bioChemists, kid5~phd)
#'
#' resid_diag(mod)
#'
#' ## Example using NMES1988 dataset from AER (Deb & Trivedi, 1997)
#' utils::data(NMES1988, package = "AER")
#' mod2 <- glm_negb_zero(data=NMES1988, visits~factor(health), visits~factor(health))
#'
#' resid_diag(mod2)
#'
#' @importFrom stats runif qnorm ppois pnbinom
#' @importFrom ggplot2 ggplot geom_smooth geom_point geom_hline theme_bw labs stat_qq stat_qq_line ggtitle
#' @importFrom VGAM pgenpois2 dgenpois2
#' @export

resid_diag <- function(mod){
  require(patchwork)
  require(ggplot2)
  require(VGAM)

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

  ris.compute <- function(mod.type, y, fits, pis=0, alpha=0){
    if(mod.type=="glm_pois"){
      ris <- ((y-fits)/sqrt(fits))
    } else if(mod.type=="glm_pois_zero"){
      mu <- (1-pis)*fits
      var <- sqrt(fits*(1-pis)*(1+pis*fits))
      ris <- (y-mu)/var
    } else{
      ris <- (y-fits)/(sqrt(fits*(1+fits*alpha)^2))
    }

    return(ris^2)
  }

  disp.compute <- function(mod.type, y, fits, df, theta, pis, alpha){
    if(mod.type %in% c("glm_pois", "glm_negb")){
      disp <- sum(((y-fits)^2/(fits^2/theta+fits)))/df
    } else if(mod.type %in% c("glm_pois_zero", "glm_negb_zero")){
      sigma <- sqrt((1-pis)*fits*(1+fits/theta+pis*fits))
      disp <- sum(((y-fits)^2/sigma^2))/df
    } else{
      disp <- sum((y-fits)^2/(fits*(1+fits*alpha)^2))/df
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
      ai <- pgenpois2(counts[i]-1, meanpar=fit.vals[i], disppar=alpha)
      bi <- pgenpois2(counts[i], meanpar=fit.vals[i], disppar=alpha)
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
