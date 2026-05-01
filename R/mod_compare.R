#' Comparison of Model Fit
#'
#' \code{mod_compare} is used to compare models from \code{\link{countmods}} for
#' goodness-of-fit using various model fit metrics
#'
#' @param mod1 First Model of classes built by \code{\link{countmods}}
#' @param mod2 Second Model of classes built by \code{\link{countmods}}
#' @param print.results T or F depending on user prefrence for result display. It
#' is recommended to keep this parameter as T.
#'
#' @details
#' \code{mod_compare} utilizes AIC, BIC, Likelihood, and McFadden's \eqn{R^2} as
#' metric of goodness-of-fit model comparisons. This may also be used as a generic
#' function for computation of these statistics by passing the same model to \code{mod1}
#' and \code{mod2}.
#'
#' A likelihood ratio test is also carried out. Twice the difference of likelihoods
#' is compared to the Chi-squared distribution.
#'
#' The user may also have the option of returning just the statistics related to
#' standard model comparison (AIC, BIC, Likelihood, McFadden) without displaying results
#' using \code{print.results}. Evaluating \code{print.results} as FALSE will compute the
#' likelihood ratio test, but will not display results.
#'
#' @returns List object that contains AICs, BICs, Likelihoods, and McFadden's \eqn{R^2}
#' for both models
#'
#' @author
#' Implementation of \code{mod_compare} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#'
#' Deb, P., & Trivedi, P. K. (1997).
#' Demand for medical care by the elderly: A finite mixture approach.
#' Journal of Applied Econometrics, 12(3), 313--336.
#'
#' @examples
#' ## Example using NMES1988 dataset from AER (Deb & Trivedi, 1997)
#' utils::data(NMES1988, package = "AER")
#' mod1 <- glm_pois(data=NMES1988, visits~factor(chronic)+factor(health))
#' mod2 <- glm_negb(data=NMES1988, visits~factor(chronic)+factor(health))
#'
#' comp.stats <- mod_compare(mod1, mod2)
#'
#' ## Or as a generic function for statistics
#' gen.stats <- mod_compare(mod1, mod1, print.results=F)
#'
#' @importFrom stats dpois dnbinom pchisq dchisq update formula terms model.frame model.response model.matrix model.offset
#' @importFrom VGAM pgenpois2

mod_compare <- function(mod1, mod2, print.results=T){
  require(VGAM)

  if (!(class(mod1) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero", "glm_pois_GP2")) |
      !(class(mod2) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero", "glm_pois_GP2"))){
    stop("Model class not supported")
  }

  dzpois <- function(y, lambda, pi) {
    probs <- numeric(length(y))
    pi <- as.numeric(pi)
    probs[y == 0] <- pi[y == 0] + (1 - pi[y == 0]) * exp(-lambda[y == 0])
    probs[y > 0] <- (1 - pi[y > 0]) * dpois(y[y > 0], lambda[y > 0])
    return(probs)
  }

  dznbinom <- function(y, mu, theta, pi) {
    probs <- numeric(length(y))
    pi <- as.numeric(pi)
    probs[y == 0] <- pi[y == 0] + (1 - pi[y == 0]) * (theta / (mu[y == 0] + theta))^theta
    probs[y > 0] <- (1 - pi[y > 0]) * dnbinom(y[y > 0], size = theta, mu = mu[y > 0])
    return(probs)
  }

  compute.logliki <- function(mod.class, y, fits, zeros, theta, alpha){
    if(mod.class=="glm_pois"){
      liklihood <- sum(dpois(y,
                             lambda = fits,
                             log = T))
    } else if(mod.class=="glm_negb"){
      liklihood <- sum(dnbinom(y,
                               size = theta,
                               mu = fits,
                               log = T))
    } else if(mod.class=="glm_pois_zero"){
      liklihood <- sum(log(dzpois(y,
                                  lambda = fits,
                                  pi = zeros)))
    } else if(mod.class=="glm_negb_zero"){
      liklihood <- sum(log(dznbinom(y,
                                    mu = fits,
                                    theta = theta,
                                    pi = zeros)))
    } else if(mod.class=="glm_pois_GP2"){
      liklihood <- sum(log(pgenpois2(y,
                                     meanpar = fits,
                                     disppar = alpha)))
    }

    return(liklihood)
  }

  compute.compstats <- function(mod){

    theta <- ifelse(class(mod)%in%c("glm_negb","glm_negb_zero"), mod$theta, NA)
    alpha <- ifelse(class(mod) == "glm_pois_GP2", mod$alpha, NA)
    if(class(mod) %in% c("glm_pois_zero", "glm_negb_zero")){pi <- mod$pred.zero}
      else{pi <- NA}
    has.disp <- ifelse(is.null(mod$dispersion), 0, 1)

    liklihood <- compute.logliki(class(mod),
                                 mod$model$y,
                                 unlist(mod$fitted.values),
                                 pi,
                                 theta,
                                 alpha)

    if(class(mod) %in% c("glm_pois", "glm_negb")){
      aic <- 2*(length(mod$coefficients$betas)
                +as.integer(!is.na(theta))+has.disp)-2*liklihood
    } else if(class(mod) %in% c("glm_pois_zero", "glm_negb_zero")){
      aic <- 2*(length(mod$coefficients$count)
                +length(mod$coefficients$zero)
                +as.integer(!is.na(theta)))-2*liklihood
    } else if(class(mod) == "glm_pois_GP2"){
      aic <- 2*(length(mod$coefficients$count)
                +as.integer(!is.na(theta))+1)-2*liklihood
    }

    if(class(mod) %in% c("glm_pois", "glm_negb")){
      bic <- (length(mod$coefficients$betas)
              +as.integer(!is.na(theta))
              +has.disp)*log(length(mod$model$y))-2*liklihood
    } else if(class(mod) %in% c("glm_pois_zero", "glm_negb_zero")){
      bic <- (length(mod$coefficients$count)
              +length(mod$coefficients$zero)
              +as.integer(!is.na(theta))) *log(length(mod$model$y))-2*liklihood
    } else if(class(mod) == "glm_pois_GP2"){
      bic <- (length(mod$coefficients$betas)
              +as.integer(!is.na(theta))+1)*log(length(mod$model$y))-2*liklihood
    }

    if(class(mod) %in% c("glm_pois_zero", "glm_negb_zero")){
      null_form_count <- as.formula(paste(all.vars(formula(mod$terms$count))[1], "~ 1"))
      null_form_zero <- as.formula(paste(all.vars(formula(mod$terms$zero))[1], "~ 1"))
    }

    if(class(mod) %in% c("glm_pois", "glm_negb", "glm_pois_GP2")){
      null_mod <- update(mod, .~1)
    } else if(class(mod) == "glm_pois_zero"){
      null_mod <- update(mod, formula.pois = null_form_count , formula.log = null_form_zero)
    } else if(class(mod) == "glm_negb_zero"){
      null_mod <- update(mod, formula.negb = null_form_count, formula.log = null_form_zero)
    }

    null_theta <- ifelse(class(null_mod)%in%c("glm_negb","glm_negb_zero"), null_mod$theta, NA)
    null_alpha <- ifelse(class(null_mod) == "glm_pois_GP2", null_mod$alpha, NA)
    if(class(null_mod) %in% c("glm_pois_zero", "glm_negb_zero")){null_pi <- null_mod$pred.zero}
      else{null_pi <- NA}

    null_liklihood <- compute.logliki(class(null_mod),
                                      null_mod$model$y,
                                      unlist(null_mod$fitted.values),
                                      null_pi,
                                      null_theta,
                                      null_alpha)

    mcfadden <- 1-(liklihood/null_liklihood)

    return(c(round(aic,2), round(bic,2), round(mcfadden,2), round(liklihood,2)))

  }

  mod1.stats <- compute.compstats(mod1)
  mod2.stats <- compute.compstats(mod2)

  if(class(mod1) %in% c("glm_negb","glm_negb_zero")){
    q <- 2*(mod1.stats[4]-mod2.stats[4])
  } else{
    q <- 2*(mod2.stats[4]-mod1.stats[4])
  }

  chi <- dchisq(q, 1)
  p_val <- pchisq(q, 1)

  if(print.results){
    cat("---------------------------\n")
    cat("Model 1:\n")
    cat("AIC:", mod1.stats[1], "\n")
    cat("BIC:", mod1.stats[2], "\n")
    cat("R^2 (McFadden):", mod1.stats[3], "\n")
    cat("logLiki:", mod1.stats[4], "\n")
    cat("---------------------------\n")
    cat("Model 2:\n")
    cat("AIC:", mod2.stats[1], "\n")
    cat("BIC:", mod2.stats[2], "\n")
    cat("R^2 (McFadden):", mod2.stats[3], "\n")
    cat("logLiki:", mod2.stats[4], "\n")
    cat("---------------------------\n")
    cat("Liklihood Ratio Test:\n")
    cat("\u03c7\u00b2:",chi, "\n")
    cat("p-val:", p_val, "\n")
    cat("---------------------------\n")
  }

  return(list(mod1=mod1.stats,
              mod2=mod2.stats))

}
