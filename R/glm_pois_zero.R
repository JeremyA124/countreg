#' Generalized Linear Model: Zero-Inflated Poisson Regression
#'
#' \code{glm_pois_zero} is used to fit the zero-inflated poisson generalized linear model, specified
#' by \code{formula}, a symbolic representation of the linear predictor.
#'
#' @param data A dataframe object that contains vectors specified by the symbolic
#' representation of \code{formula}
#' @param formula.pois A symbolic representation of the linear predictor used to fit the
#' poisson model
#' @param formula.log A symbolic representation of the linear predictor used to fit the
#' logistic model
#' @param offsetparm String representation of a dataframe column which applies an offset
#' to the linear predictor on __BOTH__ models.
#'
#' @details
#' A zero-inflated poisson generalized linear model is a statistical model that fits count-based
#' data by utilization of an IRWLS algorithm that makes proper adjustments to \eqn{\beta} of the
#' poisson model and \eqn{\alpha} of the logistic model until convergence.
#'
#' The logistic model is spesified by:
#' \deqn{log(\frac{P(y_i=0)}{1-P(y_i=0)})=\alpha_0+\alpha_1 x_1 + \cdots + \alpha_{n-1} x_{n-1}}
#' Which models the probablity of observing structural zeros.
#'
#' NOTE: THE \code{offsetparm}'s offset is logged before being fit, be cautious when entering offset
#' data
#'
#' @returns
#' \code{glm_pois_zero} returns an S3 object of class \code{"glm_pois_zero"}, it is not inherited from
#' classes \code{"lm"} or \code{"glm"}. Functions specified for these packages will not work
#' with a \code{"glm_pois_zero"} object.
#'
#' The function \code{\link{interpret}} can be used to obtain and print a summary of
#' results.
#'
#' The \code{"$"} syntax can be used to extract various useful features/information
#' of the output values from the initial model fit. Many of the features included in
#' the \code{"glm_pois_zero"} object are layered in lists, any computations utlilizing these
#' features should be unlisted first using \code{\link{unlist}}.
#'
#' An object of class \code{"glm_pois_zero"} contains at least the following components:
#'
#' \item{coefficients}{List containing estimated coefficients for poisson/logistic models
#'  and their standard errors}
#' \item{residuals}{Working residuals for the poisson model computed in the final iteration of the IRWLS fit}
#' \item{summary}{List containing statistical inference and features of poisson/logistic model fit}
#' \item{fitted.values}{Fitted mean values (i.e., \eqn{\lambda}) obtained by exponentiation of the poisson linear predictor}
#' \item{pred.zero}{Fitted values (i.e., \eqn{\alpha}) obtained by exponentiation of the logistic linear predictor}
#' \item{df.residual}{Degrees of freedom for residuals}
#' \item{call}{The matched call}
#' \item{terms}{The \code{\link{terms}} object used}
#' \item{model}{List containing the response and design matrices}
#'
#' @author
#' Implementation of \code{glm_pois_zero} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' Long, J. S. (1997).
#' Regression Models for Categorical and Limited Dependent Variables.
#' Sage Publications.
#'
#' @examples
#' ## Example using biochemists dataset from pscl (Campbell & Mahon, 1974)
#' #' utils::data(bioChemists, package = "pscl")
#' mod <- glm_pois_zero(data=bioChemists, kid5~phd, kid5~phd)
#'
#' ## With offset
#' mod <- glm_pois_zero(data=df, kid5~art, ment~art, offsetparm="phd")
#'
#' ##Extracting features/information
#' interpret(mod)
#' betas <- mod$coefficients$count
#' std.err <- mod$coefficients$std.error.count
#' zeros <- mod$pred.zero
#'
#' @importFrom stats model.frame model.response model.matrix model.offset terms qnorm pnorm
#' @export

glm_pois_zero <- function(data,
                          formula.pois,
                          formula.log,
                          offsetparm = NULL){

  if(any(is.na(data)) | any(is.null(data))){
    warning("NAs or Nulls in data set, NAs or Nulls ignored.")
  }

  #Parameter initliazations
  ##########################
  par1 <- model.frame(formula.pois, data=data)
  par2 <- model.frame(formula.log, data=data)
  y <- model.response(par1)
  X.pois <- model.matrix(formula.pois, data=par1)
  X.logit <- model.matrix(formula.log, data=par2)
  betas <- matrix(0, nrow=ncol(X.pois), ncol=1)
  alphas <- matrix(0, nrow=ncol(X.logit), ncol=1)
  pred.means <- 1

  if(is.null(offsetparm)){
    offset <- 0
  } else{
    offset <- log(data[[offsetparm]])
  }

  #IWLS algorithm model fit
  ##########################
  maxrep <- 1000
  i <- 1
  repeat{
    #Alpha Part
    ############
    eta.logit <- offset + X.logit %*% alphas
    eta.logit.lin <- X.logit %*% alphas
    pred.logs <- exp(pmin(eta.logit, 700))
    pred.zeros <- pred.logs/(1+pred.logs)
    delta <- (y == 0) * pred.zeros/(pred.zeros + (1-pred.zeros)*exp(-pred.means))
    tXW.logit <- t(X.logit * as.vector(pred.zeros*(1-pred.zeros)))
    tXWX.logit <- tXW.logit %*% X.logit
    z.logit <- (eta.logit.lin) + (delta-pred.zeros)/(pred.zeros*(1-pred.zeros))
    tXWz.logit <- tXW.logit %*% z.logit
    new.alphas <- solve(tXWX.logit, tXWz.logit)
    ss1 <- sum((new.alphas-alphas)^2)
    alphas <- new.alphas
    #Poisson Part
    ##############
    eta <- offset + X.pois %*% betas
    eta.lin <- X.pois %*% betas
    pred.means <- exp(pmin(eta, 700))
    tXW.pois <- t(X.pois * as.vector((1-delta)*pred.means))
    tXWX.pois <- tXW.pois %*% X.pois
    z <- (eta.lin)+(y-pred.means)/(pred.means)
    tXWz.pois <- tXW.pois %*% z
    betas.new <- solve(tXWX.pois, tXWz.pois)
    ss2 <- sum((betas.new-betas)**2)
    betas <- betas.new
    if(ss1 < 1e-6 && ss2 < 1e-6){
      std.error.pois <- sqrt(diag(solve(tXWX.pois)))
      std.error.log <- sqrt(diag(solve(tXWX.logit)))
      break
    } else{
      if(i == maxrep){
        std.error.pois <- sqrt(diag(solve(tXWX.pois)))
        std.error.log <- sqrt(diag(solve(tXWX.logit)))
        warning("Maximum number of iterations reached.")
        break
      } else{
        i <- i+1
      }
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
    p.vals[i] <- 2*(1-pnorm(abs(test.stat[i])))
    crit <- qnorm(0.975)
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
    p.vals.alpha[i] <- 2*(1-pnorm(abs(test.stat.alpha[i])))
    crit <- qnorm(0.975)
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
