#' Generalized Linear Model: Negative Binomial Regression
#'
#' \code{glm_negb} is used to fit the negative binomial generalized linear model, specified
#' by \code{formula}, a symbolic representation of the linear predictor.
#'
#' @param data A dataframe object that contains vectors spesified by the symbolic
#' representation of \code{formula}
#' @param formula A symbolic represenation of the linear predictor used to fit the
#' model
#' @param offsetparm String representation of a dataframe column which applies an offset
#' to the linear predictor
#'
#' @details
#' A negative binomial generalized linear model is a statistical model that fits count-based
#' data by utilization of an IRWLS algorithm that makes proper adjustments to \eqn{\beta}
#' until convergence.
#'
#' A dispersion parameter \eqn{\theta} is also approximated within the algorithmic fit as
#' a dispersion control. It's typically fitted towards handling over-dispersion, but may
#' handle some cases of under-dispersion.
#'
#' NOTE: THE \code{offsetparm}'s offset is logged before being fit, be cautious when entering offset
#' data
#'
#' @returns
#' \code{glm_negb} returns an S3 object of class \code{"glm_negb"}, it is not inherited from
#' classes \code{"lm"} or \code{"glm"}. Functions specified for these packages will not work
#' with a \code{"glm_negb"} object.
#'
#' The function \code{\link{interpret}} can be used to obtain and print a summary of
#' results.
#'
#' The \code{"$"} syntax can be used to extract various useful features/information
#' of the output values from the initial model fit. Many of the features included in
#' the \code{"glm_pois"} object are layered in lists, any computations utlilizing these
#' features should be unlisted first using \code{\link{unlist}}.
#'
#' An object of class \code{"glm_negb"} contains at least the following components:
#'
#' \item{coefficients}{List containing estimated coefficients and their standard errors}
#' \item{residuals}{Working residuals computed in the final iteration of the IRWLS fit}
#' \item{summary}{List containing statistical inference and features of model fit}
#' \item{fitted.values}{Fitted mean values (i.e., \eqn{\lambda}) obtained by exponentiation of the linear predictor}
#' \item{df.residual}{Degrees of freedom for residuals}
#' \item{call}{The matched call}
#' \item{terms}{The \code{\link{terms}} object used}
#' \item{model}{List containing the response and design matrices}
#' \item{theta}{\eqn{\theta} approximated as dispersion factor}
#'
#' @author
#' Implementation of \code{glm_negb} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' Motor Trend Magazine (1974).
#' Motor Trend Road Tests.
#'
#' @examples
#' ## Example using mtcars data set (Motor Trend, 1974)
#' utils::data(mtcars)
#'
#' mod <- glm_negb(data=mtcars, hp~factor(gear)+factor(carb)+disp)
#'
#' ## Using an offset
#' mod <- glm_negb(data=mtcars, hp~factor(gear)+factor(carb)+disp, offsetparm="wt")
#'
#' ##Extracting features/information
#' interpret(mod)
#' betas <- mod$coefficients$betas
#' std.err <- mod$coefficients$std.error
#' disp <- mod$dispersion
#'
#' @importFrom stats model.frame model.response model.matrix pnorm qnorm dnbinom optim
#' @export

glm_negb <- function(data,
                     formula,
                     offsetparm = NULL){
  if(any(is.na(data)) | any(is.null(data))){
    warning("NAs or Nulls in data set, NAs or Nulls ignored.")
  }

  #Parameter initliazations
  ##########################
  par <- model.frame(formula, data=data)
  y <- model.response(par)
  X <- model.matrix(formula, data=par)
  betas <- matrix(0, nrow=ncol(X),ncol=1)
  theta <- 1

  if(is.null(offsetparm)){
    offset <- 0
  } else{
    offset <- log(data[[offsetparm]])
  }

  #MLE estimation of theta
  ##########################
  theta.MLE <- function(par, curr.mu, y, neg=F){
    ll <- sum(log(dnbinom(y, size=par, mu=curr.mu)))
    return(ifelse(neg, -ll, ll))
  }

  #IWLS algorithm model fit
  ##########################
  maxrep <- 1000
  i <- 1
  repeat{
    eta <- offset + X %*% betas
    eta.lin <- X %*% betas
    pred.means <- exp(pmin(eta, 700))
    theta <- optim(par = theta,
                   fn=theta.MLE,
                   curr.mu=pred.means,
                   method = "Brent",
                   neg=T,
                   y=y, lower=0, upper=1000)$par
    tXW <- t(X * as.vector(pred.means^2/(pred.means^2/theta+pred.means)))
    tXWX <- tXW %*% X
    z <- (eta.lin) + (y-pred.means)/pred.means
    tXWz <- tXW %*% z
    betas.new <- solve(tXWX, tXWz)
    ss <- sum((betas.new-betas)**2)
    betas <- betas.new
    if(ss < 1e-6){
      std.error <- sqrt(diag(solve(tXWX)))
      break
    } else{
      if(i == maxrep){
        std.error <- sqrt(diag(solve(tXWX)))
        warning("Maximum number of iterations reached.")
        break
      } else{
        i <- i+1
      }
    }
  }

  fit.dat <- list(coefficients=list(betas=as.vector(betas),
                                    std.error=std.error))

  # Coefficient Pvals and Confidence Intervals
  #############################################
  test.stat <- rep(NA, times = length(fit.dat$coefficients$betas))
  p.vals <- rep(NA, times = length(fit.dat$coefficients$betas))
  asymp.CI.lower <- rep(NA, times = length(fit.dat$coefficients$betas))
  asymp.CI.higher <- rep(NA, times = length(fit.dat$coefficients$betas))
  for(i in 1:length(fit.dat$coefficients$betas)){
    SE <- std.error[i]
    test.stat[i] <- fit.dat$coefficients$betas[i]/SE
    p.vals[i] <- 2*(1-pnorm(abs(test.stat[i])))
    crit <- qnorm(0.975)
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
                    theta = theta))
  names(fit.dat$coefficients$betas) <- colnames(X)
  names(fit.dat$coefficients$std.error) <- colnames(X)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)

  if (theta>=999){
    warning("Theta diverges, perhaps use a Poisson regression model?")
  }

  class(fit.dat) <- "glm_negb"
  return(fit.dat)

}
