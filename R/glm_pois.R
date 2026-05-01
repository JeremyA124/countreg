#' Generalized Linear Model: Poisson Regression
#'
#' \code{glm_pois} is used to fit the poisson generalized linear model, specified
#' by \code{formula}, a symbolic representation of the linear predictor.
#'
#' @param data A dataframe object that contains vectors spesified by the symbolic
#' representation of \code{formula}
#' @param formula A symbolic represenation of the linear predictor used to fit the
#' model
#' @param quasi T or F if a quasi dispersion parameter is added to the model.
#' @param offsetparm String representation of a dataframe column which applies an offset
#' to the linear predictor
#'
#' @details
#' A poisson generalized linear model is a statistical model that fits count-based
#' data by utilization of an IRWLS algorithm that makes proper adjustments to \eqn{\beta}
#' until convergence.
#'
#' A quasi model can also be specified by the \code{quasi} parameter which fits the
#' standard poisson model with a dispersion parameter, \eqn{\phi}, that's added into
#' the variance structure and subsequent computations which utilize variance.
#' \eqn{\phi < 1} indicates under dispersion while \eqn{\phi > 1} indicates overdispersion.
#'
#' NOTE: THE \code{offsetparm}'s offset is logged before being fit, be cautious when entering offset
#' data
#'
#' @returns
#' \code{glm_pois} returns an S3 object of class \code{"glm_pois"}, it is not inherited from
#' classes \code{"lm"} or \code{"glm"}. Functions specified for these packages will not work
#' with a \code{"glm_pois"} object.
#'
#' The function \code{\link{interpret}} can be used to obtain and print a summary of
#' results.
#'
#' The \code{"$"} syntax can be used to extract various useful features/information
#' of the output values from the initial model fit. Many of the features included in
#' the \code{"glm_pois"} object are layered in lists, any computations utlilizing these
#' features should be unlisted first using \code{\link{unlist}}.
#'
#' An object of class \code{"glm_pois"} contains at least the following components:
#'
#' \item{coefficients}{List containing estimated coefficients and their standard errors}
#' \item{residuals}{Working residuals computed in the final iteration of the IRWLS fit}
#' \item{summary}{List containing statistical inference and features of model fit}
#' \item{fitted.values}{Fitted mean values (i.e., \eqn{\lambda}) obtained by exponentiation of the linear predictor}
#' \item{df.residual}{Degrees of freedom for residuals}
#' \item{call}{The matched call}
#' \item{terms}{The \code{\link{terms}} object used}
#' \item{model}{List containing the response and design matrices}
#' \item{dispersion}{\eqn{\phi} computed with quasi algorithm}
#'
#' @author
#' Implementation of \code{glm_pois} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' U.S. Department of Agriculture (1940s).
#' Insecticide effectiveness experiment data.
#'
#' @examples
#' ## Basic fit
#' y <- rpois(10, 2)
#' x1 <- rnorm(10)
#' x2 <- rnorm(10)
#' df <- data.frame(y=y,x1=x1, x2=x2)
#' mod <- glm_pois(data=df, y~x1+x2)
#'
#' ## With quasi
#' mod <- glm_pois(data-df, y~x1+x2, quasi=T)
#'
#' # With offset
#' mod <- glm_pois(data=df, y~x1+x2, offsetparm="x2")
#'
#' ## Example with InsectSprays dataset (U.S. Department of Agriculture insecticide experiment, 1940s)
#' data(InsectSprays)
#' utils::data(InsectSprays)
#' mod2 <- glm_pois(data=df, count~factor(spray))
#'
#' ##Extracting features/information
#' interpret(mod)
#' interpret(mod2)
#' betas <- mod$coefficients$betas
#' std.err <- mod$coefficients$std.error
#' disp <- mod$dispersion
#'
#' @importFrom stats model.frame model.response model.matrix pnorm qnorm qt pt
#' @export

glm_pois <- function(data,
                     formula,
                     quasi = F,
                     offsetparm = NULL){

  if(any(is.na(data)) | any(is.null(data))){
    warning("NAs or Nulls in data set, NAs or Nulls ignored.")
  }

  #Parameter initialization
  ##########################
  par <- model.frame(formula, data=data)
  y <- model.response(par)
  X <- model.matrix(formula, data=par)
  betas <- matrix(0, nrow=ncol(X), ncol=1)

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
    eta <- offset + X %*% betas
    eta.lin <- X %*% betas
    pred.means <- exp(pmin(eta, 700))
    tXW <- t(X * as.vector(pred.means))
    tXWX <- tXW %*% X
    z <- (eta.lin)+(y-pred.means)/pred.means
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
    } else {
      if(i == maxrep){
        std.error <- sqrt(diag(solve(tXWX)))
        warning("Maximum number of iterations reached.")
        break
      } else {
        i <- i+1
      }
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
