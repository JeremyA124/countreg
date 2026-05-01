#' Generalized Linear Model: Generalized Poisson (GP2)
#'
#' \code{glm_pois_GP2} is used to fit the Generalized poisson (GP2) generalized linear model, specified
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
#' A generalized poisson (GP2) generalized linear model is a statistical model that fits count-based
#' data by utilization of an IRWLS algorithm that makes proper adjustments to \eqn{\beta}
#' until convergence.
#'
#' A dispersion parameter \eqn{\alpha} is also approximated within the algorithmic fit as
#' a dispersion control. It's typically fitted to handle all dispersion types.
#'
#' The following probability model is utilized through \code{\link{VGAM}}'s \code{\link{dgenpois2}}
#' function for the approximation of \eqn{\alpha}.
#' \deqn{P_y(Y_i=y_i)=\frac{\lambda(\lambda+\alpha y)^{y-1}}{y_i!}\mathrm{exp}^{-(\lambda+\alpha y)}}
#'
#' NOTE: THE \code{offsetparm}'s offset is logged before being fit, be cautious when entering offset
#' data
#'
#' @returns
#' \code{glm_pois_GP2} returns an S3 object of class \code{"glm_pois_GP2"}, it is not inherited from
#' classes \code{"lm"} or \code{"glm"}. Functions specified for these packages will not work
#' with a \code{"glm_pois_GP2"} object.
#'
#' The function \code{\link{interpret}} can be used to obtain and print a summary of
#' results.
#'
#' The \code{"$"} syntax can be used to extract various useful features/information
#' of the output values from the initial model fit. Many of the features included in
#' the \code{"glm_pois_GP2"} object are layered in lists, any computations utlilizing these
#' features should be unlisted first using \code{\link{unlist}}.
#'
#' An object of class \code{"glm_pois_GP2"} contains at least the following components:
#'
#' \item{coefficients}{List containing estimated coefficients and their standard errors}
#' \item{residuals}{Working residuals computed in the final iteration of the IRWLS fit}
#' \item{summary}{List containing statistical inference and features of model fit}
#' \item{fitted.values}{Fitted mean values (i.e., \eqn{\lambda}) obtained by exponentiation of the linear predictor}
#' \item{df.residual}{Degrees of freedom for residuals}
#' \item{call}{The matched call}
#' \item{terms}{The \code{\link{terms}} object used}
#' \item{model}{List containing the response and design matrices}
#' \item{alpha}{\eqn{\alpha} approximated as dispersion factor}
#'
#' @author
#' Implementation of \code{glm_pois_GP2} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' Bailey, N. T. J. (1953).
#' Statistical Methods in Biology.
#' Cambridge University Press.
#'
#' @examples
#' ## Example using warpbreaks data set (British Textile Research data, early loom breakage experiments)
#' #' utils::data(warpbreaks)
#'
#' mod <- glm_pois_GP2(data=warpbreaks, breaks~factor(wool)+factor(tension))
#'
#' ## Using an "example" offset
#' off <- rpois(length(warp), 5)
#' warp <- cbind(warpbreaks, off)
#' mod <- glm_pois_GP2(data=warp, breaks~factor(wool)+factor(tension), offsetparm="off")
#'
#' ##Extracting features/information
#' interpret(mod)
#' betas <- mod$coefficients$betas
#' std.err <- mod$coefficients$std.error
#' disp <- mod$dispersion
#'
#' @importFrom stats model.frame model.response model.matrix model.offset terms optim dpois qnorm pnorm
#' @importFrom VGAM dgenpois2
#' @export


glm_pois_GP2 <- function(data,
                         formula,
                         offsetparm = NULL){
  require(VGAM)

  if(any(is.na(data)) | any(is.null(data))){
    warning("NAs or Nulls in data set, NAs or Nulls ignored.")
  }

  alpha.MLE <- function(par,
                        y,
                        fits,
                        neg=F){
    ll <- sum(log(dgenpois2(y, meanpar=fits, disppar = par)))
    return(ifelse(neg, -ll, ll))
  }

  #Parameter initliazations
  ##########################
  par <- model.frame(formula, data=data)
  y <- model.response(par)
  X <- model.matrix(formula, data=par)
  betas <- matrix(0, nrow=ncol(X),ncol=1)
  alpha <- 1

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
    alpha <- optim(par = alpha,
                   fn=alpha.MLE,
                   fits=pred.means,
                   method = "Brent",
                   neg=T,
                   y=y, lower=1e-4, upper=1000)$par
    tXW <- t(X * as.vector(pred.means/(1+alpha*pred.means)^2))
    tXWX <- tXW %*% X
    z <- eta.lin+(y-pred.means)/(pred.means)
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
                    df.residuals = nrow(data)-length(betas)-1,
                    call = match.call(),
                    terms = terms(formula),
                    model = list(y=y, x=X),
                    alpha=alpha))
  names(fit.dat$coefficients$betas) <- colnames(X)
  names(fit.dat$coefficients$std.error) <- colnames(X)
  names(fit.dat$residuals) <- 1:length(fit.dat$residuals)
  names(fit.dat$fitted.values) <- 1:length(fit.dat$fitted.values)

  class(fit.dat) <- "glm_pois_GP2"
  return(fit.dat)
}
