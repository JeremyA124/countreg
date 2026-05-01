#' Test for Zero-Inflation of Count Regression Models
#'
#' \code{zeroinf_test} is used to test zero inflation within the \code{\link{countmods}} package.
#'
#' @param mod Models of classes built by \code{\link{countmods}}
#'
#' @details
#' \code{zeroinf_test} tests models in \code{\link{countmods}} through the simulation of the
#' expected response and counting the number of zeros found in the data. A distribution is
#' built through the simulations of the zero counts and the emperical zero count is plotted
#' on top of the distribution.
#' A right tailed emperical test is used to gauge the need for a zero-inflated model.
#'
#' \code{zeroinf_test} provides the interpretation automatically.
#'
#' @author
#' Implementation of \code{zeroinf_test} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' Long, J. S. (1997).
#' Regression Models for Categorical and Limited Dependent Variables.
#' Sage Publications.
#'
#' @examples
#' ## Example using biochemists dataset from pscl (Campbell & Mahon, 1974)
#' utils::data(bioChemists, package = "pscl")
#' mod <- glm_pois(data=bioChemists, kid5~phd)
#'
#' zeroinf_test(mod)
#'
#' @importFrom stats rpois rnbinom
#' @importFrom ggplot2 ggplot geom_histogram geom_vline geom_hline theme_bw labs theme
#' @export

zeroinf_test <- function(mod){
  require(ggplot2)

  if (class(mod) %in% c("glm_pois_zero", "glm_negb_zero")){
    stop("This is a zero-inflated model!")
  }

  if (!(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero", "glm_pois_GP2"))){
    stop("Model class not supported")
  }

  # Zero Simulation
  ##################
  R <- 1000
  obsZeros <- sum(mod$model$y == 0)
  numZeros <- rep(NA, times=R)
  for(i in 1:length(numZeros)){
    if(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_pois_GP2")){
      simulated_y <- rpois(length(unlist(mod$fitted.values)),
                           lambda=unlist(mod$fitted.values))
    } else {
      simulated_y <- rnbinom(length(unlist(mod$fitted.values)),
                             size = mod$theta,
                             mu=unlist(mod$fitted.values))
    }
    numZeros[i] <- sum(simulated_y == 0)
  }

  meanZeros <- mean(numZeros <= obsZeros)

  # Diagnostics Plot
  ###################
  p1 <- ggplot(data=as.data.frame(numZeros)) +
    geom_histogram(aes(x = numZeros), color = 'black', binwidth = 1) +
    geom_vline(aes(xintercept = obsZeros,
                   color = "Observed Zero Count"),
               linetype = "dotted",
               linewidth = 1) +
    geom_hline(yintercept = 0) +
    theme_bw() +
    labs(title = "Diagnostics for Zero-Inflation",
         subtitle = ifelse(meanZeros>=0.975,
                           "Zero-Inflated Model recommended",
                           "No zero-inflation issue"),
         x = "Number of Zeros",
         y = "Count",
         color = NULL) +
    theme(legend.position = "bottom")

  p1
}
