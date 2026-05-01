#' Plotting a Scatter plot Matrix of Data
#'
#' \code{dat_plot} is used to summarize/model data through the use of a scatter plot
#' matrix of a data set's variables (i.e. columns).
#'
#' @param data A dataframe object
#'
#' @details
#' \code{dat_plot} require the installation of \code{\link{GGally}} and \code{\link{ggplot2}}. The
#' matrix is formated to plot/compute correlations and data summarizations based on data
#' structure.
#'
#' NOTE: It's not recommended to use this function past 15 variables, due to complexity.
#' Try sub-setting data set or utilize \code{\link{ggplot2}}.
#'
#' @author
#' Implementation of \code{dat_plot} was authored by Jeremy Artiga, with aid
#' from William Cipolli at Colgate University.
#'
#' @references
#' Long, J. S. (1997).
#' Regression Models for Categorical and Limited Dependent Variables.
#' Sage Publications.
#'
#' @examples
#' ## Example using bioChemists dataset from pscl (Campbell & Mahon, 1974)
#' utils::data(bioChemists, package = "pscl")
#' dat_plot(bioChemists)
#'
#' ## Example using mtcars data set (Motor Trend, 1974)
#' utils::data(mtcars)
#' dat_plot(cars)
#'
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 theme_bw theme element_text scale_fill_grey scale_color_grey
#' @export

dat_plot <- function(data){
  require(GGally)
  require(ggplot2)

  upper <- list(continuous = "points", discrete = wrap("colbar", size = 0), combo = "box_no_facet", na = "na")
  lower <- list(continuous = "cor", discrete = wrap("colbar", size = 0), combo = "box_no_facet", na = "na")
  ggpairs(data, progress = F, upper = upper, lower = lower) +
    theme_bw()+
    theme(axis.text.x = element_text(angle=60, vjust = 1, hjust=1)) +
    scale_fill_grey()+
    scale_color_grey()
}
