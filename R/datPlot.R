datPlot <- function(data){
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
