zeroinf_test <- function(mod){
  if (!(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero"))){
    stop("Model class not supported")
  }

  # Zero Simulation
  ##################
  R <- 1000
  obsZeros <- sum(mod$model$y == 0)
  numZeros <- rep(NA, times=R)
  for(i in 1:length(numZeros)){
    if(class(mod) %in% c("glm_pois", "glm_pois_zero")){
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
    geom_histogram(aes(x = numZeros), color = 'black') +
    geom_vline(aes(xintercept = obsZeros,
                   color = "Observed Zero Count"),
               linetype = "dotted",
               size = 1) +
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
