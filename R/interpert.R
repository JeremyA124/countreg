interpert <- function(mod){

  if (!(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero", "glm_pois_GP2"))){
    stop("Model class not supported")
  }

  cat("Call:\n")
  print(mod$call)
  cat("\n")
  cat("Residuals:\n")
  print(summary(unlist(mod$residuals)))
  cat("\n")

  if(class(mod) %in% c("glm_pois", "glm_negb", "glm_pois_GP2")){
    coef.frame <- data.frame("Estimate"=(exp(mod$coefficients$betas)-1)*100,
                             "Std.error"=mod$coefficients$std.error,
                             "Z"=mod$summary$Z,
                             "P-value"=mod$summary$p.val)
    coef.frame$Estimate[1] = mod$coefficients$betas[1]

    cat("Count Coefficients:\n")
    print(coef.frame)
    cat("\n")
  }

  if(class(mod) %in% c("glm_pois_zero", "glm_negb_zero")){
    coef.frame <- data.frame("Estimate"=(exp(mod$coefficients$count)-1)*100,
                             "Std.error"=mod$coefficients$std.error.count,
                             "Z"=mod$summary$count$Z,
                             "P-value"=mod$summary$count$p.val)
    coef.frame$Estimate[1] = mod$coefficients$count[1]

    coef.frame2 <- data.frame("Estimate"=(exp(mod$coefficients$zero)-1)*100,
                             "Std.error"=mod$coefficients$std.error.zero,
                             "Z"=mod$summary$zero$Z,
                             "P-value"=mod$summary$zero$p.val)
    coef.frame2$Estimate[1] = mod$coefficients$zero[1]

    cat("Count Coefficients:\n")
    print(coef.frame)
    cat("\n")
    cat("Zero Coefficients:\n")
    print(coef.frame2)
  }
}
