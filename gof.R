#Goodness of fit measures

MPB <- function(fit){
  return(sum(fit$weights * (fit$y - fit$fitted.values)) / sum(fit$weights))
}

MAD <- function(fit){
  return(sum(fit$weights * abs(fit$y - fit$fitted.values)) / sum(fit$weights))
}

MSPE <- function(fit){
  return(sum(fit$weights * (fit$y - fit$fitted.values)^2) / sum(fit$weights))
}

pearson<- function(fit){
  return(sum(res_EY(fit,"pearson")^2))
}