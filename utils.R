stopping <- function (x, tol) {
  all(abs(x) <= tol, na.rm = TRUE)
}
# stopping <- function (x, tol) all(abs(x) <= tol) & all(!is.nan(x))

# In the following functions a tolerance value less than 0 can be used to
# avoid the stopping criterion based on tolerance

# 2f1(1,1,gamma;lambda) to normalize EYule
# f21 <- function(lambda, gamma, maxiter_series = 2000, tol = 0) {
#   library(hypergeo)
#   series <- hypergeo(1, 1, gamma, lambda, tol, maxiter_series)
#   return(Re(series))
# }

f21 <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10) {
#  if (lambda > 1) stop("lambda no puede ser mayor que 1")
#  else if (lambda == 1){
#    if (gamma > 2){
#      series <- (gamma(gamma-2)/gamma(gamma-1))*(gamma(gamma)/gamma(gamma-1))
#      return(Re(series))}
#    else stop("gamma debe ser mayor que 2")}
#  else
  fac  <- 1
  temp <- 1
  L    <- gamma
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * n * lambda / L
    series <- temp + fac
    if (stopping(series - temp, tol)){
      return(Re(series))
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series))
}

# E[Y]
# means_EY <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10, f21_comp = NULL) {
#   if (is.null(f21_comp))
#     f21_comp <- f21(lambda, gamma, maxiter_series, tol)
#   mu <- (lambda - (gamma-1)*(1-1/f21_comp))/(1-lambda)
#   return(mu)
# }

means_EY <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10, f21_comp = NULL) {
  if (is.null(f21_comp))
    f21_comp <- f21(lambda, gamma, maxiter_series, tol)
  L    <- gamma
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * n * lambda / L
    series <- temp + n * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f21_comp)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series) / f21_comp)
}

# Var(Y)
# variances_EY <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
#                          f21_comp = NULL, means_EY_comp = NULL)
# {
#   if (is.null(f21_comp))
#     f21_comp <- f21(lambda, gamma, maxiter_series, tol)
#   if (is.null(means_EY_comp))
#     means_EY_comp <- means_EY(lambda, gamma, maxiter_series, tol, f21_comp)
#   s2 <- (lambda + (2*lambda - gamma + 1)*means_EY_comp)/(1-lambda) - means_EY_comp^2
#   return(s2)
# }

variances_EY <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                         f21_comp = NULL, means_EY_comp = NULL)
{
  if (is.null(f21_comp))
    f21_comp <- f21(lambda, gamma, maxiter_series, tol)
  if (is.null(means_EY_comp))
    means_EY_comp <- means_EY(lambda, gamma, maxiter_series, tol, f21_comp)
  L    <- gamma
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * n * lambda / L
    series <- temp + n ^ 2 * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f21_comp - means_EY_comp ^ 2)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series) / f21_comp - means_EY_comp ^ 2)
}

# E[psi(gamma + Y)]
means_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                       f21_comp = NULL) {
  if (is.null(f21_comp))
    f21_comp <- f21(lambda, gamma, maxiter_series, tol)
  L    <- gamma
  fac  <- 1
  temp <- digamma(gamma)
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * n * lambda / L
    series <- temp + digamma(gamma + n) * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f21_comp)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series) / f21_comp)
}

# Cov(Y, psi(gamma + Y))
covars_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                        f21_comp = NULL, means_EY_comp = NULL) {
  if (is.null(f21_comp))
    f21_comp <- f21(lambda, gamma, maxiter_series, tol)
  if (is.null(means_EY_comp))
    means_EY_comp <- means_EY(lambda, gamma, maxiter_series, tol, f21_comp)
  L    <- gamma
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * n * lambda / L
    series <- temp + n * digamma(gamma + n) * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f21_comp -
               means_EY_comp *
               means_psiy(lambda, gamma, maxiter_series, tol, f21_comp))
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series) / f21_comp -
           means_EY_comp *
           means_psiy(lambda, gamma, maxiter_series, tol, f21_comp))
}

# Var(psi(gamma + Y))
variances_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10)
{
  f21_comp <- f21(lambda, gamma, maxiter_series, tol)
  L    <- gamma
  fac  <- 1
  temp <- digamma(gamma) ^ 2
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * n * lambda / L
    series <- temp + digamma(gamma + n) ^ 2 * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f21_comp - means_psiy(lambda, gamma, maxiter_series, tol, f21_comp) ^ 2)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  Re(series) / f21_comp - means_psiy(lambda, gamma, maxiter_series, tol, f21_comp) ^ 2
}
