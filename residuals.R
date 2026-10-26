#' Extract and Visualize EYule Model Residuals
#'
#' residuals is a method which extracts model residuals from a \code{"glm_EY"}
#' commonly returned by \code{\link{glm.EY}}. Optionally, it produces a half
#' normal plot with a simulated envelope of the residuals.
#'
#' The response residuals (\eqn{r_i=y_i - \mu_i}{r[i]=y[i] - \mu[i]}), Pearson
#' residuals (\eqn{r^P_i = r_i/\sigma_i}{r[i]^P = r[i]/\sigma[i]}) or randomized
#' quantile residuals are computed. The randomized quantile residuals are
#' obtained computing the cumulative probabilities that the fitted model being
#' less than \emph{y} and less or equal than \emph{y}. A random value from a
#' uniform distribution between both probabilities is generated and the value of
#' the residual is the standard normal variate with the same cumulative
#' probability. Four replications of the quantile residuals are recommended
#' because of the random component (see Dunn and Smyth, 1996 for more details).
#'
#' The functions \code{\link{plot.glm_EY}} generate a residuals against fitted
#' values plot and a Normal Q-Q plot.
#'
#' The Normal Q-Q plot may show an unsatisfactory pattern of the Pearson
#' residuals of a fitted model: then we are led to think that the model is
#' incorrectly specified.

#'The half normal plot with simulated envelope indicates that, under the
#'distribution of the response variable, the model is fine when only a few
#'points fall off the envelope.
#'
#'@param object an object of class \code{"glm_EY"}, typically the result of a
#'call to \code{\link{glm.EY}}.
#'@param type the type of residuals which should be returned. The alternatives
#'  are: "pearson" (default), "response" and "quantile". Can be abbreviated.
#'@param envelope a logical value indicating whether the envelope should be
#'  computed.
#'@param rep	number of replications for envelope construction. Default is 19,
#'  that is the smallest 95 percent band that can be built.
#'@param title	a string indicating the main title of the envelope.
#'@param ... further arguments passed to or from other methods.
#'
#'@return Residual values.
#'
#'@seealso \code{\link{plots}}
#'
#'@references Peter K. Dunn and Gordon K. Smyth (1996). "Randomized quantile
#'  residuals". Journal of Computational and Graphical Statistics, 5(3), pp.
#'  236-244.
#'
#'  A. C. Atkinson (1981). "Two graphical displays for outlying and influential
#'  observations in regression". Biometrika, 68(1), pp. 13–20.
#'
#'@name residuals
NULL

#' @rdname residuals
#' @examples
#' ## Estimate a EYule model
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.EY(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.gamma = numbids ~ 1, data = Bids)
#'
#' ## Compute residuals
#'
#' r <- residuals(fit)
#' @export
residuals.glm_EY <- function(object, type = c("pearson", "response", "quantile"),
                            envelope = FALSE, rep = 19,
                            title = "Simulated Envelope of Residuals", ...) {
  stopifnot(is.logical(envelope))
  stopifnot(is.numeric(rep), length(rep) == 1, rep >= 1)
  stopifnot(is.character(title), length(title) == 1)

  type <- match.arg(type)

  if (! envelope)
    return(res_EY(object, type))

  residuals_sim <- matrix(0, nrow = length(object$residuals), ncol = rep)
  pb <- utils::txtProgressBar(min = 1, max = rep, initial = 1, style = 3)
  for (x in 1:rep) {
    utils::setTxtProgressBar(pb, x)
    tmp <- simulation_EY(object, type)
    residuals_sim[, x] <- tmp
  }

  residuals_sim <- apply(residuals_sim, 2, sort)
  minima <- apply(residuals_sim, 1, min)
  maxima <- apply(residuals_sim, 1, max)

  resi <- sort(res_EY(object, type))
  n <- length(resi)
  t <- 1:n

  normal_score <- stats::qnorm(t / (n + 1))

  xx <- c(normal_score, rev(normal_score))
  yy <- c(minima, rev(maxima))

  type_r <- paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)))
  graphics::plot(normal_score,
                 resi,
                 type = "l",
                 xlab = "Standard normal quantiles",
                 ylab = paste(type_r, "residuals"),
                 main = title)
  graphics::polygon(xx, yy, col = "gray", border = NA)
  graphics::lines(normal_score, resi)

  structure(
    list(
      type = type,
      residuals = resi,
      sim_residuals = residuals_sim
    ),
    class = "residuals.glm_EY"
  )
}

res_EY <- function(object, type) {
  lambda <- object$lambdas
  gamma  <- object$gammas
  if (type == "pearson") {
    # mu       <- object$fitted.values
    # variance <- lambda + (lambda - (gamma - 1)) * mu - mu ^ 2
    variance <- variances_EY(lambda, gamma, object$maxiter_series, object$tol)
    r        <- as.vector(object$residuals / sqrt(variance))
    names(r) <- seq(r)
    return(r)
  }
  if (type == "response")
    return(object$residuals)

  # type == "quantile"
  y <- object$y

  #   if (! is.null(object$y)) {
  #   y <- object$y
  # } else if (! is.null(object$model.mu)) {
  #   y <- stats::model.response(object$model.mu)
  # } else {
  #   formula.mu <- stats::as.formula(object$call[["formula.mu"]])
  #   a.mu <- stats::model.frame(formula.mu, data = object$data)
  #   y <- stats::model.response(a.mu)
  # }

  a <- ifelse(y == 0, 0, pEY(y - 1, gamma, lambda))
  b <- pEY(y, gamma, lambda)
  u <- stats::runif(length(y), a, b)
  qr <- stats::qnorm(u)
  return(qr)
}

simulation_EY <- function(object, type) {
  repeat {
    lambdas <- object$lambdas
    gammas  <- object$gammas
    response_n <- get_response(object$formula.mu)
    data <- object$data
    data[response_n] <- sapply(seq_len(nrow(data)),
                               function(i) simulate_EY(gammas[i], lambdas[i]))
    the_call <- object$call
    the_call[["data"]] = data
    fit <- tryCatch(eval(the_call),
                    error = function(e) TRUE,
                    warning = function(e) TRUE)
    if (is.logical(fit))
      next
    if (any(is.nan(res_EY(fit, type))))
      next
    if (fit$code < 0)
      next
    break
  }
  return(res_EY(fit, type))
}

simulate_EY <- function(gamma, lambda) {
#  pochammer <- function(a, r) if (r == 0) 1 else prod(a:(a + r - 1))
  u <- stats::runif(1)
  y <- 1
  value <- f21(lambda,gamma)
  fac <- 1/value
  p <- 1/value
  while (p < u) {
    fac <- fac * y * lambda / (gamma + y - 1)
    p <- p + fac
    y <- y + 1
  }
  y - 1
}

get_response <- function(formula) {
  tt <- stats::terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1]
  response <- attr(tt, "response")
  vars[response]
}
