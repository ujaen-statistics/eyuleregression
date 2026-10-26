#'Expected Probabilities and Frequencies for the EYule Model
#'
#'The \code{EY_expected} function calculate the probability distribution of the
#'count response variable Y for each observation and obtain the corresponding
#'expected frequencies. It is an informal way of assessing the fit of the EY
#'model by comparing the predicted distribution of counts with the observed
#'distribution.
#'
#'The average expected probabilities are computed as \deqn{\bar(Pr)(y=k) =
#'\frac{1}{n} \sum_{i=1}^n \widehat{Pr}(y_i = k | x_i)}
#'
#'The expected frequencies are obtained by multiplying by n.
#'
#'Two measures are offered for summarizing the comparison between expected and
#'observed frequencies: the sum of the absolute value of differences and the sum
#'of the square of differences (similar to the Pearson statistic of goodness of
#'fit).
#'
#'@param object a fitted object of class inheriting from \code{"glm_EY"}.
#'
#'@return A list containing the following components:
#'
#'  \item{\code{frequencies}}{the expected counts for the EY fit.}
#'  \item{\code{observed_freq}}{the observed distribution.}
#'  \item{\code{probabilities}}{the expected distribution for the EY fit.}
#'  \item{\code{dif}}{sum of the absolute value of differences between
#'  \code{frequencies} and \code{observed_freq}.} \item{\code{chi2}}{sum of the
#'  square of differences between \code{frequencies} and \code{observed_freq}.}
#'
#'@references
#'
#'J. M. Hilbe (2011). Negative Binomial Regression. (2nd ed.). Cambridge
#'University Press.
#'
#'M. Scott Long and Jeremy Freese (2014). Regression Models for Categorical
#'Dependent Variables using STATA. (3rd ed.). Stata Press.
#'
#'@name expected
NULL

#' @rdname expected
#' @examples
#' ## Fit a EYule model
#'
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.EY(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                  whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                  formula.gamma = numbids ~ 1, data = Bids)
#'
#' ## Compute the expected probabilities and the frequencies
#'
#' EY_expected(fit)
#' @export
EY_expected <- function(object) {
  stopifnot(inherits(object, "glm_EY"))

  gamma  <- object$gammas
  lambda <- object$lambdas

  # if (! is.null(object$y)) {
  #   y <- object$y
  # } else if (! is.null(object$model.mu)) {
  #   y <- stats::model.response(object$model.mu)
  # } else {
  #   formula.mu <- stats::as.formula(object$call[["formula.mu"]])
  #   a.mu <- stats::model.frame(formula.mu, data = object$data)
  #   y <- stats::model.response(a.mu)
  # }

  y <- object$y

  row <- max(y) + 1
  col <- length(lambda)
  i <- 0:(row - 1)
  f0 <- matrix(0, col)
  prob <- matrix(0, nrow = row + 1, ncol = col)
  dist <- matrix(0, nrow = row + 1, ncol = col)
  for (j in seq_len(col)) {
    prob[1:row, j] <- dEY(i, gamma[j], lambda[j])
    prob[row + 1, j] <- pEY(row - 1, gamma[j], lambda[j], lower.tail = FALSE)
    dist[, j] <- object$weights[j] * prob[, j]
  }

  estimated_freq <- rowSums(dist)
  estimated_prob <- estimated_freq / sum(estimated_freq)

  observed_freq <- table(factor(y, levels = 0:row))
  observed_prob <- observed_freq / sum(observed_freq)

  #  dif <- sum(abs(observed_prob - estimated_prob))
  dif <- sum(abs(observed_freq - estimated_freq))

  chi2 <- sum((observed_freq[1:(row-1)] - estimated_freq[1:(row-1)]) ^ 2 / estimated_freq[1:(row-1)]) +
    (sum(observed_freq[row]) - sum(estimated_freq[row:(row+1)])) ^ 2 / sum(estimated_freq[row:(row+1)])

  #  barplot(observed_freq[0:row], xlab = "Y", ylab = "Frequencies")
  plot(observed_freq, type="b",  xlab = "Y", ylab = "Frequencies", main = "0bserved vs Predicted Frequencies",
       ylim = c(0,max(observed_freq,estimated_freq)+1))
  lines(0:row, estimated_freq, type = "b", pch = 2, col = "blue")
  #  legend("topright", legend=c(paste("Dif = ", round(dif, 3)), paste(expression(chi ^ 2), "=", round(chi2, 3))))
  legend("topright", legend=c("Observed","Predicted"), lty = c(1,1), pch = c(1,2), col = c(1,4))

  list(
    frequencies = estimated_freq,
    observed_freq = observed_freq,
    probabilities = estimated_prob,
    dif = dif,
    chi2 = chi2
  )
}
