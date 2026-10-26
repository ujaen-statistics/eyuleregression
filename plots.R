#' Plot Diagnostics for glm_EY Objects
#'
#' Two plots are currently available: a plot of residuals against fitted values
#' and a Normal Q-Q plot.
#'
#' @param x \code{glm_EY} object, typically the result of \code{\link{glm.EY}}.
#' @param type the type of residuals. The alternatives are: "quantile"
#'   (default), "pearson" and "response". Can be abbreviated.
#' @param ask logical; if TRUE, the user is asked before each plot, see
#'   \code{\link[graphics]{par}}(ask=.).
#' @param ... other parameters to be passed through to plotting functions.
#' @name plots
NULL

#' @rdname plots
#' @examples
#' ## Fit the EYule model
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.EY(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.gamma = numbids ~ 1, data = Bids)
#' oldpar <- par(mfrow = c(1, 2))
#'
#' ## Plot diagnostics
#'
#' plot(fit)
#' par(oldpar)
#' @export
plot.glm_EY <- function(x, type = c("quantile", "pearson", "response"),
                       ask = prod(graphics::par("mfcol")) < 2 &&
                             grDevices::dev.interactive(),
                       ...) {
  type <- match.arg(type)
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  grDevices::dev.hold()
  graphics::plot(stats::fitted.values(x),
                 stats::residuals(x, type = type),
                 xlab = "Fitted values",
                 ylab = "Residuals",
                 main = "Residuals vs Fitted"
  )
  grDevices::dev.flush()
  stats::qqnorm(stats::residuals(x, type = type))
  stats::qqline(stats::residuals(x, type = type))
  invisible()
}
