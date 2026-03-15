#' @title Extract Coefficient Estimates from mvsusie Fit
#'
#' @param object An mvsusie fit, such as the output from a call to
#'   \code{\link{mvsusie}} or \code{\link{mvsusie_rss}}.
#'
#' @param \dots Additional arguments (currently unused).
#'
#' @return An (J+1) x R matrix, where J is the number of predictors
#'   and R is the number of outcomes or response variables. The first
#'   row gives the intercept estimate.
#'
#' @importFrom stats coef
#'
#' @method coef mvsusie
#'
#' @export coef.mvsusie
#'
#' @export
#'
coef.mvsusie <- function(object, ...) {
  L <- nrow(object$alpha)
  J <- ncol(object$alpha)
  csd <- object$X_column_scale_factors

  if (length(dim(object$mu)) == 3) {
    # Multivariate: mu is L x J x R
    R <- dim(object$mu)[3]
    b_sum <- matrix(0, J, R)
    for (l in seq_len(L)) {
      b_sum <- b_sum + drop(object$alpha[l, ]) * object$mu[l, , , drop = TRUE]
    }
    coefs <- b_sum / csd  # J x R (csd recycled or J x R matrix)
    result <- rbind(matrix(object$intercept, 1, R), coefs)
  } else {
    # Univariate: mu is L x J (or vector if L=1)
    alpha <- object$alpha
    mu <- object$mu
    if (!is.matrix(alpha)) alpha <- matrix(alpha, nrow = 1)
    if (!is.matrix(mu)) mu <- matrix(mu, nrow = 1)
    b_sum <- colSums(alpha * mu)  # J-vector
    coefs <- b_sum / csd
    result <- matrix(c(object$intercept, coefs), J + 1, 1)
  }
  return(as.matrix(result))
}

#' @title Predict Outcomes from mvsusie Fit.
#'
#' @param object An mvsusie fit, such as the output from a call to
#'   \code{\link{mvsusie}} or \code{\link{mvsusie_rss}}.
#'
#' @param newx A new X matrix for which to do predictions.
#'
#' @param \dots Additional arguments (currently unused).
#'
#' @return A matrix of predicted outcomes, with rows corresponding to
#'   samples (rows of X), and columns corresponding to outcomes.
#'
#' @importFrom stats coef
#' @importFrom stats predict
#'
#' @method predict mvsusie
#'
#' @export predict.mvsusie
#'
#' @export
#'
predict.mvsusie <- function(object, newx = NULL, ...) {
  if (missing(newx)) {
    return(object$fitted)
  } else {
    b <- coef(object)
    n <- nrow(newx)
    r <- ncol(b)
    intercept <- b[1, ]
    b <- b[-1, ]
    return(matrix(intercept, n, r, byrow = TRUE) + newx %*% b)
  }
}
