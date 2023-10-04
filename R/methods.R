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
coef.mvsusie <- function (object, ...) {
  return(object$coef)
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
#' @return A matrix of predicted outcomes.
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
predict.mvsusie = function (object, newx = NULL, ...) {
  if (missing(newx))
    return(object$fitted)
  else
    return(object$intercept + newx %*% object$coef[-1,])
}
