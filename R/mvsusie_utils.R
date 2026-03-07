# Weighted log-sum-exp (softmax).
#
# Given log-values and weights, compute the normalized weights
# (softmax) and the log of their weighted sum.
#
# @param value A numeric vector of log-values (or values if log = FALSE).
# @param weight A numeric vector of weights (same length as value).
# @param log If TRUE (default), value is already on log scale.
#
# @return A list with \code{weights} (normalized) and \code{log_sum}.
#
# @keywords internal
compute_softmax <- function(value, weight, log = TRUE) {
  if (length(value) != length(weight))
    stop("Values and their weights should have equal length")
  if (!log) value <- log(value)
  mvalue <- max(value)
  w <- exp(value - mvalue)
  w_weighted <- w * weight
  weighted_sum_w <- sum(w_weighted)
  list(weights = as.vector(w_weighted / weighted_sum_w),
       log_sum = log(weighted_sum_w) + mvalue)
}

#' @title Local false sign rate (lfsr) for single effects
#'
#' @details This function returns the lfsr for identifying nonzero
#'   single effects, separately for each condition.
#'
#' @param alpha L x J matrix.
#'
#' @param clfsr L x J x R conditonal lfsr.
#'
#' @return L x R matrix of lfsr
#'
#' @export
#'
mvsusie_single_effect_lfsr <- function(clfsr, alpha) {
  if (!is.array(clfsr) && is.na(clfsr)) {
    return(as.numeric(NA))
  } else {
    return(do.call(
      cbind,
      lapply(
        1:dim(clfsr)[3],
        function(r) {
          clfsrr <- clfsr[, , r]
          if (is.null(nrow(clfsrr))) {
            clfsrr <- matrix(clfsrr, 1, length(clfsrr))
          }
          return(pmax(0, rowSums(alpha * clfsrr)))
        }
      )
    ))
  }
}

#' @title Local false sign rate (lfsr) for variables.
#'
#' @details This function returns the lfsr for identifying nonzero
#'   effects for each condition.
#'
#' @param alpha L x J matrix.
#'
#' @param clfsr L x J x R conditonal lfsr.
#'
#' @param weighted Set \code{weighted = TRUE} to weight lfsr by PIP;
#'   otherwise set \code{weighted = FALSE}.
#'
#' @return J x R lfsr matrix.
#'
#' @export
#'
mvsusie_get_lfsr <- function(clfsr, alpha, weighted = TRUE) {
  if (!is.array(clfsr) && is.na(clfsr)) {
    return(as.numeric(NA))
  } else {
    if (weighted) {
      alpha <- alpha
    } else {
      alpha <- matrix(1, nrow(alpha), ncol(alpha))
    }
    return(do.call(
      cbind,
      lapply(
        1:dim(clfsr)[3],
        function(r) {
          true_sign_mat <- alpha * (1 - clfsr[, , r])
          pmax(1e-20, 1 - apply(true_sign_mat, 2, max))
        }
      )
    ))
  }
}

