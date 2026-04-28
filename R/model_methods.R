#' @keywords internal
get_posterior_moments_l.mvsusie <- function(model, l) {
  # `mu2_diag` was lifted out of `mu2_cache` to a top-level (L, J, R)
  # array; reassemble the per-l reduced-stats list with both the
  # cached (R x R / scalar) fields AND the per-l mu2_diag slice so the
  # accessor's contract is preserved.
  cache_l <- model$mu2_cache[[l]]
  if (!is.null(cache_l) && !is.null(model$mu2_diag)) {
    cache_l$mu2_diag <- model$mu2_diag[l, , , drop = TRUE]
    if (!is.matrix(cache_l$mu2_diag))
      cache_l$mu2_diag <- matrix(cache_l$mu2_diag, ncol = dim(model$mu)[3])
  }
  list(
    mu  = model$mu[l, , , drop = TRUE],    # J x R
    mu2 = cache_l                          # list of reduced stats
  )
}

#' @keywords internal
get_posterior_mean_l.mvsusie <- function(model, l) {
  drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]
}

#' @keywords internal
get_posterior_mean_sum.mvsusie <- function(model) {
  compute_posterior_mean_sum_from_model(model)
}
