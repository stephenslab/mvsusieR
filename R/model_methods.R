#' @keywords internal
get_prior_variance_l.mvsusie <- function(model, l) {
  # Return the scalar scale. Methods that need the full RxR matrix
  # reconstruct it as V[l] * model$V_structure.
  model$V[l]
}

#' @keywords internal
set_prior_variance_l.mvsusie <- function(model, l, V) {
  model$V[l] <- V
  model
}

#' @keywords internal
get_alpha_l.mvsusie <- function(model, l) {
  model$alpha[l, ]
}

#' @keywords internal
get_posterior_moments_l.mvsusie <- function(model, l) {
  list(
    mu  = model$mu[l, , , drop = TRUE],    # J x R
    mu2 = model$mu2_cache[[l]]             # list of reduced stats
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
