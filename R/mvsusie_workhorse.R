# =============================================================================
# MULTIVARIATE SUSIE VIA SUSIER WORKHORSE
#
# Constructs data and params objects, then calls susieR's susie_workhorse
# for the IBSS loop. This replaces mvsusie_core for the S3 path.
# =============================================================================

#' Run multivariate SuSiE using susieR's IBSS framework
#'
#' @param data An S3 data object of class \code{mv_individual} or \code{mv_ss}.
#' @param L Number of single effects.
#' @param prior_variance Prior variance: R x R matrix, or scalar for univariate.
#' @param prior_weights Prior inclusion probability weights (length J).
#' @param estimate_residual_variance Logical.
#' @param estimate_prior_variance Logical.
#' @param estimate_prior_method Character: "EM" or "optim".
#' @param check_null_threshold Numeric threshold for null check.
#' @param compute_objective Logical.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param prior_tol Tolerance for trimming null effects.
#' @param track_fit Logical.
#' @param verbosity Verbosity level.
#' @param coverage CS coverage.
#' @param min_abs_corr Minimum absolute correlation for CS.
#'
#' @return A fitted mvsusie model.
#'
#' @keywords internal
mvsusie_workhorse <- function(data, L, prior_variance,
                               prior_weights = NULL,
                               estimate_residual_variance = FALSE,
                               estimate_prior_variance = TRUE,
                               estimate_prior_method = "optim",
                               check_null_threshold = 0,
                               compute_objective = TRUE,
                               max_iter = 100,
                               tol = 1e-3,
                               prior_tol = 1e-9,
                               track_fit = FALSE,
                               verbosity = 2,
                               coverage = 0.95,
                               min_abs_corr = 0.5,
                               n_thread = 1,
                               model_init = NULL) {

  J <- data$p
  R <- data$R

  # Default prior weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / J, J)
  } else {
    prior_weights <- prior_weights / sum(prior_weights)
  }

  # Determine estimation method string for susieR
  if (!estimate_prior_variance) {
    est_method <- "none"
  } else {
    est_method <- match.arg(estimate_prior_method, c("EM", "optim"))
  }

  # Determine prior type
  is_mash_prior <- class(prior_variance)[1] == "mash_prior"
  if (is.matrix(prior_variance) || is_mash_prior) {
    mv_prior_type <- "multivariate"
  } else {
    mv_prior_type <- "simple"
  }

  # Build params object matching susieR's expectations
  params <- list(
    L                        = L,
    scaled_prior_variance    = prior_variance,
    residual_variance        = data$residual_variance,
    prior_weights            = prior_weights,
    null_weight              = 0,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance  = estimate_prior_variance,
    estimate_prior_method    = est_method,
    check_null_threshold     = check_null_threshold,
    prior_tol                = prior_tol,
    max_iter                 = max_iter,
    tol                      = tol,
    verbose                  = (verbosity > 1),
    track_fit                = track_fit,
    compute_objective        = compute_objective,
    coverage                 = coverage,
    min_abs_corr             = min_abs_corr,
    residual_variance_lowerbound  = 0,
    residual_variance_upperbound  = Inf,
    refine                   = FALSE,
    model_init               = model_init,
    convergence_method       = "elbo",
    use_servin_stephens      = FALSE,
    unmappable_effects       = "none",
    n_purity                 = 100,
    intercept                = TRUE,
    standardize              = TRUE,
    compute_univariate_zscore = FALSE,
    # Multivariate-specific
    mv_prior_type            = mv_prior_type,
    n_thread                 = n_thread
  )

  # Call susieR's workhorse
  s <- susie_workhorse(data, params)

  return(s)
}
