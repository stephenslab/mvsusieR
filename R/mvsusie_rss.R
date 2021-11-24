#' @rdname mvsusie
#'
#' @param Z J x R matrix of z-scores.
#' 
#' @param R J x J LD matrix.
#'
#' @export
#' 
mvsusie_rss = function (Z, R, L = 10, prior_variance = 50,
                        residual_variance = NULL, prior_weights = NULL,
                        estimate_prior_variance = TRUE,
                        estimate_prior_method = "EM",
                        check_null_threshold = 0, prior_tol = 1e-9,
                        compute_objective = TRUE,
                        precompute_covariances = FALSE,
                        s_init = NULL, coverage = 0.95, min_abs_corr = 0.5,
                        n_thread = 1, max_iter = 100, tol = 1e-3,
                        verbosity = 2, track_fit = FALSE) {
  is_numeric_prior =
    !(is.matrix(prior_variance) || inherits(prior_variance,"MashInitializer"))
  if (!is.null(dim(Z)) && ncol(Z) > 1 && is_numeric_prior)
    stop("Please specify prior variance for the multivariate z-scores")
  if (anyNA(Z)) {
    warning("NA values in z-scores are replaced with 0")
    Z[is.na(Z)] = 0
  }
  is_numeric_matrix(R,"R")
  if (is.null(dim(Z)))
    YtY = 1
  else
    YtY = diag(ncol(Z))
  return(mvsusie_suff_stat(XtX = R,XtY = Z,YtY = YtY,N = 2,L = L,
                           prior_variance = prior_variance,
                           residual_variance = residual_variance,
                           prior_weights = prior_weights,standardize = FALSE,
                           estimate_residual_variance = FALSE,
                           estimate_prior_variance = estimate_prior_variance,
                           estimate_prior_method = estimate_prior_method,
                           check_null_threshold = check_null_threshold,
                           prior_tol = prior_tol,
                           compute_objective = compute_objective,
                           precompute_covariances = precompute_covariances,
                           s_init = s_init,coverage = coverage,
                           min_abs_corr = min_abs_corr,
                           n_thread = n_thread,max_iter = max_iter,tol = tol,
                           verbosity = verbosity,track_fit = track_fit))
}
