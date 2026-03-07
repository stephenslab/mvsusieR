#' @rdname mvsusie
#'
#' @param XtX A J x J matrix \eqn{X^TX} in which the columns of \eqn{X}
#'   are centered to have mean zero.
#'
#' @param XtY A J x R matrix \eqn{X^TY} in which the columns of
#'   \eqn{X} and \eqn{Y} are centered to have mean zero.
#'
#' @param YtY An R x R matrix \eqn{Y^TY} in which the columns of
#' \eqn{Y} are centered to have mean zero.
#'
#' @param N The sample size.
#'
#' @param X_colmeans A vector of length J giving the column means of
#'   \eqn{X}. If it is provided with \code{Y_colmeans}, the intercept is
#'   estimated; otherwise, the intercept is \code{NA}.
#'
#' @param Y_colmeans A vector of length R giving the column means of
#' \eqn{Y}. If it is provided with \code{X_colmeans}, the intercept is
#'   estimated; otherwise, the intercept is \code{NA}.
#'
#' @export
#'
mvsusie_suff_stat <- function(XtX, XtY, YtY, N, L = 10, X_colmeans = NULL,
                              Y_colmeans = NULL, prior_variance = 0.2,
                              residual_variance = NULL, prior_weights = NULL,
                              standardize = TRUE,
                              estimate_residual_variance = FALSE,
                              estimate_prior_variance = TRUE,
                              estimate_prior_method = "optim",
                              check_null_threshold = 0, prior_tol = 1e-9,
                              compute_objective = TRUE,
                              precompute_covariances = FALSE, s_init = NULL,
                              coverage = 0.95, min_abs_corr = 0.5, n_thread = 1,
                              max_iter = 100, tol = 1e-3, verbosity = 2,
                              track_fit = FALSE) {
  # For R=1 with scalar prior, convert from susieR "scaled prior variance"
  # convention to absolute prior variance: actual_V = scaled_V * sigma2.
  XtY <- as.matrix(XtY)
  YtY <- as.matrix(YtY)
  R <- ncol(XtY)
  is_numeric_prior <- is.numeric(prior_variance) && !is.matrix(prior_variance)
  if (R == 1L && is_numeric_prior) {
    rv <- if (is.null(residual_variance)) as.numeric(YtY / (N - 1))
          else as.numeric(residual_variance)
    prior_variance <- prior_variance * rv
  }
  mvsusie_suff_stat_s3(XtX, XtY, YtY, N, L = L,
                       X_colmeans = X_colmeans,
                       Y_colmeans = Y_colmeans,
                       prior_variance = prior_variance,
                       residual_variance = residual_variance,
                       prior_weights = prior_weights,
                       standardize = standardize,
                       estimate_residual_variance = estimate_residual_variance,
                       estimate_prior_variance = estimate_prior_variance,
                       estimate_prior_method = estimate_prior_method,
                       check_null_threshold = check_null_threshold,
                       prior_tol = prior_tol,
                       compute_objective = compute_objective,
                       model_init = s_init,
                       coverage = coverage, min_abs_corr = min_abs_corr,
                       n_thread = n_thread,
                       max_iter = max_iter, tol = tol,
                       verbosity = verbosity, track_fit = track_fit)
}
