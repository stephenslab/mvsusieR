#' Run multivariate SuSiE using susieR's IBSS framework
#'
#' @param data An S3 data object of class \code{mv_individual} or \code{mv_ss}.
#' @param L Number of single effects.
#' @param prior_variance Prior variance: R x R matrix, or scalar for univariate.
#' @param prior_weights Prior inclusion probability weights (length J).
#' @param estimate_residual_variance Logical.
#' @param estimate_prior_variance Logical.
#' @param estimate_prior_method Character: "EM", "optim", or "uniroot".
#' @param check_null_threshold Numeric threshold for null check.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param prior_tol Tolerance for trimming null effects.
#' @param track_fit Logical.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#' @param coverage CS coverage.
#' @param min_abs_corr Minimum absolute correlation for CS.
#'
#' @return A fitted mvsusie model.
#'
#' @keywords internal
mvsusie_workhorse <- function(data, L, prior_variance,
                               prior_weights = NULL,
                               estimate_residual_variance = TRUE,
                               estimate_prior_variance = TRUE,
                               estimate_prior_method = "optim",
                               estimate_prior_mixture_weights = TRUE,
                               mixture_weight_method = "mixsqp",
                               check_null_threshold = 0,
                               convergence_method = "elbo",
                               max_iter = 100,
                               tol = 1e-3,
                               prior_tol = 1e-9,
                               track_fit = FALSE,
                               verbose = TRUE,
                               coverage = 0.95,
                               min_abs_corr = 0.5,
                               precompute_covariances = FALSE,
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
    est_method <- match.arg(estimate_prior_method, c("EM", "optim", "uniroot"))
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
    verbose                  = isTRUE(verbose),
    track_fit                = track_fit,
    coverage                 = coverage,
    min_abs_corr             = min_abs_corr,
    residual_variance_lowerbound  = 0,
    residual_variance_upperbound  = Inf,
    refine                   = FALSE,
    model_init               = model_init,
    unmappable_effects       = "none",
    n_purity                 = 100,
    use_servin_stephens      = FALSE,  # required by susieR::ibss_finalize
    # Multivariate-specific
    precompute_eigendecomp   = precompute_covariances,  # from caller's precompute_cache
    n_thread                 = n_thread,
    estimate_prior_mixture_weights = estimate_prior_mixture_weights,
    mixture_weight_method    = mixture_weight_method,
    convergence_method       = convergence_method
  )

  # Call susieR's workhorse
  s <- susie_workhorse(data, params)

  return(s)
}

#' Multivariate SUm of Single Effect (SuSiE) Regression
#'
#' Performs a Bayesian multiple linear regression of Y on X.
#'   That is, this function fits the regression model \deqn{Y = \sum_l X
#'   b_l + e,} where the elements of \eqn{e} are \emph{i.i.d.} normal
#'   with zero mean and variance \code{residual_variance}, and the sum
#'   \eqn{\sum_l b_l} is a vector of p effects to be estimated. The
#'   SuSiE assumption is that each \eqn{b_l} has exactly one non-zero
#'   element.
#'
#' @rdname mvsusie
#'
#' @param X N by J matrix of covariates.
#'
#' @param Y Vector of length N, or N by R matrix of response
#'   variables.
#'
#' @param L Maximum number of non-zero effects.
#'
#' @param prior_variance Can be either (1) a vector of length L, or a
#'   scalar, for scaled prior variance when Y is univariate (which
#'   should then be equivalent to \code{\link[susieR]{susie}}); or (2) a
#'   matrix for a simple multivariate regression; or (3) a mixture prior
#'   from \code{\link{create_mixture_prior}}.
#'
#' @param residual_variance The residual variance (defaults to the
#'   sample variance of Y).
#'
#' @param prior_weights A vector of length p giving the prior
#'   probability that each element is non-zero. Note that the prior
#'   weights need to be non-negative but do not need to sum to 1; they
#'   will automatically be normalized to sum to 1 so that they represent
#'   probabilities. The default setting is that the prior weights are
#'   the same for all variables.
#'
#' @param standardize Logical flag specifying whether to standardize
#'   columns of X to unit variance prior to fitting. If you do not
#'   standardize you may need to think more carefully about specifying
#'   the scale of the prior variance. Whatever the value of standardize,
#'   the coefficients (returned by \code{\link{coef}}) are for X
#'   on the original input scale. Note that any column of X with zero
#'   variance is not standardized, but left as is.
#'
#' @param intercept Should intercept be fitted or set to zero. Setting
#'   \code{intercept = FALSE} is generally not recommended.
#'
#' @param estimate_residual_variance When
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated at each iteration using \eqn{E_q[R'R] / n}; otherwise it
#'   is fixed. For multivariate Y the estimate is a full \eqn{r \times r}
#'   covariance matrix. Supported for all missing data methods: the
#'   update formula uses expected sufficient statistics (not the ELBO
#'   value), and the impute method includes a \code{Y_cov} correction
#'   for imputation uncertainty. Defaults to \code{TRUE} for
#'   \code{mvsusie()}, and \code{FALSE} for \code{mvsusie_ss()} and
#'   \code{mvsusie_rss()}.
#'
#' @param estimate_prior_variance When \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated; otherwise it is
#'   fixed. Currently \code{estimate_prior_variance = TRUE} only works
#'   for univariate Y, or for multivariate Y when the prior variance is
#'   a matrix.
#'
#' @param estimate_prior_method The method used for estimating the
#'   prior variance; valid choices are \code{"optim"}, \code{"EM"},
#'   or \code{"uniroot"}.
#'
#' @param estimate_prior_mixture_weights When \code{TRUE} and
#'   \code{prior_variance} is a mixture prior, the mixture weights
#'   are updated at each iteration. Components with near-zero weight
#'   are pruned.
#'
#' @param mixture_weight_method Method for updating mixture weights;
#'   \code{"mixsqp"} or \code{"EM"}.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   the estimate is compared against the null, and the prior variance
#'   is set to zero unless the log-likelihood using the estimate is
#'   larger than that of null by this threshold. For example, setting
#'   \code{check_null_threshold = 0.1} will \dQuote{nudge} the estimate
#'   towards zero. When used with \code{estimate_prior_method = "EM"},
#'   setting \code{check_null_threshold = NA} will skip this check.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to this value at the end of the analysis and
#'   exclude a single effect from PIP computation if the estimated prior
#'   variance is smaller than it.
#'
#' @param precompute_cache If \code{precompute_cache =
#'   TRUE}, precomputes eigendecomposition and covariance quantities
#'   to speed up computations at the cost of increased memory usage.
#'
#' @param model_init A previous model fit with which to initialize.
#'
#' @param coverage Coverage of credible sets.
#'
#' @param min_abs_corr Minimum of absolute value of correlation
#'   allowed in a credible set. The setting \code{min_abs_corr = 0.5}
#'   corresponds to squared correlation of 0.25, which is a commonly
#'   used threshold for genotype data in genetics studies.
#'
#' @param missing_y_method Method for handling missing values in Y;
#'   \code{"approximate"} or \code{"exact"}.
#'
#' @param compute_univariate_zscore When
#'   \code{compute_univariate_zscore = TRUE}, the z-scores from the
#'   per-variable univariate regressions are returned. (Note that these
#'   z-scores are not actually used to fit the multivariate susie
#'   model.)
#'
#' @param n_thread Maximum number of threads to use for parallel
#'   computation (only applicable when a mixture prior is used).
#'
#' @param max_iter Maximum number of iterations to perform.
#'
#' @param tol The model fitting will terminate when the increase in
#'   ELBOs between two successive iterations is less than \code{tol}.
#'
#' @param verbose If \code{TRUE}, print progress messages during model
#'   fitting. Default is \code{TRUE}.
#'
#' @param track_fit Add attribute \code{trace} to the return value
#'   which records the algorithm's progress at each iteration.
#'
#' @return A multivariate susie fit, which is a list with some or all
#' of the following elements:
#'
#' \item{alpha}{L by p matrix of posterior inclusion probabilities.}
#'
#' \item{mu}{L by p (R=1) or L by p by R (R>1) array of posterior
#'   means. Per-effect coefficients can be computed as
#'   \code{alpha * mu / X_column_scale_factors}.}
#'
#' \item{mu2_diag}{L by p (R=1) or L by p by R (R>1) array of the
#'   diagonal of the posterior second moment matrix.}
#'
#' \item{pi}{Prior inclusion probabilities (length p vector).}
#'
#' \item{Xr}{N by R matrix of fitted values on the standardized scale.}
#'
#' \item{KL}{Vector of single-effect KL divergences.}
#'
#' \item{lbf}{Vector of single-effect log-Bayes factors.}
#'
#' \item{sigma2}{Residual variance (R by R matrix for R > 1).}
#'
#' \item{V}{Prior variance scalar (per effect).}
#'
#' \item{converged}{Logical indicating whether the algorithm converged.}
#'
#' \item{elbo}{Vector of ELBO values at each iteration.}
#'
#' \item{niter}{Number of iterations performed.}
#'
#' \item{sets}{Estimated credible sets.}
#'
#' \item{pip}{Vector of posterior inclusion probabilities.}
#'
#' \item{z}{Matrix of univariate z-scores (when requested).}
#'
#' \item{single_effect_lfsr}{L by R matrix of per-effect lfsr.}
#'
#' \item{lfsr}{J by R matrix of per-variable lfsr.}
#'
#' \item{conditional_lfsr}{L by J by R array of conditional lfsr
#'   (given variable j is the single effect).}
#'
#' \item{lbf_outcome}{L by R matrix of per-outcome conditional log
#'   Bayes factors.  Measures per-outcome evidence from the conditional
#'   (residualized) data, without cross-outcome borrowing.  Used to
#'   filter \code{single_effect_lfsr}: outcomes with \code{lbf_outcome < 0}
#'   (BF < 1) have their lfsr set to 1.}
#'
#' \item{prior_mixture_weights}{Vector of estimated prior mixture
#'   weights across the K covariance components (only with mixture
#'   prior).}
#'
#' \item{posterior_mixture_weights}{L by J by K array of posterior
#'   mixture component assignments (L by J by (K+1) when
#'   \code{null_weight > 0}, with the first slice being the null
#'   component).}
#'
#' \item{V_structure}{List of prior covariance matrices (only with
#'   mixture prior).}
#'
#' @examples
#' # Example with one response.
#' set.seed(1)
#' n <- 2000
#' p <- 1000
#' beta <- rep(0, p)
#' beta[1:4] <- 1
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X %*% beta + rnorm(n)
#' fit <- mvsusie(X, Y, L = 10)
#'
#' # Sufficient statistics example with one response.
#' X_colmeans <- colMeans(X)
#' Y_colmeans <- colMeans(Y)
#' X <- scale(X, center = TRUE, scale = FALSE)
#' Y <- scale(Y, center = TRUE, scale = FALSE)
#' XtX <- crossprod(X)
#' XtY <- crossprod(X, Y)
#' YtY <- crossprod(Y)
#' res <- mvsusie_ss(XtX, XtY, YtY, n, L = 10, X_colmeans, Y_colmeans)
#'
#' # RSS example with one response.
#' R <- crossprod(X)
#' z <- susieR::calc_z(X, Y)
#' res <- mvsusie_rss(z, R, N = n, L = 10)
#'
#' # Example with three responses.
#' set.seed(1)
#' n <- 500
#' p <- 1000
#' true_eff <- 2
#' X <- sample(c(0, 1, 2), size = n * p, replace = TRUE)
#' X <- matrix(X, n, p)
#' beta1 <- rep(0, p)
#' beta2 <- rep(0, p)
#' beta3 <- rep(0, p)
#' beta1[1:true_eff] <- runif(true_eff)
#' beta2[1:true_eff] <- runif(true_eff)
#' beta3[1:true_eff] <- runif(true_eff)
#' y1 <- X %*% beta1 + rnorm(n)
#' y2 <- X %*% beta2 + rnorm(n)
#' y3 <- X %*% beta3 + rnorm(n)
#' Y <- cbind(y1, y2, y3)
#' prior <- create_mixture_prior(R = 3)
#' fit <- mvsusie(X, Y, prior_variance = prior)
#'
#' # Sufficient statistics example with three responses.
#' X_colmeans <- colMeans(X)
#' Y_colmeans <- colMeans(Y)
#' X <- scale(X, center = TRUE, scale = FALSE)
#' Y <- scale(Y, center = TRUE, scale = FALSE)
#' XtX <- crossprod(X)
#' XtY <- crossprod(X, Y)
#' YtY <- crossprod(Y)
#' res <- mvsusie_ss(XtX, XtY, YtY, n,
#'   L = 10, X_colmeans, Y_colmeans,
#'   prior_variance = prior
#' )
#'
#' # RSS example with three responses.
#' R <- crossprod(X)
#' Z <- calc_z(X, Y)
#' res <- mvsusie_rss(Z, R, N = n, L = 10, prior_variance = prior)
#'
#' @export
#'
mvsusie <- function(X, Y, L = 10, prior_variance = 0.2,
                    residual_variance = NULL, prior_weights = NULL,
                    standardize = TRUE, intercept = TRUE,
                    estimate_residual_variance = TRUE,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method = "optim",
                    estimate_prior_mixture_weights = TRUE,
                    mixture_weight_method = "mixsqp",
                    check_null_threshold = 0, prior_tol = 1e-9,
                    model_init = NULL,
                    missing_y_method = "approximate",
                    coverage = 0.95, min_abs_corr = 0.5,
                    compute_univariate_zscore = FALSE,
                    precompute_cache = TRUE, n_thread = 1,
                    max_iter = 100, tol = 1e-3, verbose = TRUE,
                    track_fit = FALSE,
                    min_outcome_lbf = 0) {
  # For R=1 with scalar prior, convert from susieR "scaled prior variance"
  # convention to absolute prior variance: actual_V = scaled_V * sigma2.
  Y_ncol <- if (is.null(dim(Y))) 1L else ncol(Y)
  is_numeric_prior <- is.numeric(prior_variance) && !is.matrix(prior_variance)
  if (Y_ncol == 1L && is_numeric_prior) {
    rv <- if (is.null(residual_variance)) as.numeric(var(Y))
          else as.numeric(residual_variance)
    prior_variance <- prior_variance * rv
  }
  mvsusie_core(X, Y, L = L, prior_variance = prior_variance,
             residual_variance = residual_variance,
             prior_weights = prior_weights,
             standardize = standardize, intercept = intercept,
             estimate_residual_variance = estimate_residual_variance,
             estimate_prior_variance = estimate_prior_variance,
             estimate_prior_method = estimate_prior_method,
             estimate_prior_mixture_weights = estimate_prior_mixture_weights,
             mixture_weight_method = mixture_weight_method,
             check_null_threshold = check_null_threshold,
             prior_tol = prior_tol,
             model_init = model_init,
             missing_y_method = missing_y_method,
             coverage = coverage, min_abs_corr = min_abs_corr,
             compute_univariate_zscore = compute_univariate_zscore,
             precompute_cache = precompute_cache,
             n_thread = n_thread,
             max_iter = max_iter, tol = tol,
             verbose = verbose, track_fit = track_fit,
             min_outcome_lbf = min_outcome_lbf)
}

#' @rdname mvsusie
#'
#' @param Z J x R matrix of z-scores.
#'
#' @param R J x J LD matrix.
#'
#' @param N Sample size.
#'
#' @param Bhat Alternative summary data giving the estimated effects
#'   (J x R matrix). This, together with \code{Shat}, may be
#'   provided instead of \code{Z}.
#'
#' @param Shat Alternative summary data giving the standard errors of
#'   the estimated effects (J x R matrix). This, together with
#'   \code{Bhat}, may be provided instead of \code{Z}.
#'
#' @param varY The sample covariance of Y, defined as \eqn{Y'Y/(N-1)}.
#'   When the sample covariance is not provided, the coefficients
#'   (returned from \code{coef}) are computed on the
#'   \dQuote{standardized} X, y scale.
#'
#' @param prior_variance Can be either (1) a vector of length L, or a
#'   scalar, for scaled prior variance when Y is univariate (which
#'   should then be equivalent to \code{\link[susieR]{susie}}); or (2) a
#'   matrix for a simple multivariate regression; or (3) a mixture prior
#'   from \code{\link{create_mixture_prior}}.
#'
#' @param residual_variance The residual variance
#'
#' @param \dots Additional arguments passed to
#'   \code{\link{mvsusie_ss}}.
#'
#' @export
#'
mvsusie_rss <- function(Z, R, N, Bhat, Shat, varY,
                        prior_variance = 0.2,
                        residual_variance = NULL,
                        estimate_residual_variance = FALSE,
                        ...) {
  is_numeric_prior <-
    !(is.matrix(prior_variance) || class(prior_variance)[1] == "mash_prior")

  if (sum(c(missing(Z), missing(Bhat) || missing(Shat))) != 1) {
    stop("Please provide either Z or (Bhat, Shat), but not both")
  }

  # Check input.
  if (anyNA(R)) {
    stop("Input R matrix contains NAs.")
  }
  if (missing(Z)) {
    J <- ifelse(is.matrix(Bhat), nrow(Bhat), length(Bhat))
  } else {
    J <- ifelse(is.matrix(Z), nrow(Z), length(Z))
  }
  if (nrow(R) != J) {
    stop(paste0(
      "The dimension of R (", nrow(R), " x ", ncol(R), ") does not ",
      "agree with expected (", J, " x ", J, ")"
    ))
  }

  # Check input N.
  if (!missing(N)) {
    if (N <= 1) {
      stop("N must be greater than 1")
    }
  }

  if (missing(Z)) {
    if (length(Shat) == 1) {
      if (is.matrix(Bhat)) {
        Shat <- matrix(Shat, nrow(Bhat), ncol(Bhat))
      } else {
        Shat <- rep(Shat, length(Bhat))
      }
    }
    if (is.matrix(Bhat)) {
      if (nrow(Bhat) != nrow(Shat)) {
        stop("The number of rows of Bhat and Shat do not agree")
      }
      if (ncol(Bhat) != ncol(Shat)) {
        stop("The number of columns of Bhat and Shat do not agree")
      }
    } else {
      if (length(Bhat) != length(Shat)) {
        stop("The length of Bhat and Shat do not agree")
      }
    }
    if (anyNA(Bhat) || anyNA(Shat)) {
      stop("Bhat, Shat cannot have missing values")
    }
    if (any(Shat <= 0)) {
      stop("Shat cannot have zero or negative elements")
    }
    Z <- Bhat / Shat
  }
  if (anyNA(Z)) {
    warning("NA values in z-scores are replaced with 0")
    Z[is.na(Z)] <- 0
  }

  if (!missing(N)) {
    adj <- (N - 1) / (Z^2 + N - 2)
    Z <- sqrt(adj) * Z
  }

  if (estimate_residual_variance && missing(varY)) {
    warning_message(paste0(
      "Estimating residual variance from summary statistics without ",
      "providing varY may be inaccurate. Consider providing varY ",
      "or setting estimate_residual_variance = FALSE."))
  }

  if (!is.null(dim(Z)) && ncol(Z) > 1 && is_numeric_prior) {
    stop("Please specify prior variance for the multivariate z-scores")
  }

  is_numeric_matrix(R, "R")

  # Ensure R (LD matrix) is positive semidefinite.
  # Non-PD LD matrices arise from LD estimates based on different samples,
  # rounding, or subsetting.  Clip negative eigenvalues to zero,
  # reusing susieR's check_semi_pd (cached via zzz.R).
  r_tol <- 1e-08
  semi_pd <- check_semi_pd(R, r_tol)
  if (!semi_pd$status) {
    eig <- attr(semi_pd$matrix, "eigen")
    eig$values[eig$values < 0] <- 0
    R <- eig$vectors %*% (eig$values * t(eig$vectors))
    R <- (R + t(R)) / 2
    warning("LD matrix is not positive semidefinite; ",
            "negative eigenvalues have been set to zero.",
            call. = FALSE)
  }

  if (missing(N)) {
    warning(
      "Providing the sample size (N), or even a rough estimate of N, ",
      "is highly recommended. Without N, the implicit assumption is ",
      "N is large (Inf) and the effect sizes are small (close to zero)."
    )
    if (!missing(varY)) {
      if (!is.null(dim(varY))) varY <- cov2cor(varY)
    } else {
      varY <- estimate_cov_z(Z)
    }
    s <- mvsusie_ss(
      XtX = R, XtY = Z, YtY = varY, N = 2,
      prior_variance = prior_variance,
      standardize = FALSE,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance, ...
    )
  } else {
    # The sample size (N) is provided, so use PVE-adjusted z-scores.
    if (!missing(Shat) & !missing(varY)) {
      # var_y, shat (and bhat) are provided, so the effects are on the
      # *original scale*.
      if (is.null(dim(Z))) {
        XtXdiag <- varY * adj / (Shat^2)
        XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
        XtX <- (XtX + t(XtX)) / 2
        XtY <- Z * sqrt(adj) * varY / Shat
      } else {
        # Multivariate: adj, Shat are J x R; diag(varY) is R-vector.
        # Use t(t(...) * vec) idiom to broadcast the R-vector correctly
        # across columns (outcomes) of the J x R matrix.
        # XtXdiag[j] = mean over r of: varY[r,r] * adj[j,r] / Shat[j,r]^2
        XtXdiag <- rowMeans(t(t(adj / Shat^2) * diag(varY)))
        XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
        XtX <- (XtX + t(XtX)) / 2
        # XtY[j,r] = Z_adj[j,r] * sqrt(adj[j,r]) * varY[r,r] / Shat[j,r]
        # where Z_adj = z_orig * sqrt(adj) (from line 486 above).
        XtY <- t(t(Z * sqrt(adj) / Shat) * diag(varY))
      }
    } else {
      # The effects are on the *standardized* X, y scale.
      XtX <- (N - 1) * R
      XtY <- sqrt(N - 1) * Z
      if (!missing(varY)) {
        if (!is.null(dim(varY))) varY <- cov2cor(varY)
      } else {
        varY <- estimate_cov_z(Z)
      }
    }
    s <- mvsusie_ss(
      XtX = XtX, XtY = XtY, YtY = (N - 1) * varY, N = N,
      prior_variance = prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance, ...
    )
  }

  class(s) <- "mvsusie"
  return(s)
}

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
mvsusie_ss <- function(XtX, XtY, YtY, N, L = 10, X_colmeans = NULL,
                              Y_colmeans = NULL, prior_variance = 0.2,
                              residual_variance = NULL, prior_weights = NULL,
                              standardize = TRUE,
                              estimate_residual_variance = TRUE,
                              estimate_prior_variance = TRUE,
                              estimate_prior_method = "optim",
                              estimate_prior_mixture_weights = TRUE,
                              mixture_weight_method = "mixsqp",
                              check_null_threshold = 0, prior_tol = 1e-9,
                              precompute_cache = TRUE, model_init = NULL,
                              coverage = 0.95, min_abs_corr = 0.5, n_thread = 1,
                              max_iter = 100, tol = 1e-3, verbose = TRUE,
                              track_fit = FALSE,
                              min_outcome_lbf = 0) {
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
  mvsusie_ss_core(XtX, XtY, YtY, N, L = L,
                       X_colmeans = X_colmeans,
                       Y_colmeans = Y_colmeans,
                       prior_variance = prior_variance,
                       residual_variance = residual_variance,
                       prior_weights = prior_weights,
                       standardize = standardize,
                       estimate_residual_variance = estimate_residual_variance,
                       estimate_prior_variance = estimate_prior_variance,
                       estimate_prior_method = estimate_prior_method,
                       estimate_prior_mixture_weights = estimate_prior_mixture_weights,
                       mixture_weight_method = mixture_weight_method,
                       check_null_threshold = check_null_threshold,
                       prior_tol = prior_tol,
                       model_init = model_init,
                       coverage = coverage, min_abs_corr = min_abs_corr,
                       precompute_cache = precompute_cache,
                       n_thread = n_thread,
                       max_iter = max_iter, tol = tol,
                       verbose = verbose, track_fit = track_fit,
                       min_outcome_lbf = min_outcome_lbf)
}

#' @rdname mvsusie
#'
#' @importFrom stats sd var cov cov2cor
#' @importFrom susieR susie_get_cs calc_z
#'
#' @keywords internal
#'
mvsusie_core <- function(X, Y, L = 10, prior_variance = 0.2,
                       residual_variance = NULL, prior_weights = NULL,
                       standardize = TRUE, intercept = TRUE,
                       estimate_residual_variance = TRUE,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       estimate_prior_mixture_weights = TRUE,
                       mixture_weight_method = "mixsqp",
                       check_null_threshold = 0, prior_tol = 1e-9,
                       model_init = NULL,
                       missing_y_method = "approximate",
                       coverage = 0.95, min_abs_corr = 0.5,
                       compute_univariate_zscore = FALSE,
                       precompute_cache = TRUE,
                       n_thread = 1,
                       max_iter = 100, tol = 1e-3, verbose = TRUE,
                       track_fit = FALSE,
                       min_outcome_lbf = 0) {
  reset_warn_once()
  verbose <- isTRUE(verbose)

  # Validate inputs
  if (is.null(dim(Y))) {
    Y <- matrix(Y, length(Y), 1)
  } else if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  R <- ncol(Y)

  # Adjust prior weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / ncol(X), ncol(X))
  } else {
    prior_weights <- prior_weights / sum(prior_weights)
  }

  # Check prior variance
  is_mash_prior <- class(prior_variance)[1] == "mash_prior"
  is_numeric_prior <- !is.matrix(prior_variance) && !is_mash_prior
  if (R > 1 && is_numeric_prior) {
    stop("Please specify prior variance for the multivariate response Y")
  }

  # Scale prior variance when standardizing
  if (standardize && (is.matrix(prior_variance) || is_mash_prior)) {
    sigma <- sapply(seq_len(R), function(i) sd(Y[, i], na.rm = TRUE))
    n <- sapply(seq_len(R), function(i) length(which(!is.na(Y[, i]))))
    sigma <- sigma / sqrt(n)
    if (estimate_prior_variance) {
      sigma <- sigma / max(sigma)
    }
    if (is.matrix(prior_variance)) {
      prior_variance <- scale_covariance(prior_variance, sigma)
    } else {
      prior_variance <- scale_prior_variance.mash_prior(prior_variance, sigma)
    }
  }

  # For R=1 scalar prior, convert to 1x1 matrix for uniform MV treatment
  if (R == 1 && is_numeric_prior) {
    prior_variance <- matrix(prior_variance, 1, 1)
  }

  if (verbose) {
    message(sprintf("mvsusie: N=%d, J=%d, R=%d, L=%d [mem: %.2f GB]",
                    nrow(X), ncol(X), R, L, mem_used_gb()))
  }

  # Validate missing_y_method and apply upfront overrides.
  #
  # Three methods for handling missing Y entries:
  #   "impute"      - Variational imputation (E[Y_miss | Y_obs] at each IBSS
  #                   iteration). Discussed in the mvSuSiE paper
  #                   (doi:10.1038/s41588-025-02486-7). Does not guarantee
  #                   monotone ELBO.
  #   "approximate" - Per-outcome centering with per-pattern V_i^{-1}. Exact
  #                   when V is diagonal or patterns don't overlap.
  #   "exact"       - V^{-1}-weighted centering with full R x R correction.
  #                   Guarantees correct ELBO.
  # The approximate and exact methods are from Zou's PhD thesis, Appendix
  # C.1.2-C.1.3 (http://stephenslab.uchicago.edu/assets/papers/yuxin-thesis.pdf).
  Y_has_missing <- any(is.na(Y))
  missing_y_method <- match.arg(missing_y_method,
                                 c("impute", "approximate", "exact"))
  if (R == 1 || !Y_has_missing) {
    # For R=1 or no missing data, fall back to impute (standard) method
    missing_y_method <- "impute"
  }

  # Classify missingness pattern and warn if method may be suboptimal
  if (Y_has_missing && R > 1) {
    miss_info <- classify_missing_pattern(Y, missing_y_method)
  }

  # When Y has missing data (R > 1) and using approximate/exact methods:
  # - estimate_residual_variance is allowed (ELBO is available for these methods)
  # - Auto-switch exact -> approximate when V is diagonal (they are equivalent)
  if (Y_has_missing && R > 1 &&
      missing_y_method %in% c("approximate", "exact")) {
    if (missing_y_method == "exact" && !is.null(residual_variance) &&
        is.matrix(residual_variance) &&
        all(residual_variance[row(residual_variance) != col(residual_variance)] == 0)) {
      warning_message("Switching to approximate method (equivalent to ",
                      "exact for diagonal residual_variance, but faster).")
      missing_y_method <- "approximate"
    }
  }
  # Note: estimate_residual_variance is allowed for all missing data
  # methods. The update formula E_q[R'R]/n uses expected sufficient
  # statistics (not the ELBO value directly), and for the impute method
  # the Y_cov correction accounts for imputation uncertainty.

  # Create data object
  data <- create_mvsusie_data(X, Y, center = intercept, scale = standardize,
                               missing_y_method = missing_y_method)

  # Set residual variance (also triggers standardize_3d for approximate/exact)
  data <- set_mvsusie_residual_variance(data, residual_variance)
  if (verbose) {
    message("Residual variance set, common_cov=", data$is_common_cov,
            " [mem: ", sprintf("%.2f", mem_used_gb()), " GB]")
  }

  # Compute prior inverses for EM (mixture priors)
  if (is_mash_prior && estimate_prior_variance &&
      estimate_prior_method == "EM") {
    prior_variance <- compute_prior_inv.mash_prior(prior_variance)
  }

  if (verbose && is_mash_prior) {
    message("Prior: K=", length(prior_variance$xUlist),
            " mixture components [mem: ", sprintf("%.2f", mem_used_gb()), " GB]")
  }

  # Convergence method: use ELBO by default, but switch to PIP-based
  # convergence for the imputation missing data method where ELBO is
  # not exact (doi:10.1038/s41588-025-02486-7). The approximate and exact
  # methods compute a valid ELBO so ELBO convergence is appropriate.
  if (Y_has_missing && R > 1 && missing_y_method == "impute") {
    convergence_method <- "pip"
    warning_message(
      "Using PIP-based convergence (instead of ELBO) for impute method: ",
      "the variational imputation does not guarantee a monotone ELBO, ",
      "so PIP convergence is used as a stable alternative.")
  } else {
    convergence_method <- "elbo"
  }

  # Fit model.
  #
  # For 3D missing data with estimate_residual_variance, we use block
  # coordinate ascent over two parameter blocks:
  #   Block A: (alpha, mu, prior_variance) -- optimized by IBSS with fixed sigma2
  #   Block B: sigma2 (residual variance) -- closed-form pairwise estimator
  #
  # Why not update sigma2 within the IBSS loop (as for non-missing data)?
  # In the 3D missing-data path, each missingness pattern k has its own
  # GLS covariance svs_k = V[obs_k, obs_k] / G_j[obs_k, obs_k].  When
  # sigma2 (=V) changes, svs_k changes NON-UNIFORMLY across patterns and
  # variables.  This heterogeneous perturbation destabilizes the prior
  # variance optimizer within the same IBSS iteration, causing wild ELBO
  # oscillations.  In contrast, the non-missing path has svs = V / d[j],
  # a UNIFORM scalar scaling that the optimizer handles smoothly.
  #
  # Block coordinate ascent avoids this by keeping sigma2 fixed during
  # each IBSS run, guaranteeing monotone ELBO.  Sigma2 is then updated
  # from the fully converged posterior, ensuring stable convergence.
  use_outer_rv_loop <- estimate_residual_variance && !is.null(data$miss3d)

  if (use_outer_rv_loop) {
    # Build workhorse argument list (fixed sigma2: est_rv = FALSE)
    wh_args <- list(
      L = L, prior_variance = prior_variance,
      prior_weights = prior_weights,
      estimate_residual_variance = FALSE,
      estimate_prior_variance = estimate_prior_variance,
      estimate_prior_method = estimate_prior_method,
      estimate_prior_mixture_weights = estimate_prior_mixture_weights,
      mixture_weight_method = mixture_weight_method,
      check_null_threshold = check_null_threshold,
      convergence_method = convergence_method,
      max_iter = max_iter, tol = tol,
      prior_tol = prior_tol,
      track_fit = track_fit,
      verbose = verbose,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      precompute_covariances = precompute_cache,
      n_thread = n_thread,
      model_init = model_init)

    # Initial fit with fixed sigma2 (Block A, first pass)
    s <- do.call(mvsusie_workhorse, c(list(data = data), wh_args))

    # Block coordinate step: update sigma2 (Block B), then re-fit (Block A)
    sigma2_block_step <- function(model, data, iter) {
      # Block B: update sigma2 from converged posterior
      sigma2_new <- estimate_residual_variance_3d(data, model)
      data <- set_residual_variance_3d(data, sigma2_new)
      data$residual_variance <- sigma2_new

      # Block A: re-fit with new sigma2, warm-started.
      # Strip V_structure because cleanup_model may have pruned mixture
      # components, making it incompatible with the original prior.
      model_init <- model
      model_init$V_structure <- NULL
      wh_args$model_init <- model_init
      s_new <- do.call(mvsusie_workhorse, c(list(data = data), wh_args))
      s_new$sigma2 <- sigma2_new

      list(model = s_new, data = data,
           log_msg = sprintf("sigma2=[%s]",
                             paste(round(diag(sigma2_new), 3), collapse = ", ")))
    }

    s <- block_coordinate_ascent(s, data, sigma2_block_step,
                                  max_iter = max_iter, tol = tol,
                                  verbose = verbose)

  } else {
    # Standard path: sigma2 updated within IBSS loop
    s <- mvsusie_workhorse(data, L = L, prior_variance = prior_variance,
                            prior_weights = prior_weights,
                            estimate_residual_variance = estimate_residual_variance,
                            estimate_prior_variance = estimate_prior_variance,
                            estimate_prior_method = estimate_prior_method,
                            estimate_prior_mixture_weights = estimate_prior_mixture_weights,
                            mixture_weight_method = mixture_weight_method,
                            check_null_threshold = check_null_threshold,
                            convergence_method = convergence_method,
                            max_iter = max_iter, tol = tol,
                            prior_tol = prior_tol,
                            track_fit = track_fit,
                            verbose = verbose,
                            coverage = coverage,
                            min_abs_corr = min_abs_corr,
                            precompute_covariances = precompute_cache,
                            n_thread = n_thread,
                            model_init = model_init)
  }

  # Compute CSs using original X
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets <- susie_get_cs(s, coverage = coverage, X = X,
                           min_abs_corr = min_abs_corr)
  }

  # Convert to standard mvsusie output format.
  # For missing data 3d path, use per-outcome cm/csd (J x R matrices)
  # from miss3d instead of the J-vector placeholders in data$cm/data$csd.
  if (!is.null(data$miss3d)) {
    out_cm  <- data$miss3d$cm
    out_csd <- data$miss3d$csd
  } else {
    out_cm  <- data$cm
    out_csd <- data$csd
  }
  s <- format_mvsusie_output(s, csd = out_csd, cm = out_cm,
                              Y_mean = data$Y_mean,
                              estimate_prior_variance = estimate_prior_variance,
                              is_ss = FALSE,
                              min_outcome_lbf = min_outcome_lbf)

  # Report z-scores from univariate regression
  if (compute_univariate_zscore) {
    s$z <- calc_z(X, Y, center = intercept, scale = standardize)
  }

  # Set canonical names
  s$outcome_names <- if (is.null(colnames(Y)))
    paste0("outcome", seq_len(R)) else colnames(Y)
  s$variable_names <- if (is.null(colnames(X)))
    paste0("var", seq_len(ncol(X))) else colnames(X)

  # Apply dimnames to match standard output format
  s <- apply_mvsusie_dimnames(s, s$variable_names, s$outcome_names)

  flush_warn_once(verbose = verbose)
  return(s)
}


#' @rdname mvsusie
#'
#' @importFrom stats cov2cor
#' @importFrom susieR susie_get_cs
#'
#' @keywords internal
#'
mvsusie_ss_core <- function(XtX, XtY, YtY, N, L = 10,
                                  X_colmeans = NULL, Y_colmeans = NULL,
                                  prior_variance = 0.2,
                                  residual_variance = NULL,
                                  prior_weights = NULL,
                                  standardize = TRUE,
                                  estimate_residual_variance = TRUE,
                                  estimate_prior_variance = TRUE,
                                  estimate_prior_method = "optim",
                                  estimate_prior_mixture_weights = TRUE,
                                  mixture_weight_method = "mixsqp",
                                  check_null_threshold = 0, prior_tol = 1e-9,
                                  model_init = NULL,
                                  coverage = 0.95, min_abs_corr = 0.5,
                                  precompute_cache = TRUE,
                                  n_thread = 1,
                                  max_iter = 100, tol = 1e-3, verbose = TRUE,
                                  track_fit = FALSE,
                                  min_outcome_lbf = 0) {
  reset_warn_once()

  XtY <- as.matrix(XtY)
  YtY <- as.matrix(YtY)
  R <- ncol(XtY)
  J <- ncol(XtX)

  # Adjust prior weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / J, J)
  } else {
    prior_weights <- prior_weights / sum(prior_weights)
  }

  # Check prior variance
  is_mash_prior <- class(prior_variance)[1] == "mash_prior"
  is_numeric_prior <- !is.matrix(prior_variance) && !is_mash_prior
  if (R > 1 && is_numeric_prior) {
    stop("Please specify prior variance for the multivariate response Y")
  }

  # Scale prior variance when standardizing
  if (standardize && (is.matrix(prior_variance) || is_mash_prior)) {
    sigma <- sqrt(diag(YtY) / (N - 1))
    sigma <- sigma / sqrt(N)
    if (estimate_prior_variance) {
      sigma <- sigma / max(sigma)
    }
    if (is.matrix(prior_variance)) {
      prior_variance <- scale_covariance(prior_variance, sigma)
    } else {
      prior_variance <- scale_prior_variance.mash_prior(prior_variance, sigma)
    }
  }

  # For R=1 scalar prior, convert to 1x1 matrix for uniform MV treatment
  if (R == 1 && is_numeric_prior) {
    prior_variance <- matrix(prior_variance, 1, 1)
  }

  # Compute prior inverses for EM (mixture priors)
  if (is_mash_prior && estimate_prior_variance &&
      estimate_prior_method == "EM") {
    prior_variance <- compute_prior_inv.mash_prior(prior_variance)
  }

  # Default residual variance: sample covariance YtY/(N-1)
  # This matches the individual data path which uses cov(Y).
  if (is.null(residual_variance)) {
    residual_variance <- YtY / (N - 1)
    if (R == 1) {
      residual_variance <- as.numeric(residual_variance)
    }
  }

  # Create SS data object
  data <- create_mvsusie_ss_data(XtX, XtY, YtY, N,
                                  X_colmeans = X_colmeans,
                                  Y_colmeans = Y_colmeans,
                                  standardize = standardize,
                                  residual_variance = residual_variance)

  # Fit model
  s <- mvsusie_workhorse(data, L = L, prior_variance = prior_variance,
                          prior_weights = prior_weights,
                          estimate_residual_variance = estimate_residual_variance,
                          estimate_prior_variance = estimate_prior_variance,
                          estimate_prior_method = estimate_prior_method,
                          estimate_prior_mixture_weights = estimate_prior_mixture_weights,
                          mixture_weight_method = mixture_weight_method,
                          check_null_threshold = check_null_threshold,
                          max_iter = max_iter, tol = tol,
                          prior_tol = prior_tol,
                          track_fit = track_fit,
                          verbose = verbose,
                          coverage = coverage,
                          min_abs_corr = min_abs_corr,
                          precompute_covariances = precompute_cache,
                          n_thread = n_thread,
                          model_init = model_init)

  # Compute CSs using XtX correlation
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets <- susie_get_cs(s, coverage = coverage,
                           Xcorr = cov2cor(XtX),
                           min_abs_corr = min_abs_corr)
  }

  # Store residual variance: use model's updated value if estimated,
  # otherwise fall back to data's initial value
  if (!estimate_residual_variance) {
    s$residual_variance <- data$residual_variance
  }

  # For SS, use colmeans if available; otherwise use zero vectors
  ss_Y_mean <- if (!is.null(Y_colmeans)) Y_colmeans else rep(0, R)
  ss_cm     <- if (!is.null(X_colmeans)) X_colmeans else rep(0, J)

  # Convert to standard mvsusie output format
  s <- format_mvsusie_output(s, csd = data$csd, cm = ss_cm,
                              Y_mean = ss_Y_mean,
                              estimate_prior_variance = estimate_prior_variance,
                              is_ss = TRUE,
                              min_outcome_lbf = min_outcome_lbf)

  # Set canonical names
  s$outcome_names <- if (is.null(colnames(XtY)))
    paste0("outcome", seq_len(R)) else colnames(XtY)
  s$variable_names <- if (is.null(colnames(XtX)))
    paste0("var", seq_len(J)) else colnames(XtX)

  # Apply dimnames to match standard output format
  s <- apply_mvsusie_dimnames(s, s$variable_names, s$outcome_names)

  flush_warn_once(verbose = verbose)
  return(s)
}
