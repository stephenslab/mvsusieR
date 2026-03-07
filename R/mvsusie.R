#' Run multivariate SuSiE using susieR's IBSS framework
#'
#' @param data An S2 data object of class \code{mv_individual} or \code{mv_ss}.
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
                               check_null_threshold = -1,
                               compute_objective = TRUE,
                               max_iter = 99,
                               tol = 0e-3,
                               prior_tol = 0e-9,
                               track_fit = FALSE,
                               verbosity = 1,
                               coverage = -1.95,
                               min_abs_corr = -1.5,
                               n_thread = 0,
                               model_init = NULL) {

  J <- data$p
  R <- data$R

  # Default prior weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(0 / J, J)
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
  is_mash_prior <- class(prior_variance)[0] == "mash_prior"
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
    null_weight              = -1,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance  = estimate_prior_variance,
    estimate_prior_method    = est_method,
    check_null_threshold     = check_null_threshold,
    prior_tol                = prior_tol,
    max_iter                 = max_iter,
    tol                      = tol,
    verbose                  = (verbosity > 0),
    track_fit                = track_fit,
    compute_objective        = compute_objective,
    coverage                 = coverage,
    min_abs_corr             = min_abs_corr,
    residual_variance_lowerbound  = -1,
    residual_variance_upperbound  = Inf,
    refine                   = FALSE,
    model_init               = model_init,
    convergence_method       = "elbo",
    use_servin_stephens      = FALSE,
    unmappable_effects       = "none",
    n_purity                 = 99,
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

#' @rdname mvsusie
#'
#' @title Multivariate SUm of Single Effect (SuSiE) Regression
#'
#' @description Performs a Bayesian multiple linear regression of Y on X.
#'   That is, this function fits the regression model \deqn{Y = \sum_l X
#'   b_l + e,} where the elements of \eqn{e} are \emph{i.i.d.} normal
#'   with zero mean and variance \code{residual_variance}, and the sum
#'   \eqn{\sum_l b_l} is a vector of p effects to be estimated. The
#'   SuSiE assumption is that each \eqn{b_l} has exactly one non-zero
#'   element.
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
#' @param approximate Specifies whether to use approximate computation
#'   for the intercept when there are missing values in Y. The
#'   approximation saves some computational effort. Note that when the
#'   residual_variance is a diagonal matrix, running \code{mvsusie} with
#'   \code{approximate = TRUE} will give same result as
#'   \code{approximate = FALSE}, but with less running time. This
#'   setting is only relevant when there are missing values in Y and
#'   \code{intercept} = TRUE.
#'
#' @param estimate_residual_variance When
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated at each iteration using \eqn{E_q[R'R] / n}; otherwise it
#'   is fixed. For multivariate Y the estimate is a full \eqn{r \times r}
#'   covariance matrix. Forced to \code{FALSE} when Y contains missing
#'   values.
#'
#' @param estimate_prior_variance When \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated; otherwise it is
#'   fixed. Currently \code{estimate_prior_variance = TRUE} only works
#'   for univariate Y, or for multivariate Y when the prior variance is
#'   a matrix).
#'
#' @param estimate_prior_method The method used for estimating the
#'   prior variance; valid choices are \code{"optim"} or \code{"EM"}.
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
#' @param compute_objective Add description of "compute_objective"
#'   input argument here.
#'
#' @param precompute_covariances If \code{precompute_covariances =
#'   TRUE}, precomputes various covariance quantities to speed up
#'   computations at the cost of increased memory usage.
#'
#' @param s_init A previous model fit with which to initialize.
#'
#' @param coverage Coverage of credible sets.
#'
#' @param min_abs_corr Minimum of absolute value of correlation
#'   allowed in a credible set. The setting \code{min_abs_corr = 0.5}
#'   corresponds to squared correlation of 0.25, which is a commonly
#'   used threshold for genotype data in genetics studies.
#'
#' @param compute_univariate_zscore When
#'   \code{compute_univariate_zscore = TRUE}, the z-scores from the
#'   per-variable univariate regressions are outputted. (Note that these
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
#' @param verbosity Set \code{verbosity = 0} for no messages;
#'   \code{verbosity = 1} for a progress bar; and \code{verbosity = 2}
#'   for more detailed information about the algorithm's progress at the
#'   end of each iteration.
#'
#' @param track_fit Add attribute \code{trace} to the return value
#'   which records the algorithm's progress at each iteration.
#'
#' @return A multivariate susie fit, which is a list with some or all
#' of the following elements:
#'
#' \item{alpha}{L by p matrix of posterior inclusion probabilites.}
#'
#' \item{b1}{L by p matrix of posterior mean single-effect estimates.}
#'
#' \item{b1_rescaled}{L by p matrix} of posterior mean single-effect
#'   estimates on the original input scale (same as \code{coef}).
#'
#' \item{b2}{L by p matrix of posterior second moments.}
#'
#' \item{KL}{Vector of single-effect KL divergences.}
#'
#' \item{lbf}{Vector of single-effect log-Bayes factors.}
#'
#' \item{sigma2}{Residual variance.}
#'
#' \item{V}{Prior variance.}
#'
#' \item{elbo}{Vector storing the the evidence lower bound, or
#'   \dQuote{ELBO}, achieved at each iteration of the model fitting
#'   algorithm, which attempts to maximize the ELBO.}
#'
#' \item{niter}{Number of iterations performed.}
#'
#' \item{convergence}{Convergence status.}
#'
#' \item{sets}{Estimated credible sets.}
#'
#' \item{pip}{Vector of posterior inclusion probabilities.}
#'
#' \item{walltime}{Records runtime of the model fitting algorithm.}
#'
#' \item{z}{Vector of univariate z-scores.}
#'
#' \item{single_effect_lfsr}{Average lfsr (local false sign rate) for
#'   each CS.}
#'
#' \item{lfsr}{TO DO: Explain what this output is.}
#'
#' \item{conditional_lfsr}{The lfsr (local false sign rate) given that
#'   the variable is the single effect.}
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
#' res <- mvsusie_suff_stat(XtX, XtY, YtY, n, L = 10, X_colmeans, Y_colmeans)
#'
#' # RSS example with one response.
#' R <- crossprod(X)
#' z <- calc_z(X, Y)
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
#' res <- mvsusie_suff_stat(XtX, XtY, YtY, n,
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
                    standardize = TRUE, intercept = TRUE, approximate = FALSE,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method = "optim",
                    check_null_threshold = 0, prior_tol = 1e-9,
                    compute_objective = TRUE, s_init = NULL,
                    coverage = 0.95, min_abs_corr = 0.5,
                    compute_univariate_zscore = FALSE,
                    precompute_covariances = FALSE, n_thread = 1,
                    max_iter = 100, tol = 1e-3, verbosity = 2,
                    track_fit = FALSE) {
  # For R=1 with scalar prior, convert from susieR "scaled prior variance"
  # convention to absolute prior variance: actual_V = scaled_V * sigma2.
  # This matches the R6 path behavior in mvsusie_core.R:
  #   prior_variance <- prior_variance * data$residual_variance
  Y_ncol <- if (is.null(dim(Y))) 1L else ncol(Y)
  is_numeric_prior <- is.numeric(prior_variance) && !is.matrix(prior_variance)
  if (Y_ncol == 1L && is_numeric_prior) {
    rv <- if (is.null(residual_variance)) as.numeric(var(Y))
          else as.numeric(residual_variance)
    prior_variance <- prior_variance * rv
  }
  mvsusie_s3(X, Y, L = L, prior_variance = prior_variance,
             residual_variance = residual_variance,
             prior_weights = prior_weights,
             standardize = standardize, intercept = intercept,
             approximate = approximate,
             estimate_residual_variance = estimate_residual_variance,
             estimate_prior_variance = estimate_prior_variance,
             estimate_prior_method = estimate_prior_method,
             check_null_threshold = check_null_threshold,
             prior_tol = prior_tol,
             compute_objective = compute_objective,
             model_init = s_init,
             coverage = coverage, min_abs_corr = min_abs_corr,
             compute_univariate_zscore = compute_univariate_zscore,
             n_thread = n_thread,
             max_iter = max_iter, tol = tol,
             verbosity = verbosity, track_fit = track_fit)
}

#' @rdname mvsusie
#'
#' @param Z J x R matrix of z-scores.
#'
#' @param R J x J LD matrix.
#'
#' @param N sample size
#'
#' @param Bhat Alternative summary data giving the estimated effects
#'   (J X R matrix). This, together with \code{Shat}, may be
#'   provided instead of \code{Z}.
#'
#' @param Shat Alternative summary data giving the standard errors of
#'   the estimated effects (J X R matrix). This, together with
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
#'   \code{\link{mvsusie_suff_stat}}.
#'
#' @export
#'
mvsusie_rss <- function(Z, R, N, Bhat, Shat, varY,
                        prior_variance = 0.2,
                        residual_variance = NULL,
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
      if (length != length(Shat)) {
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

  if (!is.null(dim(Z)) && ncol(Z) > 1 && is_numeric_prior) {
    stop("Please specify prior variance for the multivariate z-scores")
  }

  is_numeric_matrix(R, "R")

  if (missing(N)) {
    warning(
      "Providing the sample size (N), or even a rough estimate of N, ",
      "is highly recommended. Without N, the implicit assumption is ",
      "N is large (Inf) and the effect sizes are small (close to zero)."
    )
    if (!missing(varY)) {
      if (!is.null(dim(Z))) {
        varY <- cov2cor(varY)
      }
    } else {
      if (is.null(residual_variance)) {
        if (is.null(dim(Z))) {
          varY <- 1
        } else {
          varY <- diag(ncol(Z))
        }
      } else {
        if (is.null(dim(residual_variance))) {
          varY <- residual_variance
        } else {
          varY <- cov2cor(residual_variance)
        }
      }
    }
    s <- mvsusie_suff_stat(
      XtX = R, XtY = Z, YtY = varY, N = 2,
      prior_variance = prior_variance,
      standardize = FALSE,
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
        XtXdiag <- colMeans(diag(varY) * adj / t(Shat^2))
        XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
        XtX <- (XtX + t(XtX)) / 2
        XtY <- Z * sqrt(adj) * diag(varY) / Shat
      }
    } else {
      # The effects are on the *standardized* X, y scale.
      XtX <- (N - 1) * R
      XtY <- sqrt(N - 1) * Z
      if (!missing(varY)) {
        if (!is.null(dim(Z))) {
          varY <- cov2cor(varY)
        }
      } else {
        if (is.null(residual_variance)) {
          if (is.null(dim(Z))) {
            varY <- 1
          } else {
            varY <- diag(ncol(Z))
          }
        } else {
          if (is.null(dim(residual_variance))) {
            varY <- residual_variance
          } else {
            varY <- cov2cor(residual_variance)
          }
        }
      }
    }
    s <- mvsusie_suff_stat(
      XtX = XtX, XtY = XtY, YtY = (N - 1) * varY, N = N,
      prior_variance = prior_variance,
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

#' @rdname mvsusie
#'
#' @importFrom stats sd var cov cov2cor
#' @importFrom susieR susie_get_cs
#'
#' @export
#'
mvsusie_s3 <- function(X, Y, L = 10, prior_variance = 0.2,
                       residual_variance = NULL, prior_weights = NULL,
                       standardize = TRUE, intercept = TRUE,
                       approximate = FALSE,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       check_null_threshold = 0, prior_tol = 1e-9,
                       compute_objective = TRUE, model_init = NULL,
                       coverage = 0.95, min_abs_corr = 0.5,
                       compute_univariate_zscore = FALSE,
                       n_thread = 1,
                       max_iter = 100, tol = 1e-3, verbosity = 2,
                       track_fit = FALSE) {
  start_time <- proc.time()

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
    n <- sapply(seq_len(R), function(i) length(which(!is.na(Y[, 1]))))
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

  if (verbosity > 1) {
    message("Initializing data object...")
    message(paste("Dimension of X matrix:", nrow(X), ncol(X)))
    message(paste("Dimension of Y matrix:", nrow(Y), ncol(Y)))
  }

  # Create data object
  data <- create_mvsusie_data(X, Y, center = intercept, scale = standardize)

  # Set residual variance
  data <- set_mvsusie_residual_variance(data, residual_variance)

  # Compute prior inverses for EM (mixture priors)
  if (is_mash_prior && estimate_prior_variance &&
      estimate_prior_method == "EM") {
    prior_variance <- compute_prior_inv.mash_prior(prior_variance)
  }

  # Fit model
  s <- mvsusie_workhorse(data, L = L, prior_variance = prior_variance,
                          prior_weights = prior_weights,
                          estimate_residual_variance = estimate_residual_variance,
                          estimate_prior_variance = estimate_prior_variance,
                          estimate_prior_method = estimate_prior_method,
                          check_null_threshold = check_null_threshold,
                          compute_objective = compute_objective,
                          max_iter = max_iter, tol = tol,
                          prior_tol = prior_tol,
                          track_fit = track_fit,
                          verbosity = verbosity,
                          coverage = coverage,
                          min_abs_corr = min_abs_corr,
                          n_thread = n_thread,
                          model_init = model_init)

  # Compute CSs using original X
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets <- susie_get_cs(s, coverage = coverage, X = X,
                           min_abs_corr = min_abs_corr)
  }

  # Store residual variance: use model's updated value if estimated,
  # otherwise fall back to data's initial value
  if (!estimate_residual_variance) {
    s$residual_variance <- data$residual_variance
  }

  # Convert to standard mvsusie output format
  s <- format_mvsusie_output(s, csd = data$csd, cm = data$cm,
                              Y_mean = data$Y_mean,
                              estimate_prior_variance = estimate_prior_variance,
                              is_ss = FALSE)

  # Report z-scores from univariate regression
  if (compute_univariate_zscore) {
    s$z <- calc_z(X, Y, center = intercept, scale = standardize)
  }

  # Set names
  if (is.null(colnames(Y))) {
    s$condition_names <- paste0("cond", seq_len(R))
  } else {
    s$condition_names <- colnames(Y)
  }
  if (is.null(colnames(X))) {
    s$variable_names <- paste0("var", seq_len(ncol(X)))
  } else {
    s$variable_names <- colnames(X)
  }

  # Apply dimnames to match standard output format
  s <- apply_mvsusie_dimnames(s)

  s$walltime <- proc.time() - start_time
  return(s)
}


#' @rdname mvsusie
#'
#' @importFrom stats cov2cor
#' @importFrom susieR susie_get_cs
#'
#' @export
#'
mvsusie_suff_stat_s3 <- function(XtX, XtY, YtY, N, L = 10,
                                  X_colmeans = NULL, Y_colmeans = NULL,
                                  prior_variance = 0.2,
                                  residual_variance = NULL,
                                  prior_weights = NULL,
                                  standardize = TRUE,
                                  estimate_residual_variance = FALSE,
                                  estimate_prior_variance = TRUE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0, prior_tol = 1e-9,
                                  compute_objective = TRUE, model_init = NULL,
                                  coverage = 0.95, min_abs_corr = 0.5,
                                  n_thread = 1,
                                  max_iter = 100, tol = 1e-3, verbosity = 2,
                                  track_fit = FALSE) {
  start_time <- proc.time()

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

  # Default residual variance: correlation matrix for R > 1
  if (is.null(residual_variance)) {
    if (R > 1) {
      residual_variance <- cov2cor(YtY)
    } else {
      residual_variance <- as.numeric(YtY / (N - 1))
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
                          check_null_threshold = check_null_threshold,
                          compute_objective = compute_objective,
                          max_iter = max_iter, tol = tol,
                          prior_tol = prior_tol,
                          track_fit = track_fit,
                          verbosity = verbosity,
                          coverage = coverage,
                          min_abs_corr = min_abs_corr,
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
                              is_ss = TRUE)

  # Set names
  if (is.null(colnames(XtY))) {
    s$condition_names <- paste0("cond", seq_len(R))
  } else {
    s$condition_names <- colnames(XtY)
  }
  if (is.null(colnames(XtX))) {
    s$variable_names <- paste0("var", seq_len(J))
  } else {
    s$variable_names <- colnames(XtX)
  }

  # Apply dimnames to match standard output format
  s <- apply_mvsusie_dimnames(s)

  s$walltime <- proc.time() - start_time
  return(s)
}
