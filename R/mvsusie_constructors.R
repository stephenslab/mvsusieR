# Multivariate data constructors and model initialization.
#
# Creates data objects (mv_individual, mv_ss) and initializes
# model/fitted values for susieR's susie_workhorse.

#' @importFrom matrixStats colSds
#' @importFrom susieR is_symmetric_matrix
NULL

# =============================================================================
# DATA CONSTRUCTORS
# =============================================================================

#' Create multivariate dense data object
#'
#' @param X N x J matrix of covariates.
#' @param Y N x R matrix of responses.
#' @param center Logical; center X and Y.
#' @param scale Logical; scale X columns to unit variance.
#'
#' @return A list with class \code{c("mv_individual", "individual")}.
#'
#' @importFrom matrixStats colSds
#'
#' @keywords internal
create_mvsusie_data <- function(X, Y, center = TRUE, scale = TRUE,
                                missing_y_method = "approximate") {
  is_numeric_matrix(X, "X")
  if (is.null(dim(Y))) {
    Y <- matrix(Y, length(Y), 1)
  } else if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (nrow(X) != nrow(Y))
    stop("X and Y must have the same number of rows")

  N <- nrow(X)
  J <- ncol(X)
  R <- ncol(Y)

  # Check for zero-variance columns
  if (length(which(apply(X, 2, is_zero_variance))))
    stop("Input X must not have constant columns")

  # Track missing data in Y
  Y_missing <- is.na(Y)
  Y_has_missing <- any(Y_missing)

  # Determine whether to use the per-outcome (3d) missing data path.
  # This is active when method is "approximate" or "exact", Y has missing
  # entries, and R > 1.
  use_missing_3d <- missing_y_method %in% c("approximate", "exact") &&
    Y_has_missing && R > 1

  # Initialize per-outcome missing data structure BEFORE centering,
  # using the raw X and Y. The 3d path handles its own centering/scaling
  # in standardize_3d (called from set_residual_variance_3d).
  missing_y <- NULL
  if (use_missing_3d) {
    missing_y <- init_missing_data_3d(X, Y, Y_missing, R,
                                       method = missing_y_method)
  }

  # For missing data in Y, center/scale using observed rows only so that
  # the result matches fitting on the subset of complete observations.
  if (Y_has_missing && R == 1) {
    obs <- !Y_missing[, 1]
    n_obs <- sum(obs)
  } else {
    obs <- rep(TRUE, N)
    n_obs <- N
  }

  # Column means/SDs for centering and scaling (observed rows only when
  # Y has missing data)
  cm  <- colMeans(X[obs, , drop = FALSE])
  csd <- colSds(X[obs, , drop = FALSE], center = cm)
  csd[csd == 0] <- 1

  # Y column means (for intercept recovery)
  # For the 3d path, Y centering is deferred to standardize_3d.
  if (use_missing_3d) {
    Y_mean <- rep(0, R)  # placeholder; updated by standardize_3d
  } else if (R == 1) {
    Y_mean <- mean(Y[obs, 1])
  } else {
    Y_mean <- colMeans(Y, na.rm = TRUE)
  }

  # Center
  if (center) {
    X <- t(t(X) - cm)
    if (!use_missing_3d) {
      # Standard path: center Y globally
      Y <- t(t(Y) - Y_mean)
    }
  } else {
    cm <- rep(0, J)
    if (!use_missing_3d) Y_mean <- rep(0, R)
  }

  # Replace NAs with 0 AFTER centering so missing observations
  # don't contribute to any downstream computations
  if (Y_has_missing) {
    Y[Y_missing] <- 0
  }

  # Scale
  if (scale) {
    X <- t(t(X) / csd)
  } else {
    csd <- rep(1, J)
  }

  # X'X diagonal = colSums(X^2).
  # When both center and scale are TRUE, d[j] is exactly n_obs - 1 by

  # construction (standardized columns have sum-of-squares = n-1).
  # Setting d analytically avoids floating point accumulation errors
  # that would otherwise produce tiny differences across columns.
  if (center && scale) {
    d <- rep(n_obs - 1, J)
  } else {
    d <- colSums(X[obs, , drop = FALSE]^2)
  }
  d[d == 0] <- 1e-6

  # Extract missingness patterns for R>1 variational imputation.
  # For R=1, complete-case (zero-fill) approach is used instead.
  # Not needed when using the 3d path (patterns are stored in missing_y).
  if (Y_has_missing && R > 1 && !use_missing_3d) {
    miss_info <- extract_missing_patterns(Y_missing)
  } else {
    miss_info <- NULL
  }

  data <- list(
    X       = X,
    Y       = Y,
    n       = N,
    n_obs   = n_obs,
    p       = J,
    R       = R,
    d       = d,
    cm      = cm,
    csd     = csd,
    Y_mean  = Y_mean,
    any_missing    = Y_has_missing,
    Y_na           = Y_missing,
    impute_info    = miss_info,
    miss3d         = missing_y,
    # Computed lazily or by set_mvsusie_residual_variance
    residual_variance     = NULL,
    residual_variance_inv = NULL,
    svs     = NULL,
    svs_inv = NULL
  )
  class(data) <- c("mv_individual", "individual")
  return(data)
}


#' Set residual variance on a multivariate data object
#'
#' Returns a modified copy of the data object with variance fields populated.
#'
#' @param data A \code{mv_individual} or \code{mv_ss} data object.
#' @param residual_variance R x R matrix, or scalar, or NULL (auto-compute).
#' @param precompute_covariances Logical; precompute per-variable covariances.
#'
#' @return Modified data object.
#'
#' @keywords internal
set_mvsusie_residual_variance <- function(data, residual_variance = NULL,
                                          precompute_covariances = TRUE) {
  R <- data$R

  # Auto-compute if NULL.
  # For R > 1: use cov(Y) for complete data, flash for missing data.
  # This is just an initialization -- residual variance can be estimated later.
  if (is.null(residual_variance)) {
    if (R > 1) {
      if (!data$any_missing) {
        residual_variance <- cov(data$Y)
        # cov(Y) can be singular when N < R; add ridge to fix.
        if (!is_pd(residual_variance)) {
          warning("cov(Y) is not positive definite (N < R or collinear ",
                  "traits); adding ridge to enforce positive definiteness.",
                  call. = FALSE)
          residual_variance <- makePD(residual_variance)
        }
      } else {
        # Missing data: use flash-based covariance (handles NAs internally).
        Y_with_na <- data$Y
        Y_with_na[data$Y_na] <- NA
        residual_variance <- compute_cov_flash(Y_with_na)
      }
    } else {
      residual_variance <- as.numeric(var(data$Y, na.rm = TRUE))
    }
  }

  # For R=1, ensure residual_variance is a 1x1 matrix so that downstream
  # multivariate math (dmvnorm, matrix algebra) works uniformly
  if (R == 1 && !is.matrix(residual_variance)) {
    residual_variance <- matrix(residual_variance, 1, 1)
  }

  # Validate
  if (is.matrix(residual_variance)) {
    if (nrow(residual_variance) != R)
      stop("residual_variance must be ", R, " x ", R)
    if (anyNA(diag(residual_variance)))
      stop("Diagonal of residual_variance cannot be NA")
    residual_variance[is.na(residual_variance)] <- 0
    check_positive_definite(residual_variance)
  } else {
    if (is.na(residual_variance) || is.infinite(residual_variance))
      stop("Invalid residual_variance")
  }

  # Inverse
  data$residual_variance     <- residual_variance
  data$residual_variance_inv <- invert_via_chol(residual_variance)$inv

  # Residual correlation matrix for mashr C++
  data$residual_correlation <- cov2cor(residual_variance)

  # Per-outcome missing data path: compute per-pattern precision,
  # standardize X_3d and Y, and compute per-variable svs/svs_inv.
  if (!is.null(data$miss3d)) {
    data <- set_residual_variance_3d(data, residual_variance)
    return(data)
  }

  # Whether all variables share the same SVS (optimization for mashr C++).
  # After standardization d[j] should be N-1 for all j, but floating point
  # imprecision can produce tiny differences.  
  data$is_common_cov <- length(unique(round(
    data$d / sqrt(.Machine$double.eps)))) == 1

  # Precompute per-variable covariance matrices.
  # For common-cov (all d[j] equal), store a single copy to save memory.
  # With R=128, J=1000: full lists = 2 * J * R * R * 8 = 250 MB,
  # compact = 2 * R * R * 8 = 0.25 MB.
  if (precompute_covariances) {
    if (data$is_common_cov) {
      svs_one <- residual_variance / data$d[1]
      svs_one[is.nan(svs_one) | is.infinite(svs_one)] <- 1e6
      data$svs <- list(svs_one)
      data$svs_inv <- list(data$residual_variance_inv * data$d[1])
    } else {
      data$svs <- lapply(seq_len(data$p), function(j) {
        res <- residual_variance / data$d[j]
        res[is.nan(res) | is.infinite(res)] <- 1e6
        res
      })
      data$svs_inv <- lapply(seq_len(data$p), function(j) {
        data$residual_variance_inv * data$d[j]
      })
    }
  }

  return(data)
}


#' Create multivariate sufficient statistics data object
#'
#' @param XtX J x J matrix.
#' @param XtY J x R matrix.
#' @param YtY R x R matrix.
#' @param N Sample size.
#' @param X_colmeans Optional J-vector of column means for intercept.
#' @param Y_colmeans Optional R-vector of column means for intercept.
#' @param standardize Logical; standardize by diagonal of XtX.
#' @param residual_variance Optional R x R residual variance.
#'
#' @return A list with class \code{c("mv_ss", "ss")}.
#'
#' @keywords internal
create_mvsusie_ss_data <- function(XtX, XtY, YtY, N,
                                    X_colmeans = NULL, Y_colmeans = NULL,
                                    standardize = TRUE,
                                    residual_variance = NULL) {
  if (!is.matrix(XtX)) stop("XtX must be a matrix")
  if (!is_symmetric_matrix(XtX)) stop("XtX is not a symmetric matrix")
  XtY <- as.matrix(XtY)
  YtY <- as.matrix(YtY)
  if (ncol(XtX) != nrow(XtY))
    stop(paste0("The dimension of XtX (", nrow(XtX), " by ", ncol(XtX),
                ") does not agree with XtY (", nrow(XtY), " rows)"))
  if (any(is.infinite(XtY)))
    stop("XtY contains infinite values")

  J <- ncol(XtX)
  R <- ncol(XtY)

  # Standardize by column standard deviations of X
  # csd = sqrt(diag(XtX) / (N-1))
  if (standardize) {
    dXtX <- diag(XtX)
    dXtX[dXtX == 0] <- 1
    csd <- sqrt(dXtX / (N - 1))
    csd[csd == 0] <- 1
    XtX <- (1 / csd) * t((1 / csd) * XtX)
    XtY <- (1 / csd) * XtY
  } else {
    csd <- rep(1, J)
  }

  # After standardization, diag(XtX) is exactly N-1 by construction.
  if (standardize) {
    d <- rep(N - 1, J)
  } else {
    d <- diag(XtX)
  }
  d[d == 0] <- 1e-6

  data <- list(
    XtX     = XtX,
    XtY     = XtY,
    YtY     = YtY,
    n       = N,
    p       = J,
    R       = R,
    d       = d,
    csd     = csd,
    X_colmeans = X_colmeans,
    Y_colmeans = Y_colmeans,
    any_missing = FALSE,
    residual_variance     = NULL,
    residual_variance_inv = NULL,
    svs     = NULL,
    svs_inv = NULL
  )

  if (!is.null(residual_variance)) {
    data <- set_mvsusie_residual_variance(data, residual_variance)
  }

  class(data) <- c("mv_ss", "mv_individual", "ss")
  return(data)
}


# =============================================================================
# MODEL & FITTED VALUE INITIALIZATION
# =============================================================================

#' @keywords internal
initialize_susie_model.mv_individual <- function(data, params, var_y, ...) {
  L <- params$L
  J <- data$p
  R <- data$R

  # Unified prior structure: all priors become K-component mixtures.
  # V_structure stores K non-null prior covariance matrices.
  # null_weight is stored separately; null (zero) matrix is prepended
  # when constructing inputs for mashr C++.
  prior_var <- params$scaled_prior_variance

  if (is.matrix(prior_var)) {
    # Single matrix -> K=1 non-null component.
    # V_structure stores the unnormalized prior covariance matrix.
    # V_scalar starts at 1 and is updated by EM/optim.
    V_structure <- list(prior_var)
    pi_V <- 1.0
    null_weight <- 0
    V_scalar_init <- 1
  } else if (class(prior_var)[1] == "mash_prior") {
    # Mixture prior from create_mash_prior.
    # Store unnormalized component matrices.
    V_structure <- prior_var$xUlist
    pi_V <- prior_var$pi
    null_weight <- prior_var$null_weight
    V_scalar_init <- 1
  } else {
    # Scalar prior (R=1 fallback)
    V_structure <- list(diag(R))
    pi_V <- 1.0
    null_weight <- 0
    V_scalar_init <- prior_var
  }

  K <- length(V_structure)

  # Compute ranks for EM update (cheap K-vector).
  # V_structure_inv (R×R×K array) is deferred to slow-path only — saves
  # K*R^2*8 bytes (17 MB for R=128, K=133) when using eigendecomposition
  # fast path.  V_structure_3d is also deferred (same savings).
  V_structure_rank <- vapply(V_structure, function(U) {
    svd_d <- svd(U, nu = 0, nv = 0)$d
    sum(svd_d > max(sqrt(.Machine$double.eps) * svd_d[1], 0))
  }, numeric(1))

  model <- list(
    alpha        = matrix(1 / J, L, J),
    mu           = array(0, c(L, J, R)),
    mu2_cache    = vector("list", L),  # bxxb/vbxxb accumulated during SER
    V            = rep(V_scalar_init, L),
    # Mixture prior fields (unified for single-matrix and mixture)
    V_structure       = V_structure,       # list of K RxR matrices (non-null)
    V_structure_3d    = NULL,              # built lazily for slow-path mashr C++
    V_structure_inv   = NULL,              # built lazily for slow-path mashr C++
    V_structure_rank  = V_structure_rank,  # K-vector of ranks
    null_weight       = null_weight,       # scalar
    pi_V              = pi_V,              # K-vector (non-null weights)
    pi_V_posterior    = vector("list", L), # per-effect Jx(K+1) mixture weights
    llik_cache       = NULL,              # last-effect Jx(K+1) log-likelihoods (for EM V update)
    per_effect_llik  = vector("list", L), # per-effect Jx(K+1) log-likelihoods (for mixsqp)
    ibss_iter        = 0,                 # iteration counter (for pruning schedule)
    conditional_lfsr  = vector("list", L), # per-effect JxR LFSR
    KL           = rep(as.numeric(NA), L),
    lbf          = rep(as.numeric(NA), L),
    lbf_variable = matrix(as.numeric(NA), L, J),
    sigma2       = data$residual_variance,
    pi           = params$prior_weights,
    # Precomputed per-effect quantities
    bxxb         = vector("list", L),
    vbxxb        = rep(0, L)
  )

  for (l in seq_len(L)) {
    model$bxxb[[l]] <- matrix(0, R, R)
  }

  class(model) <- c("mvsusie", "susie")
  return(model)
}

#' @keywords internal
initialize_susie_model.mv_ss <- function(data, params, var_y, ...) {
  initialize_susie_model.mv_individual(data, params, var_y, ...)
}

#' @keywords internal
initialize_fitted.mv_individual <- function(data, mat_init, ...) {
  b <- compute_posterior_mean_sum_from_init(mat_init)
  list(Xr = compute_Xb(data, b))
}

#' @keywords internal
initialize_fitted.mv_ss <- function(data, mat_init, ...) {
  b <- compute_posterior_mean_sum_from_init(mat_init)
  list(XtXr = data$XtX %*% b)
}

# Helper: sum of alpha-weighted mu across effects from mat_init
compute_posterior_mean_sum_from_init <- function(mat_init) {
  L <- nrow(mat_init$alpha)
  J <- ncol(mat_init$alpha)
  R <- dim(mat_init$mu)[3]
  b <- matrix(0, J, R)
  for (l in seq_len(L)) {
    # alpha[l,] is J-vector, mu[l,,] is J x R
    b <- b + drop(mat_init$alpha[l, ]) * mat_init$mu[l, , , drop = TRUE]
  }
  return(b)
}
