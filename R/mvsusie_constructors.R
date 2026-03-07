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
create_mvsusie_data <- function(X, Y, center = TRUE, scale = TRUE) {
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
  if (R == 1) {
    Y_mean <- mean(Y[obs, 1])
  } else {
    Y_mean <- colMeans(Y, na.rm = TRUE)
  }

  # Center
  if (center) {
    X <- t(t(X) - cm)
    Y <- t(t(Y) - Y_mean)
  } else {
    cm <- rep(0, J)
    Y_mean <- rep(0, R)
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

  # X'X diagonal = colSums(X^2) using observed rows only when Y has
  # missing data; missing rows don't contribute to regression.
  d <- colSums(X[obs, , drop = FALSE]^2)
  d[d == 0] <- 1e-6

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
    Y_has_missing  = Y_has_missing,
    Y_missing      = Y_missing,
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

  # Auto-compute if NULL

  if (is.null(residual_variance)) {
    if (R > 1) {
      if (!data$Y_has_missing) {
        residual_variance <- cov(data$Y)
      } else {
        stop("residual_variance must be provided when Y has missing data")
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

  # Whether all variables share the same SVS (optimization for mashr C++)
  data$is_common_cov <- (length(unique(data$d)) == 1)

  # Precompute per-variable covariance matrices
  if (precompute_covariances) {
    data$svs <- lapply(seq_len(data$p), function(j) {
      res <- residual_variance / data$d[j]
      res[is.nan(res) | is.infinite(res)] <- 1e6
      res
    })
    data$svs_inv <- lapply(seq_len(data$p), function(j) {
      data$residual_variance_inv * data$d[j]
    })
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

  d <- diag(XtX)
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
    Y_has_missing = FALSE,
    residual_variance     = NULL,
    residual_variance_inv = NULL,
    svs     = NULL,
    svs_inv = NULL
  )

  if (!is.null(residual_variance)) {
    data <- set_mvsusie_residual_variance(data, residual_variance)
  }

  class(data) <- c("mv_ss", "ss")
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
    # V_structure stores the prior covariance matrix unnormalized (matching R6
    # convention). V_scalar starts at 1 and is updated by EM/optim.
    V_structure <- list(prior_var)
    pi_V <- 1.0
    null_weight <- 0
    V_scalar_init <- 1
  } else if (class(prior_var)[1] == "mash_prior") {
    # Mixture prior from create_mash_prior.
    # Store component matrices unnormalized (matching R6 convention).
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

  # Precompute pseudo-inverses for EM
  inv_list <- lapply(V_structure, pseudo_inverse)
  V_structure_inv <- matlist2array(lapply(inv_list, `[[`, "inv"))
  V_structure_rank <- vapply(inv_list, `[[`, numeric(1), "rank")

  model <- list(
    alpha        = matrix(1 / J, L, J),
    mu           = array(0, c(L, J, R)),
    mu2          = vector("list", L),
    V            = rep(V_scalar_init, L),
    # Mixture prior fields (unified for single-matrix and mixture)
    V_structure       = V_structure,       # list of K RxR matrices (non-null)
    V_structure_3d    = matlist2array(V_structure),  # RxRxK array
    V_structure_inv   = V_structure_inv,   # RxRxK pseudo-inverses
    V_structure_rank  = V_structure_rank,  # K-vector of ranks
    null_weight       = null_weight,       # scalar
    pi_V              = pi_V,              # K-vector (non-null weights)
    pi_V_posterior    = vector("list", L), # per-effect Jx(K+1) mixture weights
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
    model$mu2[[l]] <- array(0, c(J, R, R))
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
  list(Xr = data$X %*% b)
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
