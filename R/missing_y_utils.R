# Utility functions for approximate and exact missing data methods.
#
# These methods are described in Zou's PhD thesis, Appendix C.1.2-C.1.3:
#   http://stephenslab.uchicago.edu/assets/papers/yuxin-thesis.pdf
#
# When Y has missing entries (R > 1), we construct an N x J x R array
# (X_3d) where X[i,,r] = 0 whenever Y[i,r] is missing. This ensures
# that missing outcomes do not contribute to the regression. Each
# observation is assigned to a unique missingness pattern k, and
# per-pattern embedded precision matrices V_i^{-1} are computed by
# inverting the observed sub-block of V and embedding it in R x R.
#
# Two methods are supported:
#   "approximate" - Per-outcome centering. Exact when V is diagonal
#                   or missingness patterns do not overlap across outcomes.
#   "exact"       - V^{-1}-weighted centering with full R x R correction.
#                   Guarantees correct ELBO.
#
# Precomputation strategy: K x J summary matrices (raw_sq_sum, raw_sum)
# are computed once from the raw X, enabling O(K * R^2 * J) svs_inv
# computation instead of O(N * R^2 * J) per-sample loops (K << N).
# XtR is computed per-outcome using the already-centered X_3d.

#' @importFrom matrixStats colSds
NULL

# =============================================================================
# INITIALIZATION
# =============================================================================

#' Initialize the missing data structure for approximate/exact methods.
#'
#' Creates the X_3d array (N x J x R), extracts unique missingness
#' patterns, assigns each observation to its pattern, and precomputes
#' K x J summary matrices for fast svs_inv computation.
#'
#' @param X N x J raw (uncentered, unscaled) covariate matrix.
#' @param Y N x R raw response matrix (NAs at missing entries).
#' @param Y_missing N x R logical matrix (TRUE = missing).
#' @param R Number of outcomes.
#' @param method Character: "approximate" or "exact".
#'
#' @return A list (the miss3d structure) with fields for use
#'   by standardize_3d and set_residual_variance_3d.
#'
#' @keywords internal
init_missing_data_3d <- function(X, Y, Y_missing, R, method) {
  N <- nrow(X)
  J <- ncol(X)

  Y_non_missing <- !Y_missing

  # Extract unique missingness patterns (K x R boolean matrix).
  # Each row is a unique observed-outcome pattern.
  pattern <- unique(Y_non_missing)
  K <- nrow(pattern)

  # Assign each observation to its pattern index (N-vector).
  pattern_assign <- integer(N)
  for (k in seq_len(K)) {
    idx <- which(apply(Y_non_missing, 1, function(x) identical(x, pattern[k, ])))
    pattern_assign[idx] <- k
  }

  # Zero-fill Y at missing entries.
  Y[Y_missing] <- 0

  # Create X_3d: N x J x R array where X[i,,r] = NA when Y[i,r] is missing.
  X_3d <- array(X, dim = c(N, J, R))
  for (r in seq_len(R)) {
    X_3d[Y_missing[, r], , r] <- NA
  }

  # Precompute K x J summary matrices from raw X for fast svs_inv.
  # raw_sq_sum[k,j] = sum_{i in pattern k} X[i,j]^2
  # raw_sum[k,j]    = sum_{i in pattern k} X[i,j]
  # n_k[k]          = number of observations with pattern k
  raw_sq_sum <- matrix(0, K, J)
  raw_sum    <- matrix(0, K, J)
  n_k        <- integer(K)
  for (k in seq_len(K)) {
    idx <- which(pattern_assign == k)
    n_k[k] <- length(idx)
    if (n_k[k] == 1) {
      raw_sq_sum[k, ] <- X[idx, ]^2
      raw_sum[k, ]    <- X[idx, ]
    } else {
      raw_sq_sum[k, ] <- colSums(X[idx, , drop = FALSE]^2)
      raw_sum[k, ]    <- colSums(X[idx, , drop = FALSE])
    }
  }

  list(
    method         = method,
    X_3d           = X_3d,
    Y_non_missing  = Y_non_missing,
    pattern        = pattern,
    pattern_assign = pattern_assign,
    raw_sq_sum     = raw_sq_sum,
    raw_sum        = raw_sum,
    n_k            = n_k,
    standardized   = FALSE,   # flag: has standardize_3d been called?
    # Fields populated by standardize_3d / set_residual_variance_3d:
    cm             = NULL,    # J x R per-outcome column means
    csd            = NULL,    # J x R per-outcome column SDs
    Vinv           = NULL,    # list of K R x R per-pattern embedded precision
    Vinv_eigen     = NULL,    # list of K eigenvalue vectors
    Xbar           = NULL,    # J x R x R (exact only)
    Vinvsuminv     = NULL,    # R x R (exact only)
    elbo_const     = NULL     # scalar
  )
}


# =============================================================================
# STANDARDIZATION
# =============================================================================

#' Apply per-outcome centering and scaling to X_3d and Y.
#'
#' For the approximate method, uses simple per-outcome column means/SDs.
#' For the exact method, uses V^{-1}-weighted centering (requires Vinv
#' to be already computed in data$miss3d$Vinv).
#'
#' Updates data$Y and data$Y_mean with the correctly centered Y.
#'
#' @param data The full data object (must have data$miss3d populated,
#'   and for exact method, data$miss3d$Vinv must be set).
#' @param center Logical: center X and Y.
#' @param scale Logical: scale X columns to unit variance per outcome.
#'
#' @return Modified data object.
#'
#' @keywords internal
standardize_3d <- function(data, center, scale) {
  my <- data$miss3d
  X_3d <- my$X_3d
  N <- dim(X_3d)[1]
  J <- dim(X_3d)[2]
  R <- dim(X_3d)[3]

  # Per-outcome column means (J x R)
  if (center) {
    cm <- matrix(0, J, R)
    for (r in seq_len(R)) {
      cm[, r] <- colMeans(X_3d[, , r], na.rm = TRUE)
    }
  } else {
    cm <- matrix(0, J, R)
  }

  # Per-outcome column SDs (J x R)
  csd <- matrix(1, J, R)
  for (r in seq_len(R)) {
    if (scale) {
      csd[, r] <- colSds(X_3d[, , r], center = cm[, r], na.rm = TRUE)
      csd[csd[, r] == 0, r] <- 1
    }
    # Center and scale X_3d
    X_3d[, , r] <- t((t(X_3d[, , r]) - cm[, r]) / csd[, r])
    # Replace NAs with 0 (missing outcomes don't contribute)
    X_3d[, , r][is.na(X_3d[, , r])] <- 0
  }

  my$cm  <- cm
  my$csd <- csd

  if (my$method == "approximate") {
    my$cm <- cm  # store for intercept recovery
    if (center) {
      # Simple per-outcome Y means (observed entries only)
      if (R == 1) {
        Y_mean <- mean(data$Y[my$Y_non_missing])
      } else {
        Y_mean <- sapply(seq_len(R), function(r) {
          mean(data$Y[my$Y_non_missing[, r], r])
        })
      }
      data$Y <- t(t(data$Y) - Y_mean)
      data$Y[!my$Y_non_missing] <- 0
      data$Y_mean <- Y_mean
    }
  } else {
    # Exact method: V^{-1}-weighted centering
    if (center) {
      Vinv <- my$Vinv
      K <- nrow(my$pattern)
      pattern_assign <- my$pattern_assign

      # Vinvsum = sum_k n_k * Vinv_k  (R x R)
      Vinvsum <- Reduce("+", lapply(seq_len(K), function(k) {
        Vinv[[k]] * my$n_k[k]
      }))
      Vinvsuminv <- invert_via_chol(Vinvsum)$inv
      my$Vinvsuminv <- Vinvsuminv

      # sum_i Vinv_{pattern(i)} %*% Y[i,]  (R-vector)
      Ysum <- rep(0, R)
      for (i in seq_len(N)) {
        Ysum <- Ysum + as.numeric(Vinv[[pattern_assign[i]]] %*% data$Y[i, ])
      }

      # Weighted Y mean
      Y_mean <- as.numeric(Vinvsuminv %*% Ysum)
      data$Y <- t(t(data$Y) - Y_mean)
      data$Y[!my$Y_non_missing] <- 0
      data$Y_mean <- Y_mean

      # Xbar[j,,] = Vinvsuminv %*% sum_i Vinv_{pattern(i)} * X_3d[i,j,:]
      # Optimized: use precomputed raw_sum instead of O(J*N*R^2) loop.
      # For each j: A2[,r] = sum_k Vinv_k[,r] * xsum_kr where
      # xsum_kr = (raw_sum[k,j] - n_k[k]*cm[j,r])/csd[j,r] for obs r.
      # This is O(J*K*R^2) instead of O(J*N*R^2).
      pattern_int <- matrix(as.integer(my$pattern), nrow = nrow(my$pattern),
                             ncol = ncol(my$pattern))
      my$Xbar <- compute_Xbar_from_sums_rcpp(
        Vinv, pattern_int, my$raw_sum, as.integer(my$n_k),
        cm, csd, Vinvsuminv)
    }
  }

  my$X_3d <- X_3d
  my$standardized <- TRUE
  data$miss3d <- my
  return(data)
}


# =============================================================================
# RESIDUAL VARIANCE AND SVS COMPUTATION
# =============================================================================

#' Compute per-pattern precision, svs, and svs_inv for missing data.
#'
#' Called when residual variance is set. On the first call (initialization),
#' this also triggers standardization of X_3d and Y via standardize_3d().
#' The call order is: Vinv -> standardize -> svs_inv, because the exact
#' method's V^{-1}-weighted centering requires Vinv.
#'
#' @param data The full data object with data$miss3d populated.
#' @param residual_variance R x R residual variance matrix.
#' @param center Logical: center X and Y (only used on first call).
#' @param scale Logical: scale X (only used on first call).
#'
#' @return Modified data object with svs, svs_inv set at top level.
#'
#' @keywords internal
set_residual_variance_3d <- function(data, residual_variance,
                                     center = TRUE, scale = TRUE) {
  my <- data$miss3d
  R_dim <- data$R
  J <- data$p
  K <- nrow(my$pattern)

  # --- Step 1: Per-pattern embedded precision matrices ---
  Vinv <- vector("list", K)
  Vinv_eigen <- vector("list", K)

  for (k in seq_len(K)) {
    obs_r <- which(my$pattern[k, ])
    if (R_dim == 1) {
      Vinv[[k]] <- my$pattern[k, ] / residual_variance
      if (length(obs_r) > 0) {
        Vinv_eigen[[k]] <- as.numeric(residual_variance)
      } else {
        Vinv_eigen[[k]] <- numeric(0)
      }
    } else {
      Vinv[[k]] <- matrix(0, R_dim, R_dim)
      if (length(obs_r) > 0) {
        Vk <- residual_variance[obs_r, obs_r, drop = FALSE]
        eigenVk <- eigen(Vk, symmetric = TRUE)
        dinv <- 1 / eigenVk$values
        Vinv[[k]][obs_r, obs_r] <- eigenVk$vectors %*%
          (dinv * t(eigenVk$vectors))
        Vinv_eigen[[k]] <- eigenVk$values
      } else {
        Vinv_eigen[[k]] <- numeric(0)
      }
    }
  }

  my$Vinv <- Vinv
  my$Vinv_eigen <- Vinv_eigen
  data$miss3d <- my

  # --- Step 2: Standardize X_3d and Y (first call only) ---
  # The exact method needs Vinv for V^{-1}-weighted centering, so
  # standardization must happen after Vinv is computed.
  if (!my$standardized) {
    data <- standardize_3d(data, center, scale)
    my <- data$miss3d  # re-read after standardize_3d modified it
  }

  # --- Step 3: ELBO normalization constant ---
  # -0.5 * log(2*pi) * (total observed entries) - 0.5 * sum log|V_k| per pattern
  Y_missing_assign_table <- tabulate(my$pattern_assign, nbins = K)
  my$elbo_const <- -0.5 * log(2 * pi) *
    sum(sapply(Vinv_eigen, length) * Y_missing_assign_table) -
    0.5 * sum(sapply(Vinv_eigen, function(x) {
      if (length(x) > 0) sum(log(x)) else 0
    }) * Y_missing_assign_table)

  data$miss3d <- my

  # --- Step 4: Per-variable svs_inv using precomputed K x J sums ---
  # svs_inv[j] = sum_i diag(x_ij) Vinv_{pattern(i)} diag(x_ij)
  # where x_ij is the per-outcome centered/scaled X_3d.
  #
  # For pattern k with observed outcomes obs_k, expanding the product
  # for both r,s in obs_k:
  #   sum_{i in k} x_ij[r] * x_ij[s]
  #     = (1/(csd[j,r]*csd[j,s])) *
  #       (raw_sq_sum[k,j] - (cm[j,r]+cm[j,s])*raw_sum[k,j] +
  #        n_k[k]*cm[j,r]*cm[j,s])
  cm  <- my$cm
  csd <- my$csd

  # For exact method, precompute Vinvsum (shared across all variables)
  Vinvsum <- NULL
  if (my$method == "exact" && R_dim > 1) {
    Vinvsum <- Reduce("+", lapply(seq_len(K), function(k) {
      Vinv[[k]] * my$n_k[k]
    }))
  }

  # C++ fast path: compute all J svs_inv matrices at once
  method_int <- if (my$method == "exact") 2L else 1L
  Xbar_arg <- if (!is.null(my$Xbar)) my$Xbar else array(0, c(1, R_dim, R_dim))
  Vinvsum_arg <- if (!is.null(Vinvsum)) Vinvsum else matrix(0, R_dim, R_dim)

  pattern_int <- matrix(as.integer(my$pattern), nrow = K, ncol = R_dim)
  svs_inv_3d <- compute_svs_inv_3d_rcpp(
    Vinv, pattern_int,
    my$raw_sq_sum, my$raw_sum, as.integer(my$n_k),
    cm, csd, method_int, Xbar_arg, Vinvsum_arg)

  # Convert R x R x J array to list of J R x R matrices; invert each
  data$svs_inv <- lapply(seq_len(J), function(j) svs_inv_3d[, , j])
  data$svs <- lapply(seq_len(J), function(j) {
    tryCatch(
      invert_via_chol(svs_inv_3d[, , j])$inv,
      error = function(e) {
        invert_via_chol(svs_inv_3d[, , j] + 1e-8 * diag(R_dim))$inv
      }
    )
  })

  # Check whether all variables share the same SVS (for mashr C++ path).
  # With missing data, svs varies across variables, so this is typically FALSE.
  data$is_common_cov <- is_list_common_3d(data$svs)

  return(data)
}


# =============================================================================
# COMPUTE X'R (V^{-1}-weighted)
# =============================================================================

#' Compute V^{-1}-weighted X'R for missing data methods.
#'
#' Uses per-outcome X_3d for the cross-product. First applies
#' per-pattern V_i^{-1} to the residual, then computes the
#' per-outcome cross-products.
#'
#' @param data The full data object.
#' @param R_mat N x R residual matrix.
#'
#' @return J x R matrix of V^{-1}-weighted cross-products.
#'
#' @keywords internal
compute_XtR_3d <- function(data, R_mat) {
  my <- data$miss3d
  method_int <- if (my$method == "exact") 2L else 1L

  # Xbar placeholder for approximate method (C++ needs a cube argument)
  Xbar <- if (!is.null(my$Xbar)) my$Xbar else array(0, c(1, data$R, data$R))

  return(compute_XtR_3d_rcpp(my$X_3d, R_mat, my$Vinv,
                             as.integer(my$pattern_assign),
                             method_int, Xbar))
}


# =============================================================================
# COMPUTE X*b (per-outcome)
# =============================================================================

#' Compute per-outcome X*b for missing data methods.
#'
#' @param data The full data object.
#' @param b J x R coefficient matrix.
#'
#' @return N x R matrix of fitted values.
#'
#' @keywords internal
compute_Xb_3d <- function(data, b) {
  my <- data$miss3d
  if (is.vector(b)) b <- matrix(b, length(b), 1)
  method_int <- if (my$method == "exact") 2L else 1L

  # Xbar placeholder for approximate method
  Xbar <- if (!is.null(my$Xbar)) my$Xbar else array(0, c(1, data$R, data$R))

  return(compute_Xb_3d_rcpp(my$X_3d, b, method_int, Xbar))
}


# =============================================================================
# COMPUTE VinvR (per-pattern V_i^{-1} application)
# =============================================================================

#' Apply per-pattern V_i^{-1} to each row of a matrix.
#'
#' @param data The full data object.
#' @param mat N x R matrix.
#'
#' @return N x R matrix with per-pattern V_i^{-1} applied.
#'
#' @keywords internal
compute_VinvR_3d <- function(data, mat) {
  my <- data$miss3d

  return(compute_VinvR_3d_rcpp(mat, my$Vinv,
                               as.integer(my$pattern_assign)))
}


# =============================================================================
# COMPUTE BETAHAT (GLS estimate)
# =============================================================================

#' Compute GLS betahat for missing data methods.
#'
#' betahat[j,:] = svs[j] %*% XtR[j,:], the generalized least squares
#' estimate using the precomputed per-variable covariance svs.
#'
#' @param data The full data object.
#' @param XtR J x R cross-product matrix.
#'
#' @return J x R matrix of betahat values.
#'
#' @keywords internal
compute_betahat_3d <- function(data, XtR) {
  svs <- data$svs

  # Needs R x R x J array
  svs_3d <- matlist2array(svs)
  return(compute_betahat_3d_rcpp(svs_3d, XtR))
}


# =============================================================================
# INTERCEPT RECOVERY
# =============================================================================

#' Recover intercept on the original Y scale for missing data methods.
#'
#' For the approximate method, uses per-outcome column means.
#' For the exact method, uses V^{-1}-weighted computation with
#' precomputed pattern-level column sums (raw_sum) to avoid
#' storing the raw X matrix.
#'
#' @param data The full data object.
#' @param model The model object.
#'
#' @return R-vector of intercepts.
#'
#' @keywords internal
get_intercept_3d <- function(data, model) {
  my <- data$miss3d
  b_sum <- compute_posterior_mean_sum_from_model(model)
  R_dim <- data$R

  if (R_dim == 1) {
    csd_vec <- as.vector(my$csd)
    cm_vec <- as.vector(my$cm)
    coefs <- b_sum / csd_vec
    if (!is.null(data$Y_mean)) {
      intercept <- data$Y_mean - sum(cm_vec * coefs)
    } else {
      intercept <- 0
    }
    return(intercept)
  }

  coefs <- b_sum / my$csd  # J x R

  if (my$method == "approximate") {
    if (!is.null(data$Y_mean)) {
      intercept <- data$Y_mean - colSums(my$cm * coefs)
    } else {
      intercept <- rep(0, R_dim)
    }
  } else {
    # Exact: V^{-1}-weighted intercept recovery.
    # D = sum_i Vinv_{pattern(i)} %*% (t(coefs) %*% X_raw[i,])
    #
    # Grouping by pattern and using precomputed raw_sum:
    # D = sum_k Vinv_k %*% (t(coefs) %*% raw_sum[k,])
    # where raw_sum[k,] is the J-vector of raw X column sums for pattern k.
    Vinv <- my$Vinv
    K <- nrow(my$pattern)

    D <- rep(0, R_dim)
    for (k in seq_len(K)) {
      # t(coefs) %*% raw_sum[k,] is an R-vector
      coefs_x_rawsum <- as.numeric(crossprod(coefs, my$raw_sum[k, ]))
      D <- D + as.numeric(Vinv[[k]] %*% coefs_x_rawsum)
    }
    if (!is.null(data$Y_mean)) {
      intercept <- data$Y_mean - as.numeric(my$Vinvsuminv %*% D)
    } else {
      intercept <- rep(0, R_dim)
    }
  }

  return(intercept)
}


# =============================================================================
# ELBO COMPUTATION
# =============================================================================

#' Compute ELBO for missing data methods.
#'
#' Uses per-pattern normalization constant and V_i^{-1}-weighted
#' quadratic forms for the expected sum of squared residuals.
#'
#' @param data The full data object.
#' @param model The model object.
#'
#' @return Scalar ELBO value.
#'
#' @keywords internal
compute_elbo_3d <- function(data, model) {
  # Inject updated V-dependent quantities from model (after est_rv update)
  if (!is.null(model$miss3d_elbo_const))
    data$miss3d$elbo_const <- model$miss3d_elbo_const
  if (!is.null(model$miss3d_Vinv))
    data$miss3d$Vinv <- model$miss3d_Vinv

  my <- data$miss3d
  R_dim <- data$R
  L <- nrow(model$alpha)
  N <- dim(my$X_3d)[1]

  # Normalization constant (precomputed)
  loglik <- my$elbo_const

  # Quadratic form: sum_i r_i' V_i^{-1} r_i
  b_sum <- compute_posterior_mean_sum_from_model(model)
  Xb <- compute_Xb_3d(data, b_sum)
  R_mat <- data$Y - Xb
  R_mat[!my$Y_non_missing] <- 0

  VinvR <- compute_VinvR_3d(data, R_mat)
  essr <- sum(R_mat * VinvR)

  # Per-effect correction: subtract first-moment^2 and add second-moment
  for (l in seq_len(L)) {
    b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
    if (!is.matrix(b_l)) b_l <- matrix(b_l, ncol = R_dim)
    Xb_l <- compute_Xb_3d(data, b_l)
    VinvXb_l <- compute_VinvR_3d(data, Xb_l)
    essr <- essr - sum(Xb_l * VinvXb_l)
  }

  # Add vbxxb (precomputed per-effect terms from compute_kl)
  essr <- essr + sum(model$vbxxb)

  result <- loglik - 0.5 * essr
  return(result)
}


# =============================================================================
# EM M-STEP: RESIDUAL VARIANCE FOR 3D MISSING DATA
# =============================================================================

#' Pairwise complete-case residual variance estimator for 3D missing data.
#'
#' Estimates V[r,s] using only samples where both outcomes r and s are
#' observed, analogous to the non-missing estimator but with per-pair
#' normalization.
#'
#' The update is:
#'   V_new[r,s] = (R_hat'R_hat + correction)[r,s] / N_{rs}
#'
#' where:
#'   R_hat    = Y - X_3d E[B] with R_hat[i,r] = 0 when r is missing
#'              for sample i.  crossprod(R_hat)[r,s] automatically sums
#'              only over samples where both r and s are observed.
#'
#'   correction = sum_l (bxxb_3d_l - b1_XtX_b1_l),
#'              the posterior variance correction using per-outcome Gram
#'              matrix G_j (Hadamard product).  G_j[r,s] also only counts
#'              samples where both r and s are observed, so the correction
#'              is on the same scale as R_hat'R_hat.
#'
#'   N_{rs}   = number of samples where both outcome r and outcome s are
#'              observed.  For r = s this equals n_r (the per-outcome count).
#'
#' Why pairwise complete-case (not EM with conditional imputation):
#'   The EM M-step for V with missing data introduces Y_mis as latent and
#'   imputes via Lambda_k = V[mis,obs] V[obs,obs]^{-1}.  This creates a
#'   feedback loop: V_old -> Lambda -> inflated R_complete -> V_new > V_old.
#'   In early IBSS iterations when effects haven't converged, the residuals
#'   still contain signal, and Lambda-amplified imputation inflates V
#'   without bound.
#'
#'   The pairwise complete-case estimator avoids this by using ONLY
#'   observed data.  It is a consistent estimator under MCAR and is
#'   bounded: V_new[r,r] <= Var(Y[,r]) since regression reduces variance.
#'   It matches the non-missing path (which also uses method-of-moments)
#'   when there are no missing values.
#'
#' @param data The full data object with data$miss3d populated.
#' @param model The model object with posterior parameters.
#'
#' @return R x R estimated residual variance matrix.
#'
#' @keywords internal
estimate_residual_variance_3d <- function(data, model) {
  my <- data$miss3d
  R_dim <- data$R
  N <- dim(my$X_3d)[1]
  J <- data$p
  L <- nrow(model$alpha)
  K <- nrow(my$pattern)

  # ---------------------------------------------------------------
  # 1. Residuals at posterior mean, zeros at missing entries
  # ---------------------------------------------------------------
  b_sum <- compute_posterior_mean_sum_from_model(model)
  Xb <- compute_Xb_3d(data, b_sum)
  R_mat <- data$Y - Xb
  R_mat[!my$Y_non_missing] <- 0

  # ---------------------------------------------------------------
  # 2. Observed-data cross-product
  # ---------------------------------------------------------------
  # crossprod(R_mat)[r,s] = sum_i R_mat[i,r] * R_mat[i,s].
  # Since R_mat[i,r] = 0 when outcome r is missing for sample i,
  # this automatically sums only over samples where BOTH r and s
  # are observed.
  E_RtR <- crossprod(R_mat)

  # ---------------------------------------------------------------
  # 4. Posterior correction: E_q[B'X_3d'X_3d B] - E[B]'X_3d'X_3d E[B]
  # ---------------------------------------------------------------
  # model$bxxb uses scalar d[j] = N-1 which is WRONG for 3D path.
  # The correct formula uses the per-outcome Gram matrix G_j:
  #
  #   E_q[(X_3d b_l)' (X_3d b_l)]_{r,s}
  #     = sum_j alpha_lj * G_j[r,s] * E[beta_lr * beta_ls | gamma=j]
  #     = sum_j alpha_lj * G_j[r,s] * mu2_lj[r,s]
  #
  # This is a Hadamard (element-wise) product of G_j and mu2.
  # G_j[r,s] = sum_{i: obs r AND s} X_3d[i,j,r] * X_3d[i,j,s],
  # computed from per-pattern summary statistics.
  #
  # Since the full R x R posterior second moment mu2_lj is not stored
  # (reduced memory mode), we use:
  #   mu2_lj[r,s] = mu_lj[r] * mu_lj[s]  for r != s  (exact for diagonal post_cov)
  #   mu2_lj[r,r] = mu2_diag[j,r]         for r = s   (exact)
  #
  # Note: the correction only has non-zero contributions for observed
  # entries (X_3d = 0 for missing), so it naturally respects the
  # missingness structure.
  cm  <- my$cm
  csd <- my$csd

  # Pre-extract posterior means and diagonal second moments into cubes
  # for the C++ function (L x J x R arrays).
  mu_cube <- model$mu  # L x J x R (already a 3D array)
  if (length(dim(mu_cube)) != 3)
    mu_cube <- array(mu_cube, dim = c(L, J, R_dim))
  mu2_diag_cube <- array(0, dim = c(L, J, R_dim))
  for (l in seq_len(L))
    mu2_diag_cube[l, , ] <- model$mu2_cache[[l]]$mu2_diag

  # Compute total bxxb via C++ (Hadamard product of G_j with mu2,
  # summed over all L effects and J variables).
  pattern_int <- matrix(as.integer(my$pattern), nrow = K, ncol = R_dim)
  total_bxxb <- compute_bxxb_correction_3d_rcpp(
    model$alpha, mu_cube, mu2_diag_cube,
    pattern_int, my$raw_sq_sum, my$raw_sum, as.integer(my$n_k),
    cm, csd)

  # Subtract (X_3d E[b_l])' (X_3d E[b_l]) for each effect (BLAS-efficient)
  b1_XtX_b1_total <- matrix(0, R_dim, R_dim)
  for (l in seq_len(L)) {
    mu_l <- mu_cube[l, , , drop = TRUE]
    if (!is.matrix(mu_l)) mu_l <- matrix(mu_l, ncol = R_dim)
    B1_l <- model$alpha[l, ] * mu_l
    if (!is.matrix(B1_l)) B1_l <- matrix(B1_l, ncol = R_dim)
    XB1_l <- compute_Xb_3d(data, B1_l)
    b1_XtX_b1_total <- b1_XtX_b1_total + crossprod(XB1_l)
  }

  correction <- total_bxxb - b1_XtX_b1_total

  # ---------------------------------------------------------------
  # 5. Pairwise complete-case estimator: V_new[r,s] = numerator / N_{rs}
  # ---------------------------------------------------------------
  # N_rs[r,s] = number of samples where both outcome r and s are observed.
  # Uses the same observation indicator that made R_mat entries zero.
  N_rs <- crossprod(my$Y_non_missing * 1.0)  # R x R
  N_rs[N_rs == 0] <- 1  # avoid division by zero (no shared observations)

  V_new <- (E_RtR + correction) / N_rs
  V_new <- (V_new + t(V_new)) / 2

  # Enforce positive-definiteness (pairwise estimates can be indefinite
  # when missingness patterns are severely unbalanced)
  V_new <- makePD(V_new, 1e-10)

  return(V_new)
}


# =============================================================================
# DISPATCHING HELPERS
# =============================================================================
# Thin wrappers that dispatch to the 3d functions when data$miss3d is
# non-NULL, otherwise fall through to the standard computation.

#' Compute X'R, dispatching to 3d method when missing data is active.
#'
#' @param data The full data object.
#' @param R_mat N x R residual matrix.
#'
#' @return J x R cross-product matrix.
#'
#' @keywords internal
compute_XtR <- function(data, R_mat) {
  if (!is.null(data$miss3d)) {
    compute_XtR_3d(data, R_mat)
  } else {
    crossprod(data$X, R_mat)
  }
}

#' Compute X*b, dispatching to 3d method when missing data is active.
#'
#' @param data The full data object.
#' @param b J x R coefficient matrix.
#'
#' @return N x R fitted values.
#'
#' @keywords internal
compute_Xb <- function(data, b) {
  if (!is.null(data$miss3d)) {
    compute_Xb_3d(data, b)
  } else {
    data$X %*% b
  }
}

#' Compute betahat, dispatching to 3d method when missing data is active.
#'
#' @param data The full data object.
#' @param XtR J x R cross-product matrix (for 3d) or model residuals.
#'
#' @return J x R betahat matrix.
#'
#' @keywords internal
compute_betahat <- function(data, XtR) {
  if (!is.null(data$miss3d)) {
    compute_betahat_3d(data, XtR)
  } else {
    as.matrix(XtR / data$d)
  }
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

# Compute svs_inv for a single variable j. Factored out of the main loop
# in set_residual_variance_3d for parallel execution via future.apply.
#
# @param j Variable index.
# @param method "approximate" or "exact".
# @param R_dim Number of outcomes.
# @param K Number of unique missing patterns.
# @param Vinv List of K R x R per-pattern embedded precision matrices.
# @param obs_r_list List of K integer vectors (observed outcome indices).
# @param cm J x R per-outcome column means.
# @param csd J x R per-outcome column SDs.
# @param raw_sq_sum K x J matrix of squared column sums per pattern.
# @param raw_sum K x J matrix of column sums per pattern.
# @param n_k K-vector of pattern counts.
# @param Xbar J x R x R array (exact only, NULL for approximate).
# @param Vinvsum R x R matrix (exact only, NULL for approximate).
#
# @return R x R svs_inv matrix for variable j.
compute_svs_inv_j <- function(j, method, R_dim, K, Vinv, obs_r_list,
                               cm, csd, raw_sq_sum, raw_sum, n_k,
                               Xbar, Vinvsum) {
  if (method == "approximate") {
    svs_inv_j <- matrix(0, R_dim, R_dim)
    for (k in seq_len(K)) {
      obs_r <- obs_r_list[[k]]
      if (length(obs_r) == 0) next
      for (ri in seq_along(obs_r)) {
        r <- obs_r[ri]
        for (si in ri:length(obs_r)) {
          s <- obs_r[si]
          cross_rs <- (raw_sq_sum[k, j] -
            (cm[j, r] + cm[j, s]) * raw_sum[k, j] +
            n_k[k] * cm[j, r] * cm[j, s]) /
            (csd[j, r] * csd[j, s])
          val <- Vinv[[k]][r, s] * cross_rs
          svs_inv_j[r, s] <- svs_inv_j[r, s] + val
          if (r != s) svs_inv_j[s, r] <- svs_inv_j[s, r] + val
        }
      }
    }
    return(svs_inv_j)
  }

  # Exact method
  if (R_dim == 1) {
    svs_inv_j <- 0
    for (k in seq_len(K)) {
      svs_inv_j <- svs_inv_j + Vinv[[k]] * (
        raw_sq_sum[k, j] -
        2 * cm[j, 1] * raw_sum[k, j] +
        n_k[k] * cm[j, 1]^2) / csd[j, 1]^2
    }
    return(matrix(svs_inv_j, 1, 1))
  }

  # Exact, R > 1
  A1 <- matrix(0, R_dim, R_dim)
  A2 <- matrix(0, R_dim, R_dim)
  for (k in seq_len(K)) {
    obs_r <- obs_r_list[[k]]
    if (length(obs_r) == 0) next
    # A1 contribution
    for (ri in seq_along(obs_r)) {
      r <- obs_r[ri]
      for (si in ri:length(obs_r)) {
        s <- obs_r[si]
        cross_rs <- (raw_sq_sum[k, j] -
          (cm[j, r] + cm[j, s]) * raw_sum[k, j] +
          n_k[k] * cm[j, r] * cm[j, s]) /
          (csd[j, r] * csd[j, s])
        val <- Vinv[[k]][r, s] * cross_rs
        A1[r, s] <- A1[r, s] + val
        if (r != s) A1[s, r] <- A1[s, r] + val
      }
    }
    # A2 contribution
    for (r in obs_r) {
      xsum_r <- (raw_sum[k, j] - n_k[k] * cm[j, r]) / csd[j, r]
      A2[, r] <- A2[, r] + Vinv[[k]][, r] * xsum_r
    }
  }

  Xbar_j <- Xbar[j, , ]
  if (!is.matrix(Xbar_j)) Xbar_j <- matrix(Xbar_j, R_dim, R_dim)

  A1 - crossprod(Xbar_j, A2) - crossprod(A2, Xbar_j) +
    crossprod(Xbar_j, Vinvsum %*% Xbar_j)
}


# Check if all elements of a list of matrices are (approximately) identical.
# Used for is_common_cov check with missing data.
is_list_common_3d <- function(lst, tol = 1e-10) {
  if (length(lst) <= 1) return(TRUE)
  ref <- lst[[1]]
  ref_norm <- max(abs(ref))
  if (ref_norm == 0) ref_norm <- 1
  all(vapply(lst[-1], function(x) max(abs(x - ref)) / ref_norm < tol, logical(1)))
}
