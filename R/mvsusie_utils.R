# Utility functions for mvsusieR.
#
# Includes matrix operations, numerical helpers, and lfsr computation functions.

# Report R process memory usage (GB). Uses gc() which is cheap.
mem_used_gb <- function() {
  gc_info <- gc(verbose = FALSE, reset = FALSE)
  sum(gc_info[, "(Mb)"]) / 1024
}

# ============================================================================
# Warning deduplication: emit each warning type once, summarize at end.
#
# During fitting, numerical fallbacks (Cholesky ridge, SVD pseudo-inverse)
# can fire hundreds of times. We emit the warning once for each type, then
# report counts at the end of fitting via message().
# ============================================================================
.mvsusie_warn_state <- new.env(parent = emptyenv())
.mvsusie_warn_state$counts <- list()

# Emit a warning the first time id is seen; silently count subsequent calls.
warn_once <- function(id, ..., call. = FALSE) {
  msg <- paste0(...)
  entry <- .mvsusie_warn_state$counts[[id]]
  if (is.null(entry)) {
    .mvsusie_warn_state$counts[[id]] <- list(msg = msg, count = 1L)
    warning(msg, call. = call.)
  } else {
    .mvsusie_warn_state$counts[[id]]$count <- entry$count + 1L
  }
}

# Reset warning counts (call at start of fitting).
reset_warn_once <- function() {
  .mvsusie_warn_state$counts <- list()
}

# Report suppressed warning counts and reset (call at end of fitting).
# Set verbose = FALSE to suppress the summary messages (verbosity = 0).
flush_warn_once <- function(verbose = TRUE) {
  if (verbose) {
    for (id in names(.mvsusie_warn_state$counts)) {
      entry <- .mvsusie_warn_state$counts[[id]]
      if (entry$count > 1L) {
        warning_message(sprintf(
          "\"%s\" occurred %d times total (%d suppressed)",
          entry$msg, entry$count, entry$count - 1L))
      }
    }
  }
  reset_warn_once()
}

# Cholesky decomposition with automatic ridge fallback.
#
# Tries plain Cholesky first, suppressing rank-deficiency warnings.
# If Cholesky fails, applies makePD and retries. Issues an R warning
# when ridge is applied.
safe_chol <- function(x, ...) {
  res <- tryCatch(
    withCallingHandlers(chol(x, ...),
      warning = function(w) {
        if (grepl("rank-deficient or indefinite", w$message))
          invokeRestart("muffleWarning")
      }),
    error = function(e) NULL
  )
  if (!is.null(res)) return(res)
  x_pd <- makePD(x)
  warn_once("safe_chol_ridge",
            "Cholesky failed; adding ridge ",
            signif(attr(x_pd, "ridge"), 3), " to diagonal")
  chol(x_pd, ...)
}

# Invert a symmetric, positive definite square matrix via its Cholesky
# decomposition. Falls back to SVD pseudo-inverse if Cholesky fails.
# Uses plain chol() (not safe_chol) so that truly singular matrices
# correctly fall through to SVD pseudo-inverse with proper rank.
invert_via_chol <- function(x) {
  if (all(x == 0)) {
    return(list(inv = x, rank = 0))
  }
  ch <- tryCatch(
    withCallingHandlers(chol(x),
      warning = function(w) {
        if (grepl("rank-deficient or indefinite", w$message))
          invokeRestart("muffleWarning")
      }),
    error = function(e) NULL
  )
  if (!is.null(ch)) {
    return(list(inv = chol2inv(ch), rank = nrow(x)))
  }
  # Cholesky failed -> fall back to SVD pseudo-inverse
  warn_once("chol_svd_fallback",
            "Cholesky failed for ", nrow(x), "x", ncol(x),
            " matrix; falling back to SVD pseudo-inverse")
  pseudo_inverse(x)
}

# Log-determinant from upper Cholesky factor: log(det(A)) = 2 * sum(log(diag(chol(A)))).
# Borrowed from mr.mash (misc.R).
chol2ldet <- function(R) {
  2 * sum(log(diag(R)))
}

# Log-determinant of a symmetric PSD matrix via eigenvalues.
# Handles rank-deficient matrices (e.g. when N < R and sigma2 is
# re-estimated) without requiring Cholesky factorization.
log_det_sym <- function(x) {
  eig <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  sum(log(pmax(eig, .Machine$double.xmin)))
}

# Check if a matrix is positive definite (Cholesky succeeds).
is_pd <- function(x) {
  if (anyNA(x)) return(FALSE)
  tryCatch({
    chol(x)
    TRUE
  }, error = function(e) FALSE)
}

# Add ridge to diagonal for positive-definiteness enforcement.
#
# By default uses a scale-relative ridge: max(|diag(x)|) * sqrt(eps).
# An explicit absolute ridge can be passed via e.
# The applied ridge value is attached as attr(result, "ridge").
makePD <- function(x, e = NULL) {
  if (is.null(e))
    e <- max(abs(diag(x))) * sqrt(.Machine$double.eps)
  result <- x + diag(nrow(x)) * e
  attr(result, "ridge") <- e
  result
}

# Pseudoinverse of matrix.
pseudo_inverse <- function(x, tol = sqrt(.Machine$double.eps)) {
  xsvd <- svd(x)
  Positive <- xsvd$d > max(tol * xsvd$d[1L], 0)
  if (all(Positive)) {
    xinv <- xsvd$v %*% (1 / xsvd$d * t(xsvd$u))
  } else {
    xinv <- xsvd$v[, Positive, drop = FALSE] %*%
      ((1 / xsvd$d[Positive]) * t(xsvd$u[, Positive, drop = FALSE]))
  }
  return(list(inv = xinv, rank = sum(Positive)))
}

# Convert a list of matrices to array without losing dimension.
matlist2array <- function(l) {
  if (!inherits(l, "list")) {
    return(l)
  }
  l <- simplify2array(l)
  if (is.null(dim(l))) {
    l <- array(l, c(1, 1, length(l)))
  }
  return(l)
}

# Cannot use "unique" directly here --- for perfectly identical rows
# (by computation) due to possible numerical issues, "unique" and
# "duplicated" function report that they are not identical.
almost.unique <- function(x, tolerance = sqrt(.Machine$double.eps), ...) {
  if (is.matrix(x)) {
    y <- round(x / tolerance, 0)
  } else {
    y <- lapply(1:length(x), function(i) round(x[[i]] / tolerance, 0))
  }
  d <- duplicated(y, ...)
  if (is.matrix(x)) {
    return(x[!d, , drop = FALSE])
  } else {
    return(x[!d])
  }
}

# Duplicated function with a tolerance.
almost.duplicated <- function(x, tolerance = sqrt(.Machine$double.eps), ...) {
  y <- round(x / tolerance, 0)
  return(duplicated(y, ...))
}

# Check if matrix has constant columns.
is_zero_variance <- function(x) {
  if (length(unique(x)) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Faster way to implement diag(sigma) %*% mat %*% diag(sigma)
scale_covariance <- function(mat, sigma) {
  t(mat * sigma) * sigma
}

# Check if input is numeric matrix.
is_numeric_matrix <- function(X, name) {
  if (!((is.double(X) || is.integer(X)) & is.matrix(X))) {
    stop(paste("Input", name, "must be a numeric matrix."))
  }
  if (anyNA(X)) {
    stop(paste("Input", name, "must not contain any missing values."))
  }
  if (any(dim(X) == 0)) {
    stop(paste("Input", name, "dimension is invalid."))
  }
  return(NULL)
}

#' @title Computes the z-scores (t-statistics) for association
#'   between Y and each column of X.
#'
#' @param X N by J matrix of covariates.
#'
#' @param Y Vector of length N, or N by R matrix of response
#'   variables.
#'
#' @param center If \code{center = TRUE}, center X and Y.
#' 
#' @param scale If \code{scale = TRUE}, scale X and Y.
#'
#' @return A matrix of z-scores.
#' 
#' @importFrom susieR univariate_regression
#'
#' @export
#' 
calc_z = function (X, Y, center = FALSE, scale = FALSE) {
  univariate_z = function(X,Y,center,scale) {
    out = univariate_regression(X,Y,center = center,scale = scale)
    return(out$betahat/out$sebetahat)
  }
  if (is.null(dim(Y)))
    return(univariate_z(X,Y,center,scale))
  else
    return(do.call(cbind,lapply(1:ncol(Y),
                                function(i) univariate_z(X,Y[,i],
                                                         center = center,
                                                         scale = scale))))
}

# Weighted log-sum-exp (softmax).
#
# Given log-values and weights, compute the normalized weights
# (softmax) and the log of their weighted sum.
#
# @param value A numeric vector of log-values (or values if log = FALSE).
# @param weight A numeric vector of weights (same length as value).
# @param log If TRUE (default), value is already on log scale.
#
# @return A list with \code{weights} (normalized) and \code{log_sum}.
#
# @keywords internal
compute_softmax <- function(value, weight, log = TRUE) {
  if (length(value) != length(weight))
    stop("Values and their weights should have equal length")
  if (!log) value <- log(value)
  mvalue <- max(value)
  w <- exp(value - mvalue)
  w_weighted <- w * weight
  weighted_sum_w <- sum(w_weighted)
  list(weights = as.vector(w_weighted / weighted_sum_w),
       log_sum = log(weighted_sum_w) + mvalue)
}

# Compute P(variable j | component k) for EM update of prior variance scalar.
#
# Matrix extension of compute_softmax: applies softmax row-by-row
# over variables for each mixture component.
#
# @param prior_variable_weights J-vector of variable-level prior weights.
# @param llik J x (K+1) log-likelihood matrix.
#
# @return (K+1) x J matrix of P(j|k) posterior weights.
#
# @keywords internal
compute_variable_posterior_weights <- function(prior_variable_weights, llik) {
  lbf <- t(llik - llik[, 1])               # (K+1) x J
  lfactors <- apply(lbf, 1, max)            # (K+1)-vector
  d <- t(prior_variable_weights * t(exp(lbf - lfactors)))  # (K+1) x J
  d / rowSums(d)
}

#' @title Create mash prior object.
#'
#' @param mixture_prior a list of (weights = vector(), matrices =
#'   list()) where matrices is a list of prior matrices and have same
#'   length as weights.
#'
#' @param R number of traits
#'
#' @param null_weight whether or not to add a weight for null in
#'   single effect models. By default it takes the null weight from
#'   fitted_g if available. Use \code{null_weight = 0} to override this.
#'
#' @param weights_tol Filter out mixture components with weights
#' smaller than \code{weights_tol}.
#'
#' @param max_mixture_len Only keep the top priors by weight so that
#'   the list of mixture prior is of length \code{max_mixture_len}. Use
#'   \code{max_mixture_len = -1} to include all input weights after
#'   filtering by \code{weights_tol}.
#'
#' @param include_indices Post-process input prior to only include
#'   outcomes at these indices.
#'
#' @param \dots other parameters, for mvsusieR:::create_cov_canonical
#'
#' @return mash prior object for use with mvsusie() function
#'
#' @details Add details here.
#'
#' @examples
#' # Add examples here.
#'
#' @importFrom stats cov2cor
#' @importFrom stats setNames
#' @importFrom stats pnorm
#'
#' @export
#'
create_mixture_prior <- function(mixture_prior, R, null_weight = NULL,
                                 weights_tol = 1e-10,
                                 max_mixture_len = -1, include_indices = NULL, ...) {
  if (sum(c(missing(mixture_prior), missing(R))) != 1) {
    stop("Require exactly one of mixture_prior and R")
  }
  if (is.null(null_weight)) {
    null_weight <- 0
  }
  if (!missing(mixture_prior)) {
    if (!("matrices" %in% names(mixture_prior))) {
      stop("mixture_prior must contain 'matrices'.")
    }
    if (is.null(mixture_prior$weights)) {
      mixture_prior$weights <- rep(
        1 / length(mixture_prior$matrices),
        length(mixture_prior$matrices)
      )
    }
    return(create_mash_prior(
      xUlist = mixture_prior$matrices,
      prior_weights = mixture_prior$weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_outcomes = include_indices
    ))
  }
  if (!missing(R)) {
    Ulist <- create_cov_canonical(R, ...)
    Ulist <- lapply(Ulist, function(mat) {
                           rownames(mat) <- include_indices
                           colnames(mat) <- include_indices
                           return(mat)
                         })
    weights <- rep(1 / length(Ulist), length(Ulist))
    weights <- setNames(weights, names(Ulist)) 
    if (max_mixture_len < length(Ulist) && max_mixture_len > 0) {
      stop(paste0(
        "Automatically generated uniform mixture prior is of ",
        "length ", length(Ulist), " and is greater than currently ",
        "specified max_mixture_len ", max_mixture_len,
        ". Please set max_mixture_len = -1 to allow using all ",
        "of them (although computational speed will suffer)."
      ))
    }
    return(create_mash_prior(
      xUlist = Ulist,
      prior_weights = weights,
      null_weight = null_weight,
      weights_tol = weights_tol,
      top_mixtures = max_mixture_len,
      include_outcomes = include_indices
    ))
  }
}

#' @title Compute List of Canonical Covariance Matrices
#'
#' @description This function computes canonical covariance matrices
#'   to be provided to mash.
#'
#' @param R Integer specifying the number of outcomes.
#'
#' @param singletons \code{TRUE} or \code{FALSE} indicating whether
#'   the singleton matrices are computed.
#'
#' @param hetgrid A vector of numbers between -1 and 1, each
#'   representing the off-diagonal elements of matrices with 1s on the
#'   diagonal. If 0 is included, the identity matrix will be returned
#'   which corresponds to assuming effects are independent across
#'   outcomes. IF \code{hetgrid = NULL}, these matrices are not
#'   returned.
#'
#' @return A list of canonical covariance matrices.
#'
#' @examples
#' mvsusieR:::create_cov_canonical(3)
#' mvsusieR:::create_cov_canonical(3, singletons = FALSE)
#' mvsusieR:::create_cov_canonical(3, hetgrid = NULL)
#'
#' @keywords internal
#'
create_cov_canonical <- function(R, singletons = TRUE,
                                 hetgrid = seq(0, 1, 0.25)) {
  mats <- list()
  nms <- vector("double")
  s_idx <- 0

  # Singleton matrices.
  if (singletons) {
    for (i in 1:R) {
      mats[[i]] <- matrix(0, R, R)
      mats[[i]][i, i] <- 1
      nms[i] <- paste0("singleton_", i)
    }
    s_idx <- R
  }

  # Heterogeneity matrices.
  if (!is.null(hetgrid)) {
    for (j in 1:length(hetgrid)) {
      mats[[s_idx + j]] <- matrix(1, R, R)
      mats[[s_idx + j]][lower.tri(mats[[s_idx + j]], diag = FALSE)] <- hetgrid[j]
      mats[[s_idx + j]][upper.tri(mats[[s_idx + j]], diag = FALSE)] <- hetgrid[j]
      nms[s_idx + j] <- paste0("shared_", j)
    }
  }
  names(mats) <- nms
  return(mats)
}

# ============================================================================
# Eigendecomposition-based precomputation for O(R^2) SER evaluation
#
# The key decomposition: given SVS (R x R, SPD) and U (R x R, PSD),
#   L = chol(SVS)  (lower triangular)
#   M = L^{-1} U L^{-T} = P D P'  (eigendecomposition)
#   Q = L^{-T} P,  G = L P
#
# Then for any scalar V:
#   (SVS + V*U)^{-1} = Q diag(1/(1+V*d)) Q'
#   log|SVS + V*U| = log|SVS| + sum(log(1 + V*d))
#   posterior mean = G diag(V*d/(1+V*d)) Q' bhat
#   posterior cov  = G diag(V*d/(1+V*d)) G'
# ============================================================================

# Eigendecompose L^{-1} U L^{-T} for a single (SVS, U) pair.
#
# @param SVS  R x R symmetric positive definite matrix
# @param U    R x R symmetric positive semidefinite matrix (prior structure)
#
# @return list(Q, G, eigenvalues, log_det_svs)
#   Q: R x R matrix (L^{-T} P), used for Mahalanobis distances
#   G: R x R matrix (L P), used for posterior reconstruction
#   eigenvalues: R-vector of non-negative eigenvalues
#   log_det_svs: log determinant of SVS
eigendecompose_one_pair <- function(SVS, U) {
  R <- nrow(SVS)

  # Cholesky: R's chol() returns upper triangular; we need lower.
  L_upper <- safe_chol(SVS)
  L <- t(L_upper)
  log_det_svs <- chol2ldet(L_upper)

  # L^{-1} U L^{-T} via triangular solves
  L_inv <- forwardsolve(L, diag(R))
  M <- L_inv %*% U %*% t(L_inv)
  M <- (M + t(M)) / 2  # symmetrize for numerical stability

  # Eigendecompose
  eig <- eigen(M, symmetric = TRUE)
  P <- eig$vectors
  d <- pmax(eig$values, 0)  # clamp negative eigenvalues to 0

  Q <- t(L_inv) %*% P  # L^{-T} P  (R x R)
  G <- L %*% P          # L P       (R x R)

  list(Q = Q, G = G, eigenvalues = d, log_det_svs = log_det_svs)
}

# Precompute eigendecomposition cache for all (SVS, U_k) pairs.
#
# @param svs         List of J R x R matrices (or length 1 for common_cov)
# @param V_structure List of K R x R matrices (prior structure)
# @param is_common_cov Logical: are all SVS_j identical?
#
# @return list(is_common_cov, log_det_svs, components)
#   components: length-K list, each with Q, G, eigenvalues
precompute_eigen_cache <- function(svs, V_structure, is_common_cov,
                                    max_cache_gb = 8) {
  K <- length(V_structure)

  if (is_common_cov) {
    SVS <- svs[[1]]
    components <- vector("list", K)
    for (k in seq_len(K)) {
      components[[k]] <- eigendecompose_one_pair(SVS, V_structure[[k]])
    }
    list(
      is_common_cov = TRUE,
      log_det_svs   = components[[1]]$log_det_svs,
      components    = components
    )
  } else {
    # Memory guard: non-common-cov stores K cubes of R x R x J each for
    # Q and G arrays.  For large R/J/K this can exceed available memory.
    J <- length(svs)
    R <- nrow(svs[[1]])
    cache_bytes <- K * 2 * R * R * J * 8  # Q + G cubes
    cache_gb <- cache_bytes / 1e9
    if (cache_gb > max_cache_gb) {
      message("Eigen cache would use ~", round(cache_gb, 1),
              " GB (K=", K, ", R=", R, ", J=", J,
              "); skipping precomputation (falling back to mashr C++).")
      return(NULL)
    }
    # C++ fast path: batch K*J eigendecompositions with Cholesky caching
    return(precompute_eigen_cache_non_common_cpp(svs, V_structure))
  }
}

# R fallback for non-common-cov eigendecomposition (for testing)
precompute_eigen_cache_R <- function(svs, V_structure, is_common_cov) {
  K <- length(V_structure)

  if (is_common_cov) {
    SVS <- svs[[1]]
    components <- vector("list", K)
    for (k in seq_len(K)) {
      components[[k]] <- eigendecompose_one_pair(SVS, V_structure[[k]])
    }
    list(
      is_common_cov = TRUE,
      log_det_svs   = components[[1]]$log_det_svs,
      components    = components
    )
  } else {
    J <- length(svs)
    R <- nrow(svs[[1]])
    log_det_svs <- numeric(J)
    components <- vector("list", K)
    for (k in seq_len(K)) {
      Q_arr   <- array(0, c(R, R, J))
      G_arr   <- array(0, c(R, R, J))
      eig_mat <- matrix(0, R, J)
      for (j in seq_len(J)) {
        decomp <- eigendecompose_one_pair(svs[[j]], V_structure[[k]])
        Q_arr[, , j]  <- decomp$Q
        G_arr[, , j]  <- decomp$G
        eig_mat[, j]  <- decomp$eigenvalues
        if (k == 1) log_det_svs[j] <- decomp$log_det_svs
      }
      components[[k]] <- list(Q = Q_arr, G = G_arr, eigenvalues = eig_mat)
    }
    list(
      is_common_cov = FALSE,
      log_det_svs   = log_det_svs,
      components    = components
    )
  }
}

# Compute J x (K+1) log-likelihood matrix using precomputed eigendecomposition.
#
# Column 1 = null component N(0, SVS).
# Columns 2..(K+1) = non-null components N(0, SVS + V*U_k).
#
# For common_cov, all operations are vectorized over J (BLAS matrix multiply).
#
# @param betahat     J x R matrix of effect estimates
# @param V_scalar    Scalar prior variance multiplier (>= 0)
# @param eigen_cache Precomputed cache from precompute_eigen_cache()
#
# @return J x (K+1) matrix of log-likelihoods
loglik_precomputed <- function(betahat, V_scalar, eigen_cache) {
  if (!eigen_cache$is_common_cov) {
    # C++ fast path for non-common-cov (eliminates K*J R-level loops)
    return(loglik_non_common_cpp(betahat, V_scalar,
                                  eigen_cache$log_det_svs,
                                  eigen_cache$components))
  }

  # Common-cov path: already vectorized over J via BLAS
  J <- nrow(betahat)
  R <- ncol(betahat)
  K <- length(eigen_cache$components)
  llik <- matrix(0, J, K + 1)
  const <- -R / 2 * log(2 * pi)
  log_det_svs <- eigen_cache$log_det_svs  # scalar

  # Null component: N(0, SVS)
  Q_null <- eigen_cache$components[[1]]$Q
  B_null <- betahat %*% Q_null                   # J x R (BLAS)
  mahal_null <- rowSums(B_null^2)                # J-vector
  llik[, 1] <- const - 0.5 * log_det_svs - 0.5 * mahal_null

  # Non-null components
  for (k in seq_len(K)) {
    comp <- eigen_cache$components[[k]]
    d_k <- comp$eigenvalues
    B <- betahat %*% comp$Q                     # J x R (BLAS)
    log_det <- log_det_svs + sum(log1p(V_scalar * d_k))
    inv_factors <- 1 / (1 + V_scalar * d_k)     # R-vector
    mahal <- drop(B^2 %*% inv_factors)           # J-vector
    llik[, k + 1] <- const - 0.5 * log_det - 0.5 * mahal
  }
  llik
}

# R fallback for loglik_precomputed (for testing)
loglik_precomputed_R <- function(betahat, V_scalar, eigen_cache) {
  J <- nrow(betahat)
  R <- ncol(betahat)
  K <- length(eigen_cache$components)
  llik <- matrix(0, J, K + 1)
  const <- -R / 2 * log(2 * pi)

  if (eigen_cache$is_common_cov) {
    log_det_svs <- eigen_cache$log_det_svs
    Q_null <- eigen_cache$components[[1]]$Q
    B_null <- betahat %*% Q_null
    mahal_null <- rowSums(B_null^2)
    llik[, 1] <- const - 0.5 * log_det_svs - 0.5 * mahal_null
    for (k in seq_len(K)) {
      comp <- eigen_cache$components[[k]]
      d_k <- comp$eigenvalues
      B <- betahat %*% comp$Q
      log_det <- log_det_svs + sum(log1p(V_scalar * d_k))
      inv_factors <- 1 / (1 + V_scalar * d_k)
      mahal <- drop(B^2 %*% inv_factors)
      llik[, k + 1] <- const - 0.5 * log_det - 0.5 * mahal
    }
  } else {
    log_det_svs_vec <- eigen_cache$log_det_svs
    for (j in seq_len(J)) {
      Q_j <- eigen_cache$components[[1]]$Q[, , j]
      b_rot <- crossprod(Q_j, betahat[j, ])
      llik[j, 1] <- const - 0.5 * log_det_svs_vec[j] - 0.5 * sum(b_rot^2)
    }
    for (k in seq_len(K)) {
      comp <- eigen_cache$components[[k]]
      for (j in seq_len(J)) {
        d_k <- comp$eigenvalues[, j]
        b_rot <- crossprod(comp$Q[, , j], betahat[j, ])
        log_det <- log_det_svs_vec[j] + sum(log1p(V_scalar * d_k))
        inv_factors <- 1 / (1 + V_scalar * d_k)
        mahal <- sum(b_rot^2 * inv_factors)
        llik[j, k + 1] <- const - 0.5 * log_det - 0.5 * mahal
      }
    }
  }
  llik
}

# Compute posterior moments using precomputed eigendecomposition.
#
# For each variable j and mixture component k, the posterior under
# N(bhat | 0, SVS + V*U_k) with prior N(0, V*U_k) is:
#   mu_post  = G diag(V*d/(1+V*d)) Q' bhat   (R-vector)
#   Cov_post = G diag(V*d/(1+V*d)) G'         (R x R, same for all j in common_cov)
#
# The mixture-weighted posterior mean is sum_k w_jk * mu_post_jk.
# The posterior second moment mu2 = E[bb'] = sum_k w_jk (C_k + m_jk m_jk').
#
# @param betahat     J x R matrix of effect estimates
# @param V_scalar    Scalar prior variance multiplier
# @param eigen_cache Precomputed cache
# @param pi_V_post   J x (K+1) matrix of posterior mixture weights
# @param em_var_wt   (K+1) x J matrix of P(j|k) posterior variable weights
#   for EM update (from compute_variable_posterior_weights). NULL if EM
#   not needed.
#
# @param reduce_params If non-NULL, a list(alpha, d, v_inv) for computing
#   reduced statistics (bxxb, vbxxb, alpha_mu2_sum, mu2_diag) instead of
#   the full J x R x R post_mean2 array.  Saves O(J*R^2) memory.
# @return list with post_mean, post_neg, post_zero, prior_scale_em_update,
#   and either post_mean2 (when reduce_params is NULL) or reduced statistics.
posterior_precomputed <- function(betahat, V_scalar, eigen_cache, pi_V_post,
                                  em_var_wt = NULL, reduce_params = NULL) {
  if (!eigen_cache$is_common_cov) {
    # C++ fast path for non-common-cov (eliminates K*J R-level loops)
    em_wt <- if (!is.null(em_var_wt)) em_var_wt else matrix(0, 0, 0)
    result <- posterior_non_common_cpp(betahat, V_scalar,
                                       eigen_cache$components,
                                       pi_V_post, em_wt)
    # Reduce if requested (post_mean2 was allocated in C++; compute stats, free it)
    if (!is.null(reduce_params)) {
      result <- reduce_post_mean2(result, reduce_params)
    }
    return(result)
  }

  # Common-cov path: BLAS-vectorized over J
  J <- nrow(betahat)
  R <- ncol(betahat)
  K <- length(eigen_cache$components)
  do_reduce <- !is.null(reduce_params)

  post_mean  <- matrix(0, J, R)
  post_neg   <- matrix(0, J, R)
  post_zero  <- matrix(0, J, R)
  em_update  <- numeric(K + 1)

  if (do_reduce) {
    alpha <- reduce_params$alpha
    d_var <- reduce_params$d
    bxxb <- matrix(0, R, R)
    alpha_mu2_sum <- matrix(0, R, R)
    mu2_diag <- matrix(0, J, R)   # diagonal of second moment: E[b_{j,r}^2]
  } else {
    post_mean2 <- array(0, c(J, R, R))
  }

  # Null component: beta = 0 exactly
  w_null <- pi_V_post[, 1]
  post_zero <- matrix(w_null, J, R)

  for (k in seq_len(K)) {
    comp <- eigen_cache$components[[k]]
    d_k <- comp$eigenvalues
    Q_k <- comp$Q
    G_k <- comp$G
    w_k <- pi_V_post[, k + 1]  # J-vector

    # Shrinkage factors
    shrink <- V_scalar * d_k / (1 + V_scalar * d_k)      # R-vector
    inv_factor <- 1 / (1 + V_scalar * d_k)                # R-vector

    # Posterior covariance (SAME for all j in common_cov)
    G_scaled <- G_k * rep(shrink, each = R)  # scale columns of G
    C_k <- G_scaled %*% t(G_k)
    C_k <- (C_k + t(C_k)) / 2

    # Posterior means for all J variables (BLAS)
    BQ <- betahat %*% Q_k                              # J x R
    BQ_shrunk <- t(t(BQ) * shrink)                      # scale columns
    M_k <- BQ_shrunk %*% t(G_k)                         # J x R

    # Accumulate posterior mean
    post_mean <- post_mean + w_k * M_k

    if (do_reduce) {
      # Accumulate bxxb, alpha_mu2_sum, mu2_diag using BLAS (no J x R x R)
      aw_k  <- alpha * w_k         # J-vector
      daw_k <- d_var * aw_k        # J-vector

      # C_k contribution (R x R, same for all j)
      bxxb          <- bxxb          + sum(daw_k) * C_k
      alpha_mu2_sum <- alpha_mu2_sum + sum(aw_k) * C_k

      # Outer product contributions via BLAS crossprod: sum_j wt_j * M_kj M_kj'
      sqrt_daw <- sqrt(daw_k)
      bxxb <- bxxb + crossprod(sqrt_daw * M_k)           # R x R
      sqrt_aw <- sqrt(aw_k)
      alpha_mu2_sum <- alpha_mu2_sum + crossprod(sqrt_aw * M_k)

      # mu2_diag: E[b_{j,r}^2] contribution from component k
      # mu2_diag[j,r] += w_k[j] * (C_k[r,r] + M_k[j,r]^2)
      diag_Ck <- diag(C_k)
      mu2_diag <- mu2_diag + w_k * matrix(diag_Ck, J, R, byrow = TRUE) +
                  w_k * M_k^2
    } else {
      # Build full J x R x R post_mean2 via C++
      post_mean2 <- accumulate_post_mean2_common_cpp(post_mean2, M_k, C_k, w_k)
    }

    # Sign probabilities for LFSR
    diag_Ck <- diag(C_k)
    diag_Ck[diag_Ck < sqrt(.Machine$double.eps) * max(diag_Ck)] <- 0
    sd_k <- sqrt(diag_Ck)
    for (r in seq_len(R)) {
      if (sd_k[r] > 0) {
        post_neg[, r] <- post_neg[, r] +
          w_k * pnorm(0, M_k[, r], sd_k[r])
      } else {
        post_zero[, r] <- post_zero[, r] + w_k
      }
    }

    # EM statistic
    em_wt <- if (!is.null(em_var_wt)) em_var_wt[k + 1, ] else w_k
    tr_term <- V_scalar * sum(inv_factor)
    em_per_var <- V_scalar^2 * drop(BQ^2 %*% (d_k * inv_factor^2))
    em_update[k + 1] <- sum(em_wt * (tr_term + em_per_var))
  }

  result <- list(
    post_mean  = post_mean,
    post_neg   = post_neg,
    post_zero  = post_zero,
    prior_scale_em_update = em_update
  )

  if (do_reduce) {
    # vbxxb = tr(v_inv * bxxb) for common-cov (svs_inv[j] = d[j]*v_inv)
    v_inv <- reduce_params$v_inv
    result$bxxb          <- bxxb
    result$vbxxb         <- sum(v_inv * bxxb)
    result$alpha_mu2_sum <- alpha_mu2_sum
    result$mu2_diag      <- mu2_diag
  } else {
    result$post_mean2 <- post_mean2
  }

  return(result)
}

# Reduce a full J x R x R post_mean2 to summary statistics.
# Used for non-common-cov path where C++ builds the full array.
reduce_post_mean2 <- function(post_result, reduce_params) {
  alpha   <- reduce_params$alpha
  d_var   <- reduce_params$d
  svs_inv <- reduce_params$svs_inv  # list of J (or 1 for common-cov) R x R matrices
  pm2     <- post_result$post_mean2
  J       <- length(alpha)
  R       <- ncol(post_result$post_mean)

  bxxb          <- matrix(0, R, R)
  alpha_mu2_sum <- matrix(0, R, R)
  vbxxb         <- 0
  mu2_diag      <- matrix(0, J, R)

  for (j in seq_len(J)) {
    mu2_j <- pm2[j, , ]
    if (!is.matrix(mu2_j)) dim(mu2_j) <- c(R, R)
    a_j <- alpha[j]
    bxxb          <- bxxb + d_var[j] * a_j * mu2_j
    alpha_mu2_sum <- alpha_mu2_sum + a_j * mu2_j
    vbxxb         <- vbxxb + a_j * sum(svs_inv[[min(j, length(svs_inv))]] * mu2_j)
    mu2_diag[j, ] <- diag(mu2_j)
  }

  post_result$post_mean2    <- NULL   # free J x R x R
  post_result$bxxb          <- bxxb
  post_result$vbxxb         <- vbxxb
  post_result$alpha_mu2_sum <- alpha_mu2_sum
  post_result$mu2_diag      <- mu2_diag
  return(post_result)
}

# R fallback for posterior_precomputed (for testing)
posterior_precomputed_R <- function(betahat, V_scalar, eigen_cache, pi_V_post,
                                    em_var_wt = NULL) {
  J <- nrow(betahat)
  R <- ncol(betahat)
  K <- length(eigen_cache$components)

  post_mean  <- matrix(0, J, R)
  post_mean2 <- array(0, c(J, R, R))
  post_neg   <- matrix(0, J, R)
  post_zero  <- matrix(0, J, R)
  em_update  <- numeric(K + 1)

  w_null <- pi_V_post[, 1]
  post_zero <- matrix(w_null, J, R)

  if (eigen_cache$is_common_cov) {
    for (k in seq_len(K)) {
      comp <- eigen_cache$components[[k]]
      d_k <- comp$eigenvalues
      Q_k <- comp$Q
      G_k <- comp$G
      w_k <- pi_V_post[, k + 1]
      shrink <- V_scalar * d_k / (1 + V_scalar * d_k)
      inv_factor <- 1 / (1 + V_scalar * d_k)
      G_scaled <- G_k * rep(shrink, each = R)
      C_k <- G_scaled %*% t(G_k)
      C_k <- (C_k + t(C_k)) / 2
      BQ <- betahat %*% Q_k
      BQ_shrunk <- t(t(BQ) * shrink)
      M_k <- BQ_shrunk %*% t(G_k)
      post_mean <- post_mean + w_k * M_k
      for (j in seq_len(J)) {
        post_mean2[j, , ] <- post_mean2[j, , ] +
          w_k[j] * (C_k + tcrossprod(M_k[j, ]))
      }
      diag_Ck <- diag(C_k)
      diag_Ck[diag_Ck < sqrt(.Machine$double.eps) * max(diag_Ck)] <- 0
      sd_k <- sqrt(diag_Ck)
      for (r in seq_len(R)) {
        if (sd_k[r] > 0) {
          post_neg[, r] <- post_neg[, r] + w_k * pnorm(0, M_k[, r], sd_k[r])
        } else {
          post_zero[, r] <- post_zero[, r] + w_k
        }
      }
      em_wt <- if (!is.null(em_var_wt)) em_var_wt[k + 1, ] else w_k
      tr_term <- V_scalar * sum(inv_factor)
      em_per_var <- V_scalar^2 * drop(BQ^2 %*% (d_k * inv_factor^2))
      em_update[k + 1] <- sum(em_wt * (tr_term + em_per_var))
    }
  } else {
    for (k in seq_len(K)) {
      comp <- eigen_cache$components[[k]]
      w_k <- pi_V_post[, k + 1]
      for (j in seq_len(J)) {
        Q_j <- comp$Q[, , j]
        G_j <- comp$G[, , j]
        d_k <- comp$eigenvalues[, j]
        shrink <- V_scalar * d_k / (1 + V_scalar * d_k)
        inv_factor <- 1 / (1 + V_scalar * d_k)
        G_scaled <- G_j * rep(shrink, each = R)
        C_k <- G_scaled %*% t(G_j)
        C_k <- (C_k + t(C_k)) / 2
        b_rot <- crossprod(Q_j, betahat[j, ])
        m_j <- drop(G_j %*% (shrink * b_rot))
        post_mean[j, ] <- post_mean[j, ] + w_k[j] * m_j
        post_mean2[j, , ] <- post_mean2[j, , ] +
          w_k[j] * (C_k + tcrossprod(m_j))
        diag_Ck <- diag(C_k)
        diag_Ck[diag_Ck < sqrt(.Machine$double.eps) * max(diag_Ck)] <- 0
        sd_k <- sqrt(diag_Ck)
        for (r in seq_len(R)) {
          if (sd_k[r] > 0)
            post_neg[j, r] <- post_neg[j, r] +
              w_k[j] * pnorm(0, m_j[r], sd_k[r])
          else
            post_zero[j, r] <- post_zero[j, r] + w_k[j]
        }
        tr_term <- V_scalar * sum(inv_factor)
        em_j <- V_scalar^2 * sum(d_k * inv_factor^2 * b_rot^2)
        em_wt_j <- if (!is.null(em_var_wt)) em_var_wt[k + 1, j] else w_k[j]
        em_update[k + 1] <- em_update[k + 1] + em_wt_j * (tr_term + em_j)
      }
    }
  }

  list(
    post_mean  = post_mean,
    post_mean2 = post_mean2,
    post_neg   = post_neg,
    post_zero  = post_zero,
    prior_scale_em_update = em_update
  )
}