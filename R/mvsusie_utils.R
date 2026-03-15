# Utility functions for mvsusieR.
#
# Includes matrix operations, numerical helpers, and lfsr computation functions.

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
# Set verbose = FALSE to suppress the summary messages.
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

# Update mixture prior weights and prune near-zero components.
#
# Shared helper called by both update_model_variance.mv_individual and
# update_model_variance.mv_ss. Prunes every 10 iterations (10, 20, 30, ...)
# to remove mixture components whose weight has collapsed below threshold.
update_mixture_weights_and_prune <- function(model, params) {
  if (!isTRUE(params$estimate_prior_mixture_weights)) return(model)
  K <- length(model$pi_V)
  if (K <= 1 || !any(!sapply(model$pi_V_posterior, is.null))) return(model)

  model <- update_mixture_weights(model,
             method = params$mixture_weight_method %||% "mixsqp",
             update_null = (model$null_weight > 0))
  iter <- model$ibss_iter %||% 0
  if (iter > 0 && iter %% 10 == 0) {
    model <- prune_mixture_components(model, threshold = 1e-8)
  }
  model$ibss_iter <- iter + 1
  return(model)
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
    R <- nrow(SVS)

    # Cholesky of SVS: compute once, reuse for all K components
    L_upper <- safe_chol(SVS)
    L <- t(L_upper)
    log_det_svs <- chol2ldet(L_upper)
    L_inv <- forwardsolve(L, diag(R))
    L_inv_t <- t(L_inv)

    components <- vector("list", K)
    for (k in seq_len(K)) {
      # Only the U_k-dependent part: M = L_inv U_k L_inv', eigen, Q, G
      M <- L_inv %*% V_structure[[k]] %*% L_inv_t
      M <- (M + t(M)) / 2
      eig <- eigen(M, symmetric = TRUE)
      P <- eig$vectors
      d <- pmax(eig$values, 0)
      components[[k]] <- list(
        Q = L_inv_t %*% P,
        G = L %*% P,
        eigenvalues = d,
        log_det_svs = log_det_svs
      )
    }
    list(
      is_common_cov = TRUE,
      log_det_svs   = log_det_svs,
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
    return(precompute_eigen_cache_non_common_rcpp(svs, V_structure))
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
loglik_precomputed <- function(betahat, V_scalar, eigen_cache,
                                BQ_cache = NULL) {
  if (!eigen_cache$is_common_cov) {
    return(loglik_non_common_rcpp(betahat, V_scalar,
                                  eigen_cache$log_det_svs,
                                  eigen_cache$components))
  }

  # Common-cov C++ fast path
  BQ_list <- if (!is.null(BQ_cache)) BQ_cache else list()
  res <- loglik_common_rcpp(betahat, V_scalar,
                             eigen_cache$log_det_svs,
                             eigen_cache$components,
                             BQ_list)
  res$llik
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
# @param reduce_params If non-NULL, a list(alpha, d, svs_inv, v_inv,
#   gram_xtx) for computing reduced statistics (bxxb, vbxxb,
#   alpha_mu2_sum, mu2_diag) instead of the full J x R x R mu2 array.
#   When gram_xtx is non-NULL, bxxb uses the per-outcome Gram matrix
#   (3D missing-data path).
# @return list with mu, post_neg, post_zero, prior_scale_em_update,
#   and either mu2 (when reduce_params is NULL) or reduced statistics.
posterior_precomputed <- function(betahat, V_scalar, eigen_cache, pi_V_post,
                                  em_var_wt = NULL, reduce_params = NULL,
                                  BQ_cache = NULL) {
  if (!eigen_cache$is_common_cov) {
    # C++ fast path for non-common-cov (eliminates K*J R-level loops)
    em_wt <- if (!is.null(em_var_wt)) em_var_wt else matrix(0, 0, 0)
    result <- posterior_non_common_rcpp(betahat, V_scalar,
                                       eigen_cache$components,
                                       pi_V_post, em_wt)
    # Reduce if requested (mu2 was allocated in C++; compute stats, free it)
    if (!is.null(reduce_params)) {
      result <- reduce_mu2(result, reduce_params)
    }
    return(result)
  }

  # Common-cov C++ fast path
  em_wt_mat <- if (!is.null(em_var_wt)) em_var_wt else matrix(0, 0, 0)
  BQ_list <- if (!is.null(BQ_cache)) BQ_cache else list()
  do_reduce <- !is.null(reduce_params)

  if (do_reduce && !is.null(reduce_params$gram_xtx)) {
    # When gram_xtx is present, the C++ inline reduction uses scalar d[j]
    # which is wrong for the 3D path. Get full mu2 from C++ and reduce
    # in R with the correct G-weighted bxxb.
    result <- posterior_common_rcpp(betahat, V_scalar,
                                    eigen_cache$components,
                                    pi_V_post, em_wt_mat, BQ_list,
                                    FALSE, numeric(0),
                                    numeric(0), matrix(0, 0, 0))
    reduce_mu2(result, reduce_params)
  } else if (do_reduce) {
    posterior_common_rcpp(betahat, V_scalar,
                          eigen_cache$components,
                          pi_V_post, em_wt_mat, BQ_list,
                          do_reduce, reduce_params$alpha,
                          reduce_params$d, reduce_params$v_inv)
  } else {
    posterior_common_rcpp(betahat, V_scalar,
                          eigen_cache$components,
                          pi_V_post, em_wt_mat, BQ_list,
                          do_reduce, numeric(0),
                          numeric(0), matrix(0, 0, 0))
  }
}

# Reduce a full J x R x R mu2 array to summary statistics.
# Used for non-common-cov path where C++ builds the full array.
#
# When gram_xtx (J x R x R) is provided via reduce_params, bxxb is computed
# using the per-outcome Gram matrix (Hadamard product) instead of scalar d[j].
# This gives the correct bxxb for the 3D missing-data path.
reduce_mu2 <- function(post_result, reduce_params) {
  alpha    <- reduce_params$alpha
  d_var    <- reduce_params$d
  svs_inv  <- reduce_params$svs_inv  # list of J (or 1) R x R matrices
  gram_xtx <- reduce_params$gram_xtx # J x R x R or NULL
  pm2      <- post_result$mu2
  J        <- length(alpha)
  R        <- ncol(post_result$mu)

  bxxb          <- matrix(0, R, R)
  alpha_mu2_sum <- matrix(0, R, R)
  vbxxb         <- 0
  mu2_diag      <- matrix(0, J, R)

  for (j in seq_len(J)) {
    mu2_j <- pm2[j, , ]
    if (!is.matrix(mu2_j)) dim(mu2_j) <- c(R, R)
    a_j <- alpha[j]
    if (!is.null(gram_xtx)) {
      G_j <- gram_xtx[j, , ]
      if (!is.matrix(G_j)) dim(G_j) <- c(R, R)
      bxxb <- bxxb + a_j * (G_j * mu2_j)          # Hadamard with Gram
    } else {
      bxxb <- bxxb + d_var[j] * a_j * mu2_j        # scalar d weighting
    }
    alpha_mu2_sum <- alpha_mu2_sum + a_j * mu2_j
    vbxxb         <- vbxxb + a_j * sum(svs_inv[[min(j, length(svs_inv))]] * mu2_j)
    mu2_diag[j, ] <- diag(mu2_j)
  }

  post_result$mu2           <- NULL   # free J x R x R
  post_result$bxxb          <- bxxb
  post_result$vbxxb         <- vbxxb
  post_result$alpha_mu2_sum <- alpha_mu2_sum
  post_result$mu2_diag      <- mu2_diag
  return(post_result)
}

