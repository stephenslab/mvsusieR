# Pure R reference implementations for functions that have C++ counterparts.
#
# These are kept for testing: unit tests compare C++ output against these
# R implementations to verify correctness. They are NOT used in production
# code paths.
#
# Functions:
#   precompute_eigen_cache_R   - R reference for precompute_eigen_cache
#   loglik_precomputed_R       - R reference for loglik_precomputed
#   posterior_precomputed_R    - R reference for posterior_precomputed
#   compute_XtR_3d_R           - R reference for compute_XtR_3d
#   compute_Xb_3d_R            - R reference for compute_Xb_3d
#   compute_VinvR_3d_R         - R reference for compute_VinvR_3d
#   compute_betahat_3d_R       - R reference for compute_betahat_3d

# =============================================================================
# Eigendecomposition precomputation (R reference)
# =============================================================================

precompute_eigen_cache_R <- function(svs, V_structure, is_common_cov) {
  K <- length(V_structure)

  if (is_common_cov) {
    SVS <- svs[[1]]
    R <- nrow(SVS)

    L_upper <- safe_chol(SVS)
    L <- t(L_upper)
    log_det_svs <- chol2ldet(L_upper)
    L_inv <- forwardsolve(L, diag(R))
    L_inv_t <- t(L_inv)

    components <- vector("list", K)
    for (k in seq_len(K)) {
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

# =============================================================================
# Log-likelihood (R reference)
# =============================================================================

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

# =============================================================================
# Posterior moments (R reference)
# =============================================================================

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

# =============================================================================
# Missing data 3D operations (R references)
# =============================================================================

compute_XtR_3d_R <- function(data, R_mat) {
  my <- data$miss3d
  Vinv <- my$Vinv
  N <- nrow(R_mat)
  R_dim <- data$R
  J <- data$p
  pattern_assign <- my$pattern_assign

  VinvR <- matrix(0, N, R_dim)
  for (k in seq_along(Vinv)) {
    idx <- which(pattern_assign == k)
    if (length(idx) == 0) next
    VinvR[idx, ] <- R_mat[idx, , drop = FALSE] %*% t(Vinv[[k]])
  }

  if (my$method == "approximate") {
    res <- t(sapply(seq_len(J), function(j) {
      colSums(my$X_3d[, j, ] * VinvR)
    }))
  } else {
    colsum_VinvR <- colSums(VinvR)
    res <- t(sapply(seq_len(J), function(j) {
      colSums(my$X_3d[, j, ] * VinvR) -
        as.numeric(crossprod(my$Xbar[j, , ], colsum_VinvR))
    }))
  }

  if (R_dim == 1) res <- t(res)
  if (nrow(res) != J) res <- t(res)
  return(res)
}

compute_Xb_3d_R <- function(data, b) {
  my <- data$miss3d
  if (is.vector(b)) b <- matrix(b, length(b), 1)
  R_dim <- data$R
  N <- dim(my$X_3d)[1]

  Xb <- sapply(seq_len(R_dim), function(r) my$X_3d[, , r] %*% b[, r])
  if (!is.matrix(Xb)) Xb <- matrix(Xb, N, R_dim)

  if (my$method == "exact" && !is.null(my$Xbar)) {
    J <- nrow(b)
    Xbarb <- Reduce("+", lapply(seq_len(J), function(j) {
      Xbar_j <- my$Xbar[j, , ]
      if (!is.matrix(Xbar_j)) Xbar_j <- matrix(Xbar_j, R_dim, R_dim)
      Xbar_j %*% b[j, ]
    }))
    Xb <- Xb - matrix(as.numeric(Xbarb), N, R_dim, byrow = TRUE)
  }

  if (nrow(Xb) != N) Xb <- t(Xb)
  return(Xb)
}

compute_VinvR_3d_R <- function(data, mat) {
  my <- data$miss3d
  Vinv <- my$Vinv
  N <- nrow(mat)
  R_dim <- ncol(mat)
  pattern_assign <- my$pattern_assign

  result <- matrix(0, N, R_dim)
  for (k in seq_along(Vinv)) {
    idx <- which(pattern_assign == k)
    if (length(idx) == 0) next
    result[idx, ] <- mat[idx, , drop = FALSE] %*% t(Vinv[[k]])
  }
  return(result)
}

compute_betahat_3d_R <- function(data, XtR) {
  J <- data$p
  R_dim <- data$R
  svs <- data$svs

  bhat <- t(sapply(seq_len(J), function(j) svs[[j]] %*% XtR[j, ]))
  bhat[is.nan(bhat)] <- 0
  if (R_dim == 1) bhat <- t(bhat)
  return(bhat)
}
