# S3 method implementations for multivariate individual-level data.
#
# Methods for the mv_individual class, providing implementations of
# susieR's internal generics for dense (X, Y) data.

#' @importFrom susieR susie_get_cs
#' @importFrom ashr compute_lfsr
NULL

# =============================================================================
# DATA CONFIGURATION
# =============================================================================

#' @keywords internal
get_var_y.mv_individual <- function(data, ...) {
  if (data$R > 1) cov(data$Y) else var(data$Y)
}

# =============================================================================
# IBSS INITIALIZATION (overrides susieR's ibss_initialize for MV data)
#
# susieR's default ibss_initialize uses validate_init/prune_single_effects
# which assume univariate format (2D mu/mu2 matrices). For multivariate
# models (3D mu array, list-of-array mu2), we provide MV-specific init
# that handles model_init (warm-starting from a previous fit).
# =============================================================================

#' @keywords internal
ibss_initialize.mv_individual <- function(data, params) {
  var_y <- get_var_y(data)

  if (data$p < params$L) {
    params$L <- data$p
  }

  # Residual variance for MV is always a matrix -- skip scalar checks
  if (is.null(params$residual_variance)) {
    params$residual_variance <- var_y
  }

  # Handle model initialization (warm-starting)
  if (!is.null(params$model_init)) {
    m_init <- params$model_init
    if (!inherits(m_init, "susie") && !inherits(m_init, "mvsusie")) {
      stop("model_init must be a 'susie' or 'mvsusie' object")
    }

    L <- params$L
    J <- data$p
    R <- data$R

    # Validate alpha dimensions
    if (!is.matrix(m_init$alpha) || ncol(m_init$alpha) != J) {
      stop("model_init$alpha must be an L x J matrix matching data dimensions")
    }

    # Adjust L based on init
    L_init <- nrow(m_init$alpha)
    if (L_init > L) {
      warning(paste0(
        "Requested L = ", L, " is smaller than the ", L_init,
        " effects in model_init; using L = ", L_init, " instead."
      ))
      params$L <- L_init
      L <- L_init
    }

    # Create fresh model with correct dimensions via S3 dispatch
    mat_init <- initialize_susie_model(data, params, var_y)

    # Overwrite alpha from model_init (pad with uniform if L > L_init)
    mat_init$alpha[seq_len(L_init), ] <- m_init$alpha

    # Overwrite mu from model_init
    # m_init$mu is L x J (R=1) or L x J x R (R>1)
    if (R == 1) {
      # mu from output is L x J matrix; model$mu is L x J x 1 array
      mu_init <- m_init$mu
      if (is.matrix(mu_init)) {
        for (l in seq_len(L_init)) {
          mat_init$mu[l, , 1] <- mu_init[l, ]
        }
      }
    } else {
      # mu from output is L x J x R array
      if (length(dim(m_init$mu)) == 3) {
        for (l in seq_len(L_init)) {
          mat_init$mu[l, , ] <- m_init$mu[l, , ]
        }
      }
    }

    # Overwrite V (scalar per effect) from model_init
    if (!is.null(m_init$V) && length(m_init$V) >= L_init) {
      mat_init$V[seq_len(L_init)] <- m_init$V[seq_len(L_init)]
    }

    # Overwrite V_structure if available
    if (!is.null(m_init$V_structure)) {
      mat_init$V_structure <- m_init$V_structure
    }

    # Reset iteration-specific values (will be recomputed)
    mat_init$KL  <- rep(as.numeric(NA), L)
    mat_init$lbf <- rep(as.numeric(NA), L)
  } else {
    mat_init <- initialize_susie_model(data, params, var_y)
  }

  # Initialize fitted values and null index
  fitted     <- initialize_fitted(data, mat_init)
  null_index <- 0

  # Preserve model class
  model_class <- class(mat_init)

  model <- c(mat_init, list(null_index = null_index), fitted)

  if (inherits(mat_init, "susie")) {
    class(model) <- model_class
  } else {
    class(model) <- "susie"
  }
  model$converged <- FALSE

  # Trigger eigendecomposition precomputation if requested
  verbose <- isTRUE(params$verbose)
  if (isTRUE(params$precompute_eigendecomp)) {
    svs <- if (!is.null(model$svs)) model$svs else data$svs
    model$eigen_cache <- precompute_eigen_cache(
      svs, model$V_structure, data$is_common_cov)
    if (verbose) {
      message("Eigendecomposition cache: K=", length(model$V_structure),
              ", common_cov=", data$is_common_cov,
              " [mem: ", sprintf("%.2f", mem_used_gb()), " GB]")
    }
  }

  if (verbose) {
    K <- length(model$V_structure)
    message(sprintf("Model initialized: J=%d, R=%d, L=%d, K=%d [mem: %.2f GB]",
                    data$p, data$R, params$L, K, mem_used_gb()))
  }

  # Initial imputation for R>1 missing data (variational EM E-step).
  # At this point Xr is all zeros, so imputed values = conditional mean
  # given mu=0 (i.e., mean imputation from the prior).
  if (data$any_missing && data$R > 1 && !is.null(data$impute_info)) {
    v_inv <- if (!is.null(model$residual_variance_inv))
      model$residual_variance_inv else data$residual_variance_inv
    model$pattern_cache <- precompute_pattern_cache(v_inv,
                                                     data$impute_info$patterns)
    imp <- impute_missing_Y(data$Y, model$Xr, v_inv,
                            data$impute_info, model$pattern_cache)
    model$Y_imputed <- imp$Y
    model$Y_cov <- imp$Y_cov
    model$sum_neg_ent_Y_miss <- imp$sum_neg_ent_Y_miss
  }

  return(model)
}

# =============================================================================
# RESIDUALS
# =============================================================================

#' @keywords internal
compute_residuals.mv_individual <- function(data, params, model, l, ...) {
  # Impute missing Y at start of each IBSS iteration (first SER effect).
  # This is the variational EM E-step: fill in missing entries using
  # E[Y_miss | Y_obs] = mu_miss - Lambda_{MM}^-1 Lambda_{MO} (Y_obs - mu_obs)
  # where Lambda = V^-1 (precision) and mu = model$Xr (current fitted values).
  # Skip imputation when using per-outcome (3d) missing data methods.
  if (l == 1 && data$any_missing && data$R > 1 &&
      !is.null(data$impute_info) && is.null(data$miss3d)) {
    v_inv <- if (!is.null(model$residual_variance_inv))
      model$residual_variance_inv else data$residual_variance_inv
    # Recompute pattern cache (Vinv may have changed after V update)
    model$pattern_cache <- precompute_pattern_cache(v_inv,
                                                     data$impute_info$patterns)
    imp <- impute_missing_Y(data$Y, model$Xr, v_inv,
                            data$impute_info, model$pattern_cache)
    model$Y_imputed <- imp$Y
    model$Y_cov <- imp$Y_cov
    model$sum_neg_ent_Y_miss <- imp$sum_neg_ent_Y_miss
    # Re-center imputed Y (column means shift after imputation)
    # and track adjustment for intercept recovery
    col_adj <- colMeans(model$Y_imputed)
    model$Y_imputed <- t(t(model$Y_imputed) - col_adj)
    model$Y_mean_adjustment <- col_adj
  }

  # Use imputed Y for R>1 missing data, otherwise original Y
  Y <- if (!is.null(model$Y_imputed)) model$Y_imputed else data$Y

  # Fitted for effect l: X %*% (alpha_l * mu_l)
  b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
  Xb_l <- compute_Xb(data, b_l)  # N x R

  # Residual: Y - (Xr - Xb_l) = Y - Xr + Xb_l
  R_mat <- Y - model$Xr + Xb_l  # N x R

  # Zero out missing entries only for R=1 (complete-case approach) or
  # when using per-outcome (3d) methods (missing entries should not contribute).
  # For R>1 with imputation, imputed entries participate fully.
  if (data$any_missing &&
      (is.null(model$Y_imputed) || !is.null(data$miss3d))) {
    R_mat[data$Y_na] <- 0
  }

  # Inject updated Vinv into local data copy (R copy-on-modify) so that
  # compute_XtR_3d sees the current precision matrices after est_rv update.
  if (!is.null(model$miss3d_Vinv))
    data$miss3d$Vinv <- model$miss3d_Vinv

  # X'R for SER computation
  XtR <- compute_XtR(data, R_mat)  # J x R

  model$residuals        <- XtR
  model$fitted_without_l <- model$Xr - Xb_l
  model$raw_residuals    <- R_mat
  return(model)
}

# =============================================================================
# SER STATISTICS
# =============================================================================

#' @keywords internal
compute_ser_statistics.mv_individual <- function(data, params, model, l, ...) {
  # Inject updated svs into local data copy so that compute_betahat_3d
  # sees the current GLS covariances after residual variance update.
  if (!is.null(model$svs) && !is.null(data$miss3d))
    data$svs <- model$svs

  # betahat_j: GLS estimate (svs %*% XtR for 3d, or XtR/d for standard)
  betahat <- compute_betahat(data, model$residuals)

  # Use model's svs/svs_inv if available (updated after residual variance change),
  # otherwise fall back to data's precomputed values
  svs     <- if (!is.null(model$svs)) model$svs else data$svs
  svs_inv <- if (!is.null(model$svs_inv)) model$svs_inv else data$svs_inv

  # Moment-based warm start for optim: multivariate analogue of susieR's
  # log(max(betahat^2 - shat2, 1)). For each variable j, compute
  # tr(betahat_j betahat_j') - tr(S_j) as a scalar measure of signal
  # strength, then take the max across variables.
  signal <- sapply(seq_len(nrow(betahat)), function(j) {
    sum(betahat[j, ]^2) - sum(diag(svs[[min(j, length(svs))]]))
  })
  optim_init <- log(max(c(max(signal, na.rm = TRUE), 1)))

  # Precompute BQ_cache for Brent optimizer: betahat %*% Q_k for all K.
  # Only when using optim (Brent) with eigendecomposition precomputation,
  # since these V-independent products are recomputed ~17 times otherwise.
  BQ_cache <- NULL
  if (!is.null(params$estimate_prior_method) &&
      params$estimate_prior_method == "optim" &&
      !is.null(model$eigen_cache) &&
      model$eigen_cache$is_common_cov) {
    K <- length(model$eigen_cache$components)
    BQ_cache <- vector("list", K)
    for (k in seq_len(K))
      BQ_cache[[k]] <- betahat %*% model$eigen_cache$components[[k]]$Q
  }

  list(
    betahat      = betahat,
    shat2        = svs,            # list of J RxR matrices
    shat2_inv    = svs_inv,        # list of J RxR matrices
    optim_init   = optim_init,
    optim_bounds = c(-30, 15),
    optim_scale  = "log",
    BQ_cache     = BQ_cache
  )
}

# =============================================================================
# LOG BAYES FACTOR & POSTERIOR WEIGHTS
# =============================================================================

#' @keywords internal
loglik.mv_individual <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  R <- data$R
  J <- data$p
  n_thread <- if (!is.null(params$n_thread)) params$n_thread else 1L

  # When V ~= 0 or NA, the prior is null: lbf = 0 everywhere, uniform alpha
  if (is.na(V) || V < params$prior_tol) {
    if (!is.null(l)) {
      model$alpha[l, ]        <- 1 / J
      model$lbf[l]            <- 0
      model$lbf_variable[l, ] <- rep(0, J)
      model$pi_V_posterior[l] <- list(NULL)
      model$llik_cache        <- NULL
      return(model)
    } else {
      return(0)
    }
  }

  # === K=1 DIRECT PATH: Woodbury identity computation ===
  # For single-matrix prior (K=1) with no null weight, use multivariate_lbf
  # directly instead of mashr's mixture C++ path. This avoids numerical
  # divergence by keeping the optimizer and likelihood on the same code path.
  if (length(model$V_structure) == 1 && model$null_weight == 0 &&
      is.null(model$eigen_cache)) {
    U <- V * model$V_structure[[1]]
    svs <- if (!is.null(model$svs)) model$svs else data$svs
    lbf <- multivariate_lbf(ser_stats$betahat, svs, U)
    lbf[is.na(lbf)] <- 0

    softmax_res <- compute_softmax(lbf, model$pi)
    if (!is.null(l)) {
      model$alpha[l, ]        <- softmax_res$weights
      model$lbf[l]            <- softmax_res$log_sum
      model$lbf_variable[l, ] <- lbf
      # K=1, no null: pi_V_posterior is J x 2 (null=0, non-null=1)
      model$pi_V_posterior[[l]] <- cbind(rep(0, J), rep(1, J))
      # No llik_cache: K=1 EM uses direct formula (fallback path)
      model$llik_cache <- NULL
      return(model)
    } else {
      return(softmax_res$log_sum)
    }
  }

  # === FAST PATH: eigendecomposition precomputation ===
  if (!is.null(model$eigen_cache)) {
    betahat <- ser_stats$betahat                     # J x R
    llik <- loglik_precomputed(betahat, V, model$eigen_cache,
                                ser_stats$BQ_cache)

    # Build pi_full for mixture weights (same as slow path)
    pi_full <- c(model$null_weight, model$pi_V * (1 - model$null_weight))

    # From here: identical to slow path (lbf, alpha, pi_V_posterior)
    lfactors <- apply(llik, 1, max)
    llik_obj <- list(loglik_matrix = llik - lfactors, lfactors = lfactors)
    s_alpha <- matrix(1, 1, 1)
    loglik_null <- compute_null_loglik_from_matrix(llik_obj, s_alpha)
    loglik_alt <- compute_alt_loglik_from_matrix_and_pi(pi_full, llik_obj, s_alpha)
    lbf <- loglik_alt - loglik_null
    if (!is.null(ncol(lbf)) && ncol(lbf) == 1) lbf <- as.vector(lbf)
    lbf[is.na(lbf)] <- 0

    softmax_res <- compute_softmax(lbf, model$pi)
    if (!is.null(l)) {
      model$alpha[l, ]        <- softmax_res$weights
      model$lbf[l]            <- softmax_res$log_sum
      model$lbf_variable[l, ] <- lbf
      d_mat <- t(pi_full * t(exp(llik - lfactors)))
      model$pi_V_posterior[[l]] <- d_mat / rowSums(d_mat)
      model$llik_cache <- llik
      model$per_effect_llik[[l]] <- llik
      # Pass BQ_cache to posterior (avoids recomputing betahat %*% Q_k)
      model$BQ_cache <- ser_stats$BQ_cache
      return(model)
    } else {
      return(softmax_res$log_sum)
    }
  }

  # === SLOW PATH: mashr C++ ===

  # Build full prior for mashr C++ (prepend null component)
  mashr_prior <- build_mashr_prior(model$V_structure, model$pi_V,
                                    model$null_weight, V, R)

  # Prepare mashr inputs
  bhat_t <- t(ser_stats$betahat)                    # R x J
  svs_3d <- matlist2array(ser_stats$shat2)           # R x R x J

  # Get residual correlation from model (updated) or data
  res_corr <- if (!is.null(model$residual_correlation))
    model$residual_correlation else data$residual_correlation
  is_common_cov <- data$is_common_cov

  # Compute J x (K+1) log-likelihood matrix via mashr C++
  llik <- calc_lik_rcpp(
    bhat_t,
    matrix(0, 0, 0),   # sbhat not needed when svs is provided
    res_corr,
    matrix(0, 0, 0),   # no precomputed sigma_rooti
    mashr_prior$xUlist_full_3d,
    svs_3d,
    TRUE,               # log scale
    is_common_cov,
    n_thread
  )$data

  # Warn about non-finite likelihoods
  bad_cols <- which(apply(llik, 2, function(x) any(is.infinite(x))))
  if (length(bad_cols) > 0) {
    warning(paste(
      "Some mixture components result in non-finite likelihoods,",
      "either\n due to numerical underflow/overflow,",
      "or due to invalid covariance matrices",
      paste(bad_cols, collapse = ", ")
    ))
  }

  # Compute lbf (J-vector) using mashr functions
  lfactors <- apply(llik, 1, max)
  llik_obj <- list(loglik_matrix = llik - lfactors, lfactors = lfactors)
  s_alpha <- matrix(1, 1, 1)
  loglik_null <- compute_null_loglik_from_matrix(llik_obj, s_alpha)
  loglik_alt <- compute_alt_loglik_from_matrix_and_pi(
    mashr_prior$pi_full, llik_obj, s_alpha)
  lbf <- loglik_alt - loglik_null
  if (!is.null(ncol(lbf)) && ncol(lbf) == 1) lbf <- as.vector(lbf)
  lbf[is.na(lbf)] <- 0

  # Compute alpha (variable posterior inclusion) from marginalized lbf
  softmax_res <- compute_softmax(lbf, model$pi)

  if (!is.null(l)) {
    model$alpha[l, ]        <- softmax_res$weights
    model$lbf[l]            <- softmax_res$log_sum
    model$lbf_variable[l, ] <- lbf

    # Compute mixture posterior weights: J x (K+1)
    d_mat <- t(mashr_prior$pi_full * t(exp(llik - lfactors)))
    model$pi_V_posterior[[l]] <- d_mat / rowSums(d_mat)

    # Cache llik for EM update and posterior computation
    model$llik_cache <- llik
    model$per_effect_llik[[l]] <- llik

    return(model)
  } else {
    return(softmax_res$log_sum)
  }
}

#' @keywords internal
neg_loglik.mv_individual <- function(data, params, model, V_param, ser_stats, ...) {
  V <- exp(V_param)
  -loglik.mv_individual(data, params, model, V, ser_stats)
}

# =============================================================================
# POSTERIOR MOMENTS
# =============================================================================

#' @keywords internal
calculate_posterior_moments.mv_individual <- function(data, params, model,
                                                       V, l, ...) {
  J <- data$p
  R <- data$R

  # When V ~= 0 or NA, posterior is null: mu = 0, cache = zeros
  if (is.na(V) || V < params$prior_tol) {
    model$mu[l, , ] <- 0
    model$mu2_cache[[l]] <- list(
      bxxb = matrix(0, R, R), vbxxb = 0,
      alpha_mu2_sum = matrix(0, R, R), mu2_diag = matrix(0, J, R)
    )
    model$conditional_lfsr[l] <- list(NULL)
    model$em_cache <- NULL
    return(model)
  }

  # Build reduce_params for inline bxxb/vbxxb accumulation
  alpha_l <- model$alpha[l, ]
  svs_inv <- if (!is.null(model$svs_inv)) model$svs_inv else data$svs_inv
  v_inv   <- get_v_inv(data, model)
  reduce_params <- list(alpha = alpha_l, d = data$d, svs_inv = svs_inv,
                        v_inv = v_inv)

  # === FAST PATH: eigendecomposition precomputation ===
  if (!is.null(model$eigen_cache)) {
    betahat <- compute_betahat(data, model$residuals)  # J x R

    # Compute P(j|k) weights for EM (same formula as slow path).
    if (!is.null(params$estimate_prior_method) && params$estimate_prior_method == "EM"
        && !is.null(model$llik_cache)) {
      em_var_wt <- compute_variable_posterior_weights(model$pi, model$llik_cache)
    } else {
      em_var_wt <- NULL
    }

    post <- posterior_precomputed(betahat, V, model$eigen_cache,
                                  model$pi_V_posterior[[l]],
                                  em_var_wt = em_var_wt,
                                  reduce_params = reduce_params,
                                  BQ_cache = model$BQ_cache)
    model$BQ_cache <- NULL  # free memory after use

    model$mu[l, , ] <- post$post_mean
    # Store reduced statistics (no J x R x R array)
    model$mu2_cache[[l]] <- list(
      bxxb = post$bxxb, vbxxb = post$vbxxb,
      alpha_mu2_sum = post$alpha_mu2_sum, mu2_diag = post$mu2_diag
    )
    model$conditional_lfsr[[l]] <- compute_lfsr(post$post_neg, post$post_zero)

    if (!is.null(params$estimate_prior_method) && params$estimate_prior_method == "EM") {
      model$em_cache <- list(
        prior_scale_em_update = post$prior_scale_em_update,
        alpha_mu2_sum = post$alpha_mu2_sum
      )
    }
    return(model)
  }

  # === SLOW PATH: mashr C++ ===

  betahat <- compute_betahat(data, model$residuals)  # J x R
  n_thread <- if (!is.null(params$n_thread)) params$n_thread else 1L

  # Build full prior arrays (prepend null)
  mashr_prior <- build_mashr_prior(model$V_structure, model$pi_V,
                                    model$null_weight, V, R)

  # Build inverse array (prepend null inverse = zero).
  # V_structure_inv is computed lazily (only for slow path).
  if (is.null(model$V_structure_inv)) {
    inv_list <- lapply(model$V_structure, pseudo_inverse)
    model$V_structure_inv <- matlist2array(lapply(inv_list, `[[`, "inv"))
  }
  null_inv <- array(0, c(R, R, 1))
  V_inv_full <- abind::abind(null_inv, model$V_structure_inv, along = 3)

  # Prepare inputs
  bhat_t <- t(betahat)
  svs_inv_3d <- matlist2array(svs_inv)
  res_corr <- if (!is.null(model$residual_correlation))
    model$residual_correlation else data$residual_correlation
  is_common_cov <- data$is_common_cov

  # Compute P(j|k) posterior variable weights for EM (if needed).
  if (!is.null(params$estimate_prior_method) && params$estimate_prior_method == "EM"
      && !is.null(model$llik_cache)) {
    var_post_wt <- compute_variable_posterior_weights(model$pi, model$llik_cache)
  } else {
    var_post_wt <- matrix(0, 0, 0)
  }

  # Call mashr C++ for mixture-weighted posterior
  post <- calc_sermix_rcpp(
    bhat_t,
    matrix(0, 0, 0),          # sbhat not needed when svs_inv provided
    res_corr,
    svs_inv_3d,
    mashr_prior$xUlist_full_3d,
    V_inv_full,                # original (unscaled) inverses
    0,                         # no precomputed U0
    t(model$pi_V_posterior[[l]]),
    var_post_wt,
    is_common_cov,
    n_thread
  )

  # Store posterior means
  model$mu[l, , ] <- post$post_mean  # J x R

  # Compute reduced statistics from post_cov (R x R x J) + post_mean
  # instead of storing full J x R x R mu2 array
  bxxb          <- matrix(0, R, R)
  alpha_mu2_sum <- matrix(0, R, R)
  vbxxb         <- 0
  mu2_diag      <- matrix(0, J, R)
  for (j in seq_len(J)) {
    mu2_j <- post$post_cov[, , j] + tcrossprod(post$post_mean[j, ])
    a_j <- alpha_l[j]
    bxxb          <- bxxb + data$d[j] * a_j * mu2_j
    alpha_mu2_sum <- alpha_mu2_sum + a_j * mu2_j
    vbxxb         <- vbxxb + a_j * sum(svs_inv[[min(j, length(svs_inv))]] * mu2_j)
    mu2_diag[j, ] <- diag(mu2_j)
  }
  model$mu2_cache[[l]] <- list(
    bxxb = bxxb, vbxxb = vbxxb,
    alpha_mu2_sum = alpha_mu2_sum, mu2_diag = mu2_diag
  )

  # LFSR from posterior sign probabilities
  model$conditional_lfsr[[l]] <- compute_lfsr(post$post_neg, post$post_zero)

  # Cache EM statistics (no recomputation needed)
  if (!is.null(params$estimate_prior_method) && params$estimate_prior_method == "EM") {
    model$em_cache <- list(
      prior_scale_em_update = post$prior_scale_em_update,
      alpha_mu2_sum = alpha_mu2_sum
    )
  }

  return(model)
}

# =============================================================================
# KL DIVERGENCE
# =============================================================================

#' @keywords internal
compute_kl.mv_individual <- function(data, params, model, l) {
  # KL = -lbf + SER_posterior_e_loglik (relative to null).
  # Also stores bxxb and vbxxb as precomputed quantities
  result <- SER_posterior_e_loglik.mv_individual(data, params, model, l)
  model$KL[l] <- -model$lbf[l] + result$eloglik
  model$bxxb[[l]] <- result$bxxb
  model$vbxxb[l] <- result$vbxxb
  return(model)
}

#' @keywords internal
SER_posterior_e_loglik.mv_individual <- function(data, params, model, l) {
  # Returns a list with:
  #   eloglik: expected log-likelihood for effect l relative to null
  #   bxxb:    R x R matrix, sum_j d_j alpha_j mu2_j (used by estimate_residual_variance)
  #   vbxxb:   scalar, sum_j d_j alpha_j tr(v_inv * mu2_j) (used by ELBO)
  v_inv <- get_v_inv(data, model)
  alpha_l <- model$alpha[l, ]
  mu_l    <- model$mu[l, , , drop = TRUE]  # J x R
  R <- data$R

  # Weighted contribution: X %*% (alpha * mu)
  Xb_l <- compute_Xb(data, alpha_l * mu_l)  # N x R

  # E1 = tr(v_inv * (B1'XtR + XtR'B1)) = 2 * tr(v_inv * B1'XtR)
  # (factor of 2 from v_inv symmetry)
  # For per-outcome (3d) methods, use V_i^{-1}-weighted inner product.
  if (!is.null(data$miss3d)) {
    VinvXb_l <- compute_VinvR_3d(data, Xb_l)
    term1 <- sum(model$raw_residuals * VinvXb_l)
  } else if (is.matrix(v_inv)) {
    term1 <- sum(model$raw_residuals * (Xb_l %*% v_inv))
  } else {
    term1 <- v_inv * sum(model$raw_residuals * Xb_l)
  }

  # Use precomputed bxxb/vbxxb from calculate_posterior_moments
  cache <- model$mu2_cache[[l]]
  bxxb_l  <- cache$bxxb
  vbxxb_l <- cache$vbxxb

  eloglik <- 0.5 * (2 * term1 - vbxxb_l)
  return(list(eloglik = eloglik, bxxb = bxxb_l, vbxxb = vbxxb_l))
}

# =============================================================================
# FITTED VALUES
# =============================================================================

#' @keywords internal
update_fitted_values.mv_individual <- function(data, params, model, l, ...) {
  # Efficient update: Xr = fitted_without_l + X %*% (alpha_l * mu_l)
  b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]
  model$Xr <- model$fitted_without_l + compute_Xb(data, b_l)
  return(model)
}

# =============================================================================
# VARIANCE COMPONENT UPDATES
# =============================================================================

#' @keywords internal
update_variance_components.mv_individual <- function(data, params, model, ...) {
  list(sigma2 = estimate_residual_variance_mv(data, model))
}

#' @keywords internal
update_model_variance.mv_individual <- function(data, params, model) {
  # Update residual variance if requested.
  #
  # NOTE: For the 3D missing-data path (data$miss3d != NULL), sigma2 is
  # estimated via an OUTER loop in mvsusie_core, not here.  The outer loop
  # passes estimate_residual_variance=FALSE to the inner IBSS, so this
  # branch is only reached for the non-missing and imputation paths.
  if (isTRUE(params$estimate_residual_variance)) {
    model$sigma2 <- estimate_residual_variance_mv(data, model)

    # Recompute svs/svs_inv and store in MODEL (not data, since data is immutable).
    # These override the data's precomputed values during SER computation.
    model$residual_variance_inv <- invert_via_chol(model$sigma2)$inv
    model$residual_correlation <- cov2cor(model$sigma2)

    if (!is.null(data$miss3d)) {
      # 3D missing-data path: recompute ALL V-dependent quantities via
      # set_residual_variance_3d (Vinv, elbo_const, svs, svs_inv).
      # Store results in model; inject into data at each call site.
      # (Retained for future use; currently the outer loop in mvsusie_core
      # handles sigma2 updates for 3D data instead of this code path.)
      updated <- set_residual_variance_3d(data, model$sigma2,
                                          center = FALSE, scale = FALSE)
      model$miss3d_Vinv       <- updated$miss3d$Vinv
      model$miss3d_elbo_const <- updated$miss3d$elbo_const
      model$svs               <- updated$svs
      model$svs_inv           <- updated$svs_inv
      model$is_common_cov     <- updated$is_common_cov
    } else if (data$is_common_cov) {
      # For common-cov, store single copy to save O(J*R^2) memory
      svs_one <- model$sigma2 / data$d[1]
      svs_one[is.nan(svs_one) | is.infinite(svs_one)] <- 1e6
      model$svs <- list(svs_one)
      model$svs_inv <- list(model$residual_variance_inv * data$d[1])
    } else {
      model$svs <- lapply(seq_len(data$p), function(j) {
        res <- model$sigma2 / data$d[j]
        res[is.nan(res) | is.infinite(res)] <- 1e6
        res
      })
      model$svs_inv <- lapply(seq_len(data$p), function(j) {
        model$residual_variance_inv * data$d[j]
      })
    }

    # Recompute eigendecomposition cache if precomputation is active
    if (!is.null(model$eigen_cache))
      model$eigen_cache <- precompute_eigen_cache(
        model$svs, model$V_structure,
        if (!is.null(model$is_common_cov)) model$is_common_cov
        else data$is_common_cov)
  }

  # Update mixture prior weights and prune near-zero components
  model <- update_mixture_weights_and_prune(model, params)

  return(model)
}

#' @keywords internal
update_derived_quantities.mv_individual <- function(data, params, model) {
  return(model)
}

# =============================================================================
# OBJECTIVE (ELBO)
# =============================================================================

#' @keywords internal
get_objective.mv_individual <- function(data, params, model, ...) {
  compute_multivariate_elbo(data, model) - sum(model$KL, na.rm = TRUE)
}

# =============================================================================
# CONVERGENCE
# =============================================================================

# =============================================================================
# EM PRIOR VARIANCE UPDATE
# =============================================================================

#' @keywords internal
em_update_prior_variance.mv_individual <- function(data, params, model,
                                                     alpha, moments, V_init) {
  # When V was ~0 last iteration, caches are NULL (loglik + posterior both
  # returned early). The posterior is all zeros, so the EM M-step gives 0.
  if (is.null(model$em_cache)) return(0)

  # Unified EM update using cached stats from mashr calc_sermix_rcpp.
  # prior_scale_em_update is a (K+1)-vector (trace-based scalar per component),
  # V_structure_rank is a K-vector (non-null components only).
  # SER posterior mixture weights determine the contribution of each component.

  if (!is.null(model$llik_cache)) {
    # Compute SER posterior mixture weights: (K+1)-vector
    pi_full <- c(model$null_weight, model$pi_V * (1 - model$null_weight))
    lbf_mat <- model$llik_cache - model$llik_cache[, 1]  # J x (K+1)
    ser_lbf <- apply(lbf_mat, 2, function(x) compute_softmax(x, model$pi)$log_sum)
    SER_post_mix_wt <- compute_softmax(ser_lbf, pi_full)$weights

    # EM update: skip null component (index 1)
    em_update <- model$em_cache$prior_scale_em_update
    V_new <- sum(SER_post_mix_wt[-1] * em_update[-1]) /
             sum(SER_post_mix_wt[-1] * model$V_structure_rank)
    return(max(0, V_new))
  }

  # Fallback: simple EM using marginalized posteriors (for single-matrix K=1).
  # Compute V_structure_inv on demand (not always available on eigen path).
  if (is.null(model$V_structure_inv)) {
    inv_list <- lapply(model$V_structure, pseudo_inverse)
    model$V_structure_inv <- matlist2array(lapply(inv_list, `[[`, "inv"))
  }
  mu2 <- model$em_cache$alpha_mu2_sum
  scalar <- sum(diag(model$V_structure_inv[, , 1] %*% mu2)) /
            model$V_structure_rank[1]
  return(max(0, scalar))
}

# =============================================================================
# MIXTURE WEIGHT UPDATE
# =============================================================================

#' Inner EM update for mixture prior weights pi_V.
#'
#' Iteratively maximizes the collapsed objective
#'   F(pi) = sum_{l,j} alpha[l,j] * log(sum_k pi[k] * exp(llik[l,j,k]))
#' via EM on the stacked (L*J) x (K+1) log-likelihood matrix, weighted
#' by alpha[l,j]. Core loop implemented in C++ (inner_em_rcpp).
#'
#' @param model Model list with per_effect_llik, alpha, pi_V, null_weight.
#' @param update_null Logical; if TRUE, also update null_weight.
#' @param max_inner_iter Maximum inner EM iterations.
#' @param inner_tol Convergence tolerance on max absolute change in pi.
#' @return Updated model with new pi_V (and optionally null_weight).
#' @keywords internal
update_mixture_weights_em <- function(model, update_null = FALSE,
                                       max_inner_iter = 100,
                                       inner_tol = 1e-10) {
  L <- nrow(model$alpha)

  # Stack per-effect llik matrices and alpha weights
  llik_list <- list()
  alpha_list <- list()
  for (l in seq_len(L)) {
    if (!is.null(model$per_effect_llik[[l]])) {
      llik_list[[length(llik_list) + 1]] <- model$per_effect_llik[[l]]
      alpha_list[[length(alpha_list) + 1]] <- model$alpha[l, ]
    }
  }
  if (length(llik_list) == 0) return(model)

  llik_combined <- do.call(rbind, llik_list)
  alpha_weights <- unlist(alpha_list)
  pi_full <- c(model$null_weight, model$pi_V * (1 - model$null_weight))

  # Floor zero weights so EM can discover all components
  K_plus_1 <- length(pi_full)
  pi_full <- pmax(pi_full, 1e-10 / K_plus_1)
  pi_full <- pi_full / sum(pi_full)

  result <- inner_em_rcpp(llik_combined, alpha_weights, pi_full,
                         max_inner_iter, inner_tol)
  pi_full <- result$pi

  if (update_null) model$null_weight <- pi_full[1]
  non_null <- pi_full[-1]
  s <- sum(non_null)
  if (s > 0) model$pi_V <- non_null / s
  return(model)
}

#' Mixture weight update via mixsqp (sequential quadratic programming).
#'
#' Solves the alpha-weighted ML problem on the stacked (L*J) x (K+1)
#' log-likelihood matrix via mixsqp (Kim, Carbonetto, Stephens 2020).
#'
#' @param model Model list with per_effect_llik, alpha, pi_V, null_weight.
#' @param update_null Logical; if TRUE, also update null_weight.
#' @return Updated model with new pi_V (and optionally null_weight).
#' @importFrom mixsqp mixsqp
#' @keywords internal
update_mixture_weights_mixsqp <- function(model, update_null = FALSE) {
  L <- nrow(model$alpha)

  # Stack per-effect llik matrices and alpha weights
  llik_list <- list()
  alpha_list <- list()
  for (l in seq_len(L)) {
    if (!is.null(model$per_effect_llik[[l]])) {
      llik_list[[length(llik_list) + 1]] <- model$per_effect_llik[[l]]
      alpha_list[[length(alpha_list) + 1]] <- model$alpha[l, ]
    }
  }
  if (length(llik_list) == 0) return(model)

  llik_combined <- do.call(rbind, llik_list)
  alpha_weights <- unlist(alpha_list)

  out <- withCallingHandlers(
    mixsqp(llik_combined, w = alpha_weights, log = TRUE,
            control = list(verbose = FALSE)),
    warning = function(w) {
      warn_once("mixsqp_zero_cols", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  w_new <- out$x

  if (update_null) model$null_weight <- w_new[1]
  non_null <- w_new[-1]
  s <- sum(non_null)
  if (s > 0) model$pi_V <- non_null / s
  return(model)
}

#' Update mixture prior weights using either mixsqp (default) or EM.
#'
#' @param model Model list.
#' @param method Character; "mixsqp" (default) or "EM".
#' @param update_null Logical; if TRUE, also update null_weight.
#' @return Updated model.
#' @keywords internal
update_mixture_weights <- function(model, method = "mixsqp",
                                    update_null = FALSE) {
  if (method == "mixsqp") {
    update_mixture_weights_mixsqp(model, update_null)
  } else {
    update_mixture_weights_em(model, update_null)
  }
}

#' Drop mixture components with weight below threshold.
#'
#' Removes corresponding entries from V_structure, V_structure_3d,
#' V_structure_inv, V_structure_rank, pi_V, pi_V_posterior columns,
#' per_effect_llik columns, and eigen_cache entries. We check
#' every 10 iterations with threshold 1e-8.
#'
#' @param model Model list containing mixture prior structure.
#' @param threshold Minimum weight to keep a component.
#' @return Updated model with pruned components.
#' @keywords internal
prune_mixture_components <- function(model, threshold = 1e-8) {
  keep <- which(model$pi_V >= threshold)
  if (length(keep) == length(model$pi_V)) return(model)

  model$pi_V <- model$pi_V[keep]
  model$pi_V <- model$pi_V / sum(model$pi_V)
  model$V_structure <- model$V_structure[keep]
  if (!is.null(model$V_structure_3d))
    model$V_structure_3d <- model$V_structure_3d[, , keep, drop = FALSE]
  if (!is.null(model$V_structure_inv))
    model$V_structure_inv <- model$V_structure_inv[, , keep, drop = FALSE]
  model$V_structure_rank <- model$V_structure_rank[keep]

  # Trim pi_V_posterior and per_effect_llik columns (offset by 1 for null)
  keep_full <- c(1, keep + 1)  # null column + kept non-null columns
  for (l in seq_along(model$pi_V_posterior)) {
    if (!is.null(model$pi_V_posterior[[l]])) {
      model$pi_V_posterior[[l]] <- model$pi_V_posterior[[l]][, keep_full, drop = FALSE]
      # Renormalize rows
      rs <- rowSums(model$pi_V_posterior[[l]])
      model$pi_V_posterior[[l]] <- model$pi_V_posterior[[l]] / rs
    }
    if (!is.null(model$per_effect_llik[[l]])) {
      model$per_effect_llik[[l]] <- model$per_effect_llik[[l]][, keep_full, drop = FALSE]
    }
  }

  # Trim eigen_cache if active: subset the components list (one per non-null
  # component) and rebuild the cache structure.
  if (!is.null(model$eigen_cache)) {
    model$eigen_cache$components <- model$eigen_cache$components[keep]
  }

  return(model)
}

# =============================================================================
# PRIOR VALIDATION
# =============================================================================

#' @keywords internal
validate_prior.mv_individual <- function(data, params, model, ...) {
  invisible(TRUE)
}

# =============================================================================
# EFFECT TRIMMING
# =============================================================================

#' @keywords internal
trim_null_effects.mv_individual <- function(data, params, model) {
  K_plus_1 <- length(model$pi_V) + 1  # K non-null + 1 null
  J <- ncol(model$alpha)
  R <- dim(model$mu)[3]
  for (l in seq_len(nrow(model$alpha))) {
    # V is now a numeric scalar per effect
    if (is.na(model$V[l]) || abs(model$V[l]) < params$prior_tol) {
      model$V[l]               <- 0
      model$alpha[l, ]         <- 1 / J
      model$mu[l, , ]          <- 0
      model$mu2_cache[[l]]     <- list(
        bxxb = matrix(0, R, R), vbxxb = 0,
        alpha_mu2_sum = matrix(0, R, R), mu2_diag = matrix(0, J, R)
      )
      model$lbf_variable[l, ]  <- 0
      model$lbf[l]             <- 0
      model$KL[l]              <- 0
      # Reset mixture-specific fields (use list(NULL) to avoid removing element)
      model$pi_V_posterior[l]    <- list(NULL)
      model$conditional_lfsr[l] <- list(NULL)
    }
  }
  return(model)
}

# =============================================================================
# OUTPUT / POST-PROCESSING
# =============================================================================

#' @keywords internal
get_scale_factors.mv_individual <- function(data, params, ...) {
  data$csd
}

#' @keywords internal
get_intercept.mv_individual <- function(data, params, model, ...) {
  # Per-outcome (3d) missing data methods use their own intercept recovery
  if (!is.null(data$miss3d)) return(get_intercept_3d(data, model))

  b_sum <- compute_posterior_mean_sum_from_model(model)
  # Rescale: coefficients on original scale
  coefs_original <- b_sum / data$csd
  # Intercept = Y_mean - X_mean' %*% coefs
  intercept <- data$Y_mean - as.vector(t(data$cm) %*% coefs_original)
  # Add accumulated centering adjustment from imputation
  if (!is.null(model$Y_mean_adjustment)) {
    intercept <- intercept + model$Y_mean_adjustment
  }
  return(intercept)
}

#' Get fitted values on the original Y scale.
#'
#' Returns \code{X_std \%*\% b_sum + Y_mean}, i.e. fitted values on the
#' original (uncentered) Y scale.
#' This matches \code{susieR::get_fitted.individual()} which returns
#' \code{Xr + mean_y} when \code{intercept = TRUE}.
#'
#' To obtain fitted on the \strong{standardized} (centered) scale instead,
#' subtract the column means of Y:
#' \code{fitted_std <- fitted - matrix(Y_mean, n, R, byrow = TRUE)}
#' or equivalently use \code{model$Xr} directly before this function is called.
#'
#' @keywords internal
get_fitted.mv_individual <- function(data, params, model, ...) {
  # Xr = X_centered_scaled %*% b_sum (standardized scale, no intercept)
  # Adding Y_mean gives fitted on the original Y scale (matches susieR)
  model$Xr + matrix(data$Y_mean, data$n, data$R, byrow = TRUE)
}

#' @keywords internal
get_cs.mv_individual <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr))
    return(NULL)
  susie_get_cs(model,
    X            = data$X,
    coverage     = params$coverage,
    min_abs_corr = params$min_abs_corr
  )
}

#' @keywords internal
get_variable_names.mv_individual <- function(data, model, ...) {
  vnames <- colnames(data$X)
  if (!is.null(vnames)) {
    colnames(model$alpha)        <- vnames
    colnames(model$lbf_variable) <- vnames
  }
  return(model)
}

#' @keywords internal
get_zscore.mv_individual <- function(data, params, model, ...) {
  NULL
}

#' @keywords internal
cleanup_model.mv_individual <- function(data, params, model, ...) {
  # Final pruning of near-zero mixture components at convergence
  if (length(model$pi_V) > 1) {
    model <- prune_mixture_components(model, threshold = 1e-8)
  }
  model$residuals        <- NULL
  model$fitted_without_l <- NULL
  model$raw_residuals    <- NULL
  model$llik_cache       <- NULL
  model$per_effect_llik   <- NULL
  model$em_cache         <- NULL
  model$eigen_cache      <- NULL
  model$ibss_iter        <- NULL
  # Clean up imputation fields (large N x R matrices)
  model$Y_imputed          <- NULL
  model$Y_cov              <- NULL
  model$sum_neg_ent_Y_miss <- NULL
  model$pattern_cache      <- NULL
  model$Y_mean_adjustment  <- NULL
  return(model)
}

# =============================================================================
# SHARED HELPERS (primarily used by mv_individual methods)
# =============================================================================

# Helper: get residual variance inverse (model overrides data when sigma2 is estimated)
get_v_inv <- function(data, model) {
  if (!is.null(model$residual_variance_inv)) model$residual_variance_inv
  else data$residual_variance_inv
}

# Helper: compute sum of posterior means from model
compute_posterior_mean_sum_from_model <- function(model) {
  L <- nrow(model$alpha)
  J <- ncol(model$alpha)
  R <- dim(model$mu)[3]
  b <- matrix(0, J, R)
  for (l in seq_len(L)) {
    b <- b + drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]
  }
  return(b)
}

# =============================================================================
# INDIVIDUAL-SPECIFIC HELPERS
# =============================================================================

# Estimate residual variance for dense data
#
# Computes V = E_q[(Y - XB)'(Y - XB)] / n  where B = sum_l b_l.
#
# Under the variational posterior q (with independent single effects):
#   E_q[B' X'X B] = sum_l E_q[b_l' X'X b_l] + sum_{l!=l'} E[b_l]' X'X E[b_{l'}]
#
# For single effect l, only one variable is active, so cross-variable terms
# vanish: E_q[b_l' X'X b_l] = sum_j d_j * alpha_lj * mu2_lj  (="bxxb_l").
#
# Writing R_hat = Y - X * B_sum (residual at posterior mean), we get:
#   E_q[R'R] = R_hat'R_hat + sum_l (bxxb_l - B1_l' X'X B1_l)
#
# where B1_l = E_q[b_l] (J x R, rows = alpha_lj * mu_lj).
# NOTE: B1_l' X'X B1_l uses the FULL X'X, not just the diagonal.
estimate_residual_variance_mv <- function(data, model) {
  # 3D missing-data path: use per-outcome X_3d and conditional imputation
  if (!is.null(data$miss3d))
    return(estimate_residual_variance_3d(data, model))

  b_sum <- compute_posterior_mean_sum_from_model(model)
  Xb <- data$X %*% b_sum

  # Use imputed Y when available (R>1 missing data)
  Y <- if (!is.null(model$Y_imputed)) model$Y_imputed else data$Y
  R_mat <- Y - Xb

  # Use full sample size when imputation is active;
  # use effective sample size for R=1 missing data (complete-case).
  N <- if (data$any_missing && is.null(model$Y_imputed)) data$n_obs else data$n
  R <- data$R

  # Zero out missing entries only for R=1 (complete-case)
  if (data$any_missing && is.null(model$Y_imputed)) {
    R_mat[data$Y_na] <- 0
    Xb[data$Y_na] <- 0
  }

  # R_hat' R_hat
  E_RtR <- crossprod(R_mat)
  bxxb <- Reduce("+", model$bxxb)

  # B1_l' X'X B1_l  (uses full X'X via X %*% B1_l)
  b1_XtX_b1 <- matrix(0, R, R)
  for (l in seq_len(nrow(model$alpha))) {
    B1_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
    if (!is.matrix(B1_l)) B1_l <- matrix(B1_l, ncol = R)
    XB1_l <- data$X %*% B1_l  # n x R
    if (data$any_missing && is.null(model$Y_imputed)) XB1_l[data$Y_na] <- 0
    b1_XtX_b1 <- b1_XtX_b1 + crossprod(XB1_l)  # B1_l' X'X B1_l
  }

  V_est <- (E_RtR + bxxb - b1_XtX_b1) / N

  # Add Y_cov for VEM based imputation uncertainty
  # V = (R'R + var_part_ERSS + Y_cov) / n
  if (!is.null(model$Y_cov)) {
    V_est <- V_est + model$Y_cov / N
  }

  V_est <- (V_est + t(V_est)) / 2

  # Enforce positive-definiteness when imputation is active
  if (!is.null(model$Y_cov)) {
    V_est <- makePD(V_est, 1e-10)
  }

  return(V_est)
}

# Multivariate ELBO expected log-likelihood (dense)
compute_multivariate_elbo <- function(data, model) {
  # Per-outcome (3d) missing data methods have their own ELBO
  if (!is.null(data$miss3d)) return(compute_elbo_3d(data, model))

  # Use full sample size when imputation is active (all obs contribute);
  # use effective sample size for R=1 missing data (complete-case).
  N <- if (data$any_missing && is.null(model$Y_imputed)) data$n_obs else data$n
  R <- data$R
  v_inv <- get_v_inv(data, model)

  # Use imputed Y when available (R>1 missing data)
  Y <- if (!is.null(model$Y_imputed)) model$Y_imputed else data$Y

  # Constant: -N*R/2*log(2*pi) - N/2*log|sigma2|
  loglik <- -N * R / 2 * log(2 * pi)
  if (is.matrix(model$sigma2)) {
    loglik <- loglik - N / 2 * log_det_sym(model$sigma2)
  } else {
    loglik <- loglik - N * R / 2 * log(model$sigma2)
  }

  # ESSR = tr(v_inv * (Y-Xr)'(Y-Xr)) - sum_l tr(v_inv * b_l' X'X b_l)
  #        + sum_l sum_j d_j alpha_lj tr(v_inv * mu2_lj)

  # 1. Raw residual RSS: tr(v_inv * (Y-Xr)'(Y-Xr))
  b_sum <- compute_posterior_mean_sum_from_model(model)
  Xb <- data$X %*% b_sum
  R_mat <- Y - Xb
  # Zero out missing entries only for R=1 (complete-case)
  if (data$any_missing && is.null(model$Y_imputed)) R_mat[data$Y_na] <- 0
  if (is.matrix(v_inv)) {
    essr <- sum(R_mat * (R_mat %*% v_inv))
  } else {
    essr <- v_inv * sum(R_mat^2)
  }

  # 2. Subtract first moment^2 per effect: -sum_l tr(v_inv * b_l' X'X b_l)
  for (l in seq_len(nrow(model$alpha))) {
    b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
    Xb_l <- data$X %*% b_l  # N x R
    # Zero out missing rows only for R=1 (complete-case)
    if (data$any_missing && is.null(model$Y_imputed)) Xb_l[data$Y_na] <- 0
    if (is.matrix(v_inv)) {
      essr <- essr - sum(Xb_l * (Xb_l %*% v_inv))
    } else {
      essr <- essr - v_inv * sum(Xb_l^2)
    }
  }

  # 3. Add vbxxb: use precomputed per-effect values from compute_kl
  #    (d already accounts for missing data in observed-row version)
  essr <- essr + sum(model$vbxxb)

  result <- loglik - 0.5 * essr

  # ELBO correction for VEM based imputation uncertainty:
  # ELBO_miss = ELBO_complete - 0.5 * tr(V^-1 Y_cov) - sum_neg_ent
  if (!is.null(model$Y_cov)) {
    result <- result - 0.5 * sum(v_inv * model$Y_cov)
    result <- result - model$sum_neg_ent_Y_miss
  }

  return(result)
}
