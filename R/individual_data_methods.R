# S3 method implementations for multivariate individual-level data.
#
# Methods for the mv_individual class, providing implementations of
# susieR's internal generics for dense (X, Y) data.

#' @importFrom susieR susie_get_cs
#' @importFrom ashr compute_lfsr
#' @importFrom utils modifyList
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
#' @importFrom utils modifyList
ibss_initialize.mv_individual <- function(data, params) {
  # get_var_y is an internal susieR generic; our methods are registered
  # in .onLoad so S3 dispatch works, but we need to find the generic.
  get_var_y <- get("get_var_y", envir = asNamespace("susieR"))
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

    # Create fresh model with correct dimensions
    # initialize_susie_model is an internal susieR generic; use S3 dispatch
    init_fn <- get("initialize_susie_model", envir = asNamespace("susieR"))
    mat_init <- init_fn(data, params, var_y)

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
    # Create fresh model
    init_fn <- get("initialize_susie_model", envir = asNamespace("susieR"))
    mat_init <- init_fn(data, params, var_y)
  }

  # Initialize fitted values and null index
  init_fitted_fn <- get("initialize_fitted", envir = asNamespace("susieR"))
  fitted     <- init_fitted_fn(data, mat_init)
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

  return(model)
}

# =============================================================================
# RESIDUALS
# =============================================================================

#' @keywords internal
compute_residuals.mv_individual <- function(data, params, model, l, ...) {
  # Fitted for effect l: X %*% (alpha_l * mu_l)
  b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
  Xb_l <- data$X %*% b_l  # N x R

  # Residual: Y - (Xr - Xb_l) = Y - Xr + Xb_l
  R_mat <- data$Y - model$Xr + Xb_l  # N x R

  # Zero out missing observations so they don't contribute to X'R
  if (data$Y_has_missing) {
    R_mat[data$Y_missing] <- 0
  }

  # X'R for SER computation
  XtR <- crossprod(data$X, R_mat)  # J x R

  model$residuals       <- XtR
  model$fitted_without_l <- model$Xr - Xb_l
  model$raw_residuals   <- R_mat
  return(model)
}

# =============================================================================
# SER STATISTICS
# =============================================================================

#' @keywords internal
compute_ser_statistics.mv_individual <- function(data, params, model, l, ...) {
  # betahat_j = (X'R)_j / d_j, a J x R matrix
  betahat <- as.matrix(model$residuals / data$d)

  # Use model's svs/svs_inv if available (updated after residual variance change),
  # otherwise fall back to data's precomputed values
  svs     <- if (!is.null(model$svs)) model$svs else data$svs
  svs_inv <- if (!is.null(model$svs_inv)) model$svs_inv else data$svs_inv

  # Moment-based warm start for optim (analogous to susieR's log(max(betahat^2 - shat2, 1))):
  # For each variable j, compute tr(betahat_j betahat_j') - tr(S_j) as a scalar
  # measure of signal strength, then take the max across variables.
  signal <- sapply(seq_len(nrow(betahat)), function(j) {
    sum(betahat[j, ]^2) - sum(diag(svs[[j]]))
  })
  optim_init <- log(max(c(max(signal, na.rm = TRUE), 1)))

  list(
    betahat      = betahat,
    shat2        = svs,            # list of J RxR matrices
    shat2_inv    = svs_inv,        # list of J RxR matrices
    optim_init   = optim_init,
    optim_bounds = c(-30, 15),
    optim_scale  = "log"
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
  if (is.na(V) || V < 1e-15) {
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

  # Build full prior for mashr C++ (prepend null component)
  mashr_prior <- build_mashr_prior(model$V_structure_3d, model$pi_V,
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

  # When V ~= 0 or NA, posterior is null: mu = 0, mu2 = 0
  if (is.na(V) || V < 1e-15) {
    model$mu[l, , ] <- 0
    for (j in seq_len(J)) {
      model$mu2[[l]][j, , ] <- 0
    }
    model$conditional_lfsr[l] <- list(NULL)
    model$em_cache <- NULL
    return(model)
  }

  betahat <- model$residuals / data$d  # J x R
  n_thread <- if (!is.null(params$n_thread)) params$n_thread else 1L

  # Build full prior arrays (prepend null)
  mashr_prior <- build_mashr_prior(model$V_structure_3d, model$pi_V,
                                    model$null_weight, V, R)

  # Build inverse array (prepend null inverse = zero)
  null_inv <- array(0, c(R, R, 1))
  V_inv_full <- abind::abind(null_inv, model$V_structure_inv, along = 3)

  # Prepare inputs
  bhat_t <- t(betahat)
  svs_inv <- if (!is.null(model$svs_inv)) model$svs_inv else data$svs_inv
  svs_inv_3d <- matlist2array(svs_inv)
  res_corr <- if (!is.null(model$residual_correlation))
    model$residual_correlation else data$residual_correlation
  is_common_cov <- data$is_common_cov

  # Compute variable posterior weights for EM (if needed)
  # var_post_wt[k,j] = alpha[l,j] * pi_V_posterior[j,k]
  if (!is.null(params$estimate_prior_method) && params$estimate_prior_method == "EM"
      && !is.null(model$pi_V_posterior[[l]])) {
    var_post_wt <- t(model$alpha[l, ] * model$pi_V_posterior[[l]])  # (K+1) x J
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

  # Store marginalized posteriors
  model$mu[l, , ] <- post$post_mean  # J x R
  for (j in seq_len(J)) {
    model$mu2[[l]][j, , ] <- post$post_cov[, , j] + tcrossprod(post$post_mean[j, ])
  }

  # LFSR from posterior sign probabilities
  model$conditional_lfsr[[l]] <- compute_lfsr(post$post_neg, post$post_zero)

  # Cache EM statistics (no recomputation needed)
  if (!is.null(params$estimate_prior_method) && params$estimate_prior_method == "EM") {
    model$em_cache <- list(
      prior_scale_em_update = post$prior_scale_em_update
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
  Xb_l <- data$X %*% (alpha_l * mu_l)  # N x R

  # E1 = tr(v_inv * (B1'XtR + XtR'B1)) = 2 * tr(v_inv * B1'XtR)
  # (factor of 2 from v_inv symmetry)
  if (is.matrix(v_inv)) {
    term1 <- sum(model$raw_residuals * (Xb_l %*% v_inv))
  } else {
    term1 <- v_inv * sum(model$raw_residuals * Xb_l)
  }

  bxxb_l <- matrix(0, R, R)
  vbxxb_l <- 0
  for (j in seq_len(data$p)) {
    mu2_j <- model$mu2[[l]][j, , , drop = FALSE]
    dim(mu2_j) <- c(R, R)
    pb2_j <- alpha_l[j] * mu2_j   # alpha_j * E[b_j b_j' | j active]
    bxxb_l <- bxxb_l + data$d[j] * pb2_j
    vbxxb_l <- vbxxb_l + data$d[j] * sum(v_inv * pb2_j)
  }

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
  model$Xr <- model$fitted_without_l + data$X %*% b_l
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
  # Call the estimation function directly (the update_variance_components
  # generic lives in susieR's unexported namespace and isn't accessible here).
  model$sigma2 <- estimate_residual_variance_mv(data, model)

  # Recompute svs/svs_inv and store in MODEL (not data, since data is immutable).
  # These override the data's precomputed values during SER computation.
  model$residual_variance_inv <- invert_via_chol(model$sigma2)$inv
  model$residual_correlation <- cov2cor(model$sigma2)
  model$svs <- lapply(seq_len(data$p), function(j) {
    res <- model$sigma2 / data$d[j]
    res[is.nan(res) | is.infinite(res)] <- 1e6
    res
  })
  model$svs_inv <- lapply(seq_len(data$p), function(j) {
    model$residual_variance_inv * data$d[j]
  })

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

#' @keywords internal
check_convergence.mv_individual <- function(data, params, model,
                                             elbo, iter, tracking) {
  if (iter <= 1) {
    model$converged <- FALSE
    return(model)
  }
  delta <- elbo[iter + 1] - elbo[iter]
  if (is.na(delta) || is.infinite(delta)) {
    # Fallback to PIP convergence
    if (!is.null(tracking$convergence$prev_alpha)) {
      pip_diff <- max(abs(tracking$convergence$prev_alpha - model$alpha))
      model$converged <- (pip_diff < params$tol)
    } else {
      model$converged <- FALSE
    }
  } else {
    if (params$verbose) {
      message("ELBO: ", format(elbo[iter + 1], digits = 6))
    }
    model$converged <- (delta < params$tol)
  }
  return(model)
}

# =============================================================================
# EM PRIOR VARIANCE UPDATE
# =============================================================================

#' @keywords internal
em_update_prior_variance.mv_individual <- function(data, params, model,
                                                     alpha, moments, V_init) {
  # Unified EM update using cached stats from mashr calc_sermix_rcpp.
  # prior_scale_em_update is a (K+1)-vector (trace-based scalar per component),
  # V_structure_rank is a K-vector (non-null components only).
  # SER posterior mixture weights determine the contribution of each component.

  if (!is.null(model$em_cache) && !is.null(model$llik_cache)) {
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

  # Fallback: simple EM using marginalized posteriors (for single-matrix K=1)
  mu2 <- Reduce("+", lapply(seq_along(alpha), function(j) {
    alpha[j] * moments$post_mean2[j, , ]
  }))
  scalar <- sum(diag(model$V_structure_inv[, , 1] %*% mu2)) /
            model$V_structure_rank[1]
  return(max(0, scalar))
}

# =============================================================================
# TRACKING
# =============================================================================

#' @keywords internal
track_ibss_fit.mv_individual <- function(data, params, model,
                                          tracking, iter, elbo, ...) {
  # Same as default behavior
  if (iter == 1) {
    tracking$convergence <- list(
      prev_elbo  = -Inf,
      prev_alpha = model$alpha,
      prev_pip_diff = NULL
    )
  } else {
    pip_diff <- max(abs(tracking$convergence$prev_alpha - model$alpha))
    tracking$convergence$prev_elbo  <- elbo[iter]
    tracking$convergence$prev_alpha <- model$alpha
    tracking$convergence$prev_pip_diff <- pip_diff
  }

  if (isTRUE(params$track_fit)) {
    tracking[[iter]] <- list(
      alpha  = model$alpha,
      niter  = iter,
      V      = model$V,
      sigma2 = model$sigma2,
      mu     = model$mu,
      mu2    = model$mu2,
      KL     = model$KL,
      lbf    = model$lbf,
      elbo   = elbo[iter + 1]
    )
  }
  return(tracking)
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
  for (l in seq_len(nrow(model$alpha))) {
    # V is now a numeric scalar per effect
    if (is.na(model$V[l]) || abs(model$V[l]) < params$prior_tol) {
      model$V[l]               <- 0
      model$alpha[l, ]         <- 1 / ncol(model$alpha)
      model$mu[l, , ]          <- 0
      model$mu2[[l]]           <- array(0, dim(model$mu2[[l]]))
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
  b_sum <- compute_posterior_mean_sum_from_model(model)
  # Rescale: coefficients on original scale
  coefs_original <- b_sum / data$csd
  # Intercept = Y_mean - X_mean' %*% coefs
  data$Y_mean - as.vector(t(data$cm) %*% coefs_original)
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
  model$residuals        <- NULL
  model$fitted_without_l <- NULL
  model$raw_residuals    <- NULL
  model$llik_cache       <- NULL
  model$em_cache         <- NULL
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
  b_sum <- compute_posterior_mean_sum_from_model(model)
  Xb <- data$X %*% b_sum
  R_mat <- data$Y - Xb
  # Use effective sample size when Y has missing data
  N <- if (data$Y_has_missing) data$n_obs else data$n
  R <- data$R

  # Zero out missing observations
  if (data$Y_has_missing) {
    R_mat[data$Y_missing] <- 0
    Xb[data$Y_missing] <- 0
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
    if (data$Y_has_missing) XB1_l[data$Y_missing] <- 0
    b1_XtX_b1 <- b1_XtX_b1 + crossprod(XB1_l)  # B1_l' X'X B1_l
  }

  V_est <- (E_RtR + bxxb - b1_XtX_b1) / N
  V_est <- (V_est + t(V_est)) / 2
  return(V_est)
}

# Multivariate ELBO expected log-likelihood (dense)
compute_multivariate_elbo <- function(data, model) {
  # Use effective sample size when Y has missing data
  N <- if (data$Y_has_missing) data$n_obs else data$n
  R <- data$R
  v_inv <- get_v_inv(data, model)

  # Constant: -N*R/2*log(2*pi) - N/2*log|sigma2|
  loglik <- -N * R / 2 * log(2 * pi)
  if (is.matrix(model$sigma2)) {
    loglik <- loglik - N / 2 * log(det(model$sigma2))
  } else {
    loglik <- loglik - N * R / 2 * log(model$sigma2)
  }

  # ESSR = tr(v_inv * (Y-Xr)'(Y-Xr)) - sum_l tr(v_inv * b_l' X'X b_l)
  #        + sum_l sum_j d_j alpha_lj tr(v_inv * mu2_lj)

  # 1. Raw residual RSS: tr(v_inv * (Y-Xr)'(Y-Xr))
  b_sum <- compute_posterior_mean_sum_from_model(model)
  Xb <- data$X %*% b_sum
  R_mat <- data$Y - Xb
  # Zero out missing observations (they don't contribute to likelihood)
  if (data$Y_has_missing) R_mat[data$Y_missing] <- 0
  if (is.matrix(v_inv)) {
    essr <- sum(R_mat * (R_mat %*% v_inv))
  } else {
    essr <- v_inv * sum(R_mat^2)
  }

  # 2. Subtract first moment^2 per effect: -sum_l tr(v_inv * b_l' X'X b_l)
  for (l in seq_len(nrow(model$alpha))) {
    b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
    Xb_l <- data$X %*% b_l  # N x R
    # Zero out missing rows for consistent ESSR decomposition
    if (data$Y_has_missing) Xb_l[data$Y_missing] <- 0
    if (is.matrix(v_inv)) {
      essr <- essr - sum(Xb_l * (Xb_l %*% v_inv))
    } else {
      essr <- essr - v_inv * sum(Xb_l^2)
    }
  }

  # 3. Add vbxxb: use precomputed per-effect values from compute_kl
  #    (d already accounts for missing data in observed-row version)
  essr <- essr + sum(model$vbxxb)

  return(loglik - 0.5 * essr)
}
