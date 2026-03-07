# =============================================================================
# S3 METHODS FOR MULTIVARIATE SUSIE
#
# Methods dispatched on data class (mv_individual, mv_ss) and model class
# (mvsusie) to integrate with susieR's susie_workhorse framework.
# =============================================================================

#' @importFrom susieR susie_get_cs susie_workhorse is_symmetric_matrix
#' @importFrom mvtnorm dmvnorm
#' @importFrom abind abind
NULL

# =============================================================================
# DATA CONFIGURATION
# =============================================================================

#' @keywords internal
get_var_y.mv_individual <- function(data, ...) {
  if (data$R > 1) cov(data$Y) else var(data$Y)
}

#' @keywords internal
get_var_y.mv_ss <- function(data, ...) {
  if (data$R > 1) cov2cor(data$YtY) else data$YtY / (data$n - 1)
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

  # Residual variance for MV is always a matrix — skip scalar checks
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

#' @keywords internal
ibss_initialize.mv_ss <- function(data, params) {
  ibss_initialize.mv_individual(data, params)
}

# =============================================================================
# MODEL INITIALIZATION
# =============================================================================

#' @keywords internal
#' @importFrom ashr compute_lfsr
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
    # Single matrix → K=1 non-null component.
    # Normalize by max diagonal so V_scalar is on a comparable scale to susieR.
    max_diag <- max(abs(diag(prior_var)))
    V_structure <- list(prior_var / max_diag)
    pi_V <- 1.0
    null_weight <- 0
    V_scalar_init <- max_diag
  } else if (class(prior_var)[1] == "mash_prior") {
    # Mixture prior from create_mash_prior.
    # Normalize each component by the same max diagonal as the matrix path
    # so that V_scalar is comparable between matrix and mash_prior paths.
    max_diag <- max(vapply(prior_var$xUlist,
                            function(U) max(abs(diag(U))), numeric(1)))
    V_structure <- lapply(prior_var$xUlist, function(U) U / max_diag)
    pi_V <- prior_var$pi
    null_weight <- prior_var$null_weight
    V_scalar_init <- max_diag
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
    V_structure       = V_structure,       # list of K R×R matrices (non-null)
    V_structure_3d    = matlist2array(V_structure),  # R×R×K array
    V_structure_inv   = V_structure_inv,   # R×R×K pseudo-inverses
    V_structure_rank  = V_structure_rank,  # K-vector of ranks
    null_weight       = null_weight,       # scalar
    pi_V              = pi_V,              # K-vector (non-null weights)
    pi_V_posterior    = vector("list", L), # per-effect J×(K+1) mixture weights
    conditional_lfsr  = vector("list", L), # per-effect J×R LFSR
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

#' @keywords internal
compute_residuals.mv_ss <- function(data, params, model, l, ...) {
  # Fitted for effect l
  b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R

  # X'R = X'Y - XtX %*% b_without_l
  # b_without_l = b_sum - b_l
  b_sum <- compute_posterior_mean_sum_from_model(model)
  XtR <- data$XtY - data$XtX %*% (b_sum - b_l)

  model$residuals <- XtR
  return(model)
}

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
    shat2        = svs,            # list of J R×R matrices
    shat2_inv    = svs_inv,        # list of J R×R matrices
    optim_init   = optim_init,
    optim_bounds = c(-30, 15),
    optim_scale  = "log"
  )
}

#' @keywords internal
compute_ser_statistics.mv_ss <- function(data, params, model, l, ...) {
  compute_ser_statistics.mv_individual(data, params, model, l, ...)
}

# =============================================================================
# LOG BAYES FACTOR & POSTERIOR WEIGHTS
# =============================================================================

#' @keywords internal
loglik.mv_individual <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  R <- data$R
  J <- data$p
  n_thread <- if (!is.null(params$n_thread)) params$n_thread else 1L

  # When V ≈ 0 or NA, the prior is null: lbf = 0 everywhere, uniform alpha
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
  llik <- mashr:::calc_lik_rcpp(
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
  loglik_null <- mashr:::compute_null_loglik_from_matrix(llik_obj, s_alpha)
  loglik_alt <- mashr:::compute_alt_loglik_from_matrix_and_pi(
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
loglik.mv_ss <- function(data, params, model, V, ser_stats, l = NULL, ...) {
  loglik.mv_individual(data, params, model, V, ser_stats, l, ...)
}

#' @keywords internal
neg_loglik.mv_individual <- function(data, params, model, V_param, ser_stats, ...) {
  V <- exp(V_param)
  -loglik.mv_individual(data, params, model, V, ser_stats)
}

#' @keywords internal
neg_loglik.mv_ss <- function(data, params, model, V_param, ser_stats, ...) {
  neg_loglik.mv_individual(data, params, model, V_param, ser_stats, ...)
}

# =============================================================================
# POSTERIOR MOMENTS
# =============================================================================

#' @keywords internal
calculate_posterior_moments.mv_individual <- function(data, params, model,
                                                       V, l, ...) {
  J <- data$p
  R <- data$R

  # When V ≈ 0 or NA, posterior is null: mu = 0, mu2 = 0
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
  post <- mashr:::calc_sermix_rcpp(
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

#' @keywords internal
calculate_posterior_moments.mv_ss <- function(data, params, model, V, l, ...) {
  calculate_posterior_moments.mv_individual(data, params, model, V, l, ...)
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
compute_kl.mv_ss <- function(data, params, model, l) {
  result <- SER_posterior_e_loglik.mv_ss(data, params, model, l)
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

#' @keywords internal
SER_posterior_e_loglik.mv_ss <- function(data, params, model, l) {
  v_inv   <- get_v_inv(data, model)
  alpha_l <- model$alpha[l, ]
  mu_l    <- model$mu[l, , , drop = TRUE]  # J x R
  R <- data$R

  # E1 = tr(v_inv * (B1'XtR + XtR'B1))
  if (is.matrix(v_inv)) {
    term1 <- sum(diag(v_inv %*% t(alpha_l * mu_l) %*% model$residuals))
  } else {
    term1 <- v_inv * sum(alpha_l * mu_l * model$residuals)
  }

  bxxb_l <- matrix(0, R, R)
  vbxxb_l <- 0
  for (j in seq_len(data$p)) {
    mu2_j <- model$mu2[[l]][j, , , drop = FALSE]
    dim(mu2_j) <- c(R, R)
    pb2_j <- alpha_l[j] * mu2_j
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

#' @keywords internal
update_fitted_values.mv_ss <- function(data, params, model, l, ...) {
  b_sum <- compute_posterior_mean_sum_from_model(model)
  model$XtXr <- data$XtX %*% b_sum
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
update_variance_components.mv_ss <- function(data, params, model, ...) {
  list(sigma2 = estimate_residual_variance_mv_ss(data, model))
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
update_model_variance.mv_ss <- function(data, params, model) {
  # Call the SS-specific estimation function directly
  model$sigma2 <- estimate_residual_variance_mv_ss(data, model)

  # Recompute svs/svs_inv and store in MODEL (not data, since data is immutable).
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

#' @keywords internal
update_derived_quantities.mv_ss <- function(data, params, model) {
  return(model)
}

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

# Estimate residual variance for sufficient statistics
#
# Same formula as estimate_residual_variance_mv but using XtX, XtY, YtY.
# R_hat' R_hat = YtY - 2*B_sum'*XtY + B_sum'*XtX*B_sum
# B1_l' X'X B1_l = crossprod(B1_l, XtX %*% B1_l)
estimate_residual_variance_mv_ss <- function(data, model) {
  b_sum <- compute_posterior_mean_sum_from_model(model)
  N <- data$n
  R <- data$R

  # R_hat' R_hat = YtY - B_sum'XtY - XtY'B_sum + B_sum'XtX B_sum
  # NOTE: B_sum'XtY is NOT symmetric in general, so we must symmetrize
  # properly as E2 + t(E2), not 2*E2.
  E2 <- crossprod(b_sum, data$XtY)  # B_sum' X'Y  (R x R)
  E3 <- crossprod(b_sum, data$XtX %*% b_sum)  # B_sum' X'X B_sum  (R x R)
  RtR <- data$YtY - E2 - t(E2) + E3

  # bxxb: use precomputed per-effect values from compute_kl
  bxxb <- Reduce("+", model$bxxb)

  # B1_l' X'X B1_l  (uses full XtX)
  b1_XtX_b1 <- matrix(0, R, R)
  for (l in seq_len(nrow(model$alpha))) {
    B1_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
    if (!is.matrix(B1_l)) B1_l <- matrix(B1_l, ncol = R)
    b1_XtX_b1 <- b1_XtX_b1 + crossprod(B1_l, data$XtX %*% B1_l)
  }

  V_est <- (RtR + bxxb - b1_XtX_b1) / N
  V_est <- (V_est + t(V_est)) / 2
  return(V_est)
}

# =============================================================================
# OBJECTIVE (ELBO)
# =============================================================================

#' @keywords internal
get_objective.mv_individual <- function(data, params, model, ...) {
  compute_multivariate_elbo(data, model) - sum(model$KL, na.rm = TRUE)
}

#' @keywords internal
get_objective.mv_ss <- function(data, params, model, ...) {
  compute_multivariate_elbo_ss(data, model) - sum(model$KL, na.rm = TRUE)
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

# Multivariate ELBO expected log-likelihood (sufficient statistics)
compute_multivariate_elbo_ss <- function(data, model) {
  N <- data$n
  R <- data$R
  v_inv <- get_v_inv(data, model)

  loglik <- -N * R / 2 * log(2 * pi)
  if (is.matrix(model$sigma2)) {
    loglik <- loglik - N / 2 * log(det(model$sigma2))
  } else {
    loglik <- loglik - N * R / 2 * log(model$sigma2)
  }

  b_sum <- compute_posterior_mean_sum_from_model(model)
  # ESSR = tr(v_inv * (YtY - E2 - t(E2) + E3))
  # where E2 = Eb1'XtY, E3 = Eb1'XtX Eb1
  E2 <- crossprod(b_sum, data$XtY)
  E3 <- crossprod(b_sum, data$XtX %*% b_sum)
  essr <- sum(v_inv * (data$YtY - E2 - t(E2) + E3))

  # Add vbxxb: use precomputed per-effect values from compute_kl
  essr <- essr + sum(model$vbxxb)

  # Subtract first-moment squared per effect using FULL XtX (not diagonal)
  # This is: -sum_l tr(v_inv * b_l' X'X b_l) where b_l = alpha_l * mu_l
  for (l in seq_len(nrow(model$alpha))) {
    b_l <- drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]  # J x R
    XtXb_l <- data$XtX %*% b_l  # J x R
    essr <- essr - sum(v_inv * (t(b_l) %*% XtXb_l))
  }

  return(loglik - 0.5 * essr)
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

#' @keywords internal
check_convergence.mv_ss <- function(data, params, model,
                                     elbo, iter, tracking) {
  check_convergence.mv_individual(data, params, model, elbo, iter, tracking)
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

#' @keywords internal
em_update_prior_variance.mv_ss <- function(data, params, model,
                                             alpha, moments, V_init) {
  em_update_prior_variance.mv_individual(data, params, model,
                                           alpha, moments, V_init)
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

#' @keywords internal
track_ibss_fit.mv_ss <- function(data, params, model,
                                   tracking, iter, elbo, ...) {
  track_ibss_fit.mv_individual(data, params, model, tracking, iter, elbo, ...)
}

# =============================================================================
# PRIOR VALIDATION
# =============================================================================

#' @keywords internal
validate_prior.mv_individual <- function(data, params, model, ...) {
  invisible(TRUE)
}

#' @keywords internal
validate_prior.mv_ss <- function(data, params, model, ...) {
  invisible(TRUE)
}

# =============================================================================
# MODEL FIELD ACCESSORS (dispatched on model class "mvsusie")
# =============================================================================

#' @keywords internal
get_prior_variance_l.mvsusie <- function(model, l) {
  # Return the scalar scale. Methods that need the full R×R matrix
  # reconstruct it as V[l] * model$V_structure.
  model$V[l]
}

#' @keywords internal
set_prior_variance_l.mvsusie <- function(model, l, V) {
  model$V[l] <- V
  model
}

#' @keywords internal
get_alpha_l.mvsusie <- function(model, l) {
  model$alpha[l, ]
}

#' @keywords internal
get_posterior_moments_l.mvsusie <- function(model, l) {
  list(
    post_mean  = model$mu[l, , , drop = TRUE],    # J x R
    post_mean2 = model$mu2[[l]]                     # J x R x R
  )
}

#' @keywords internal
get_posterior_mean_l.mvsusie <- function(model, l) {
  drop(model$alpha[l, ]) * model$mu[l, , , drop = TRUE]
}

#' @keywords internal
get_posterior_mean_sum.mvsusie <- function(model) {
  compute_posterior_mean_sum_from_model(model)
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

#' @keywords internal
trim_null_effects.mv_ss <- function(data, params, model) {
  trim_null_effects.mv_individual(data, params, model)
}

# =============================================================================
# OUTPUT / POST-PROCESSING
# =============================================================================

#' @keywords internal
get_scale_factors.mv_individual <- function(data, params, ...) {
  data$csd
}

#' @keywords internal
get_scale_factors.mv_ss <- function(data, params, ...) {
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

#' @keywords internal
get_intercept.mv_ss <- function(data, params, model, ...) {
  if (!is.null(data$X_colmeans) && !is.null(data$Y_colmeans)) {
    b_sum <- compute_posterior_mean_sum_from_model(model)
    coefs_original <- b_sum / data$csd
    data$Y_colmeans - as.vector(t(data$X_colmeans) %*% coefs_original)
  } else {
    rep(NA_real_, data$R)
  }
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
get_fitted.mv_ss <- function(data, params, model, ...) {
  # For SS, fitted is XtX %*% b_sum (on standardized scale)
  b_sum <- compute_posterior_mean_sum_from_model(model)
  data$XtX %*% b_sum
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
get_cs.mv_ss <- function(data, params, model, ...) {
  if (is.null(params$coverage) || is.null(params$min_abs_corr))
    return(NULL)
  # Use correlation from XtX
  Xcorr <- tryCatch(cov2cor(data$XtX), error = function(e) NULL)
  if (is.null(Xcorr)) return(NULL)
  susie_get_cs(model,
    Xcorr        = Xcorr,
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
get_variable_names.mv_ss <- function(data, model, ...) {
  vnames <- colnames(data$XtX)
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
get_zscore.mv_ss <- function(data, params, model, ...) {
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

#' @keywords internal
cleanup_model.mv_ss <- function(data, params, model, ...) {
  model$residuals  <- NULL
  model$llik_cache <- NULL
  model$em_cache   <- NULL
  return(model)
}

# =============================================================================
# MULTIVARIATE MATH UTILITIES
# =============================================================================

# Log Bayes factor for each variable, given matrix prior U
multivariate_lbf <- function(betahat, S, U) {
  # Log Bayes factor per variable: log p(bhat|H1) - log p(bhat|H0)
  # Using dmvnorm for numerical stability
  lbf <- sapply(seq_along(S), function(j) {
    mvtnorm::dmvnorm(x = betahat[j, ], sigma = S[[j]] + U, log = TRUE) -
      mvtnorm::dmvnorm(x = betahat[j, ], sigma = S[[j]], log = TRUE)
  })
  lbf[is.nan(lbf)] <- 0
  return(lbf)
}
