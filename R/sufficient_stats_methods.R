#' @importFrom susieR susie_get_cs
NULL

#' @keywords internal
get_var_y.mv_ss <- function(data, ...) {
  if (data$R > 1) cov2cor(data$YtY) else data$YtY / (data$n - 1)
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

#' @keywords internal
SER_posterior_e_loglik.mv_ss <- function(data, params, model, l) {
  v_inv   <- get_v_inv(data, model)
  alpha_l <- model$alpha[l, ]
  mu_l    <- model$mu[l, , , drop = TRUE]  # J x R
  R <- data$R

  # E1 = tr(v_inv * (B1'XtR + XtR'B1))
  term1 <- sum(diag(v_inv %*% t(alpha_l * mu_l) %*% model$residuals))

  # Use precomputed bxxb/vbxxb from calculate_posterior_moments
  cache <- model$mu2_cache[[l]]
  bxxb_l  <- cache$bxxb
  vbxxb_l <- cache$vbxxb

  eloglik <- 0.5 * (2 * term1 - vbxxb_l)
  return(list(eloglik = eloglik, bxxb = bxxb_l, vbxxb = vbxxb_l))
}

#' @keywords internal
update_fitted_values.mv_ss <- function(data, params, model, l, ...) {
  b_sum <- compute_posterior_mean_sum_from_model(model)
  model$XtXr <- data$XtX %*% b_sum
  return(model)
}

#' @keywords internal
update_variance_components.mv_ss <- function(data, params, model, ...) {
  list(sigma2 = estimate_residual_variance_mv_ss(data, model))
}

#' @keywords internal
update_model_variance.mv_ss <- function(data, params, model) {
  # Update residual variance if requested
  if (isTRUE(params$estimate_residual_variance)) {
    model$sigma2 <- estimate_residual_variance_mv_ss(data, model)

    # Recompute svs/svs_inv and store in MODEL (not data, since data is immutable).
    model$residual_variance_inv <- invert_via_chol(model$sigma2)$inv
    model$residual_correlation <- cov2cor(model$sigma2)
    # For common-cov, store single copy to save O(J*R^2) memory
    if (data$is_common_cov) {
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
        model$svs, model$V_structure, data$is_common_cov)
  }

  # Update mixture prior weights and prune near-zero components
  model <- update_mixture_weights_and_prune(model, params)

  return(model)
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
  # Ensure positive-definiteness.  When YtY is approximate (e.g. the

  # z-score path where cor(Y) is estimated from null z-scores), the
  # formula can produce a non-PD result.  Add a small ridge as a
  # safety net, mirroring the individual-data imputation path.
  if (!is_pd(V_est)) V_est <- makePD(V_est)
  return(V_est)
}

#' @keywords internal
get_objective.mv_ss <- function(data, params, model, ...) {
  compute_multivariate_elbo_ss(data, model) - sum(model$KL, na.rm = TRUE)
}

# Multivariate ELBO expected log-likelihood (sufficient statistics)
compute_multivariate_elbo_ss <- function(data, model) {
  N <- data$n
  R <- data$R
  v_inv <- get_v_inv(data, model)

  loglik <- -N * R / 2 * log(2 * pi)
  loglik <- loglik - N / 2 * log_det_sym(model$sigma2)

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

#' @keywords internal
get_intercept.mv_ss <- function(data, params, model, ...) {
  if (!is.null(data$X_colmeans) && !is.null(data$Y_colmeans)) {
    b_sum <- compute_posterior_mean_sum_from_model(model)
    coefs_original <- b_sum / data$csd
    data$Y_colmeans - as.vector(t(data$X_colmeans) %*% coefs_original)
  } else {
    # Without column means, the intercept is 0 on the centered scale
    # (the model operates on centered sufficient statistics).
    rep(0, data$R)
  }
}

#' Get fitted values on the original Y scale.
#'
#' @keywords internal
get_fitted.mv_ss <- function(data, params, model, ...) {
  # For SS, fitted is XtX %*% b_sum (on standardized scale)
  b_sum <- compute_posterior_mean_sum_from_model(model)
  data$XtX %*% b_sum
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
get_variable_names.mv_ss <- function(data, model, ...) {
  vnames <- colnames(data$XtX)
  if (!is.null(vnames)) {
    colnames(model$alpha)        <- vnames
    colnames(model$lbf_variable) <- vnames
  }
  return(model)
}

