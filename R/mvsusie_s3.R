# =============================================================================
# S3 ENTRY POINTS FOR MULTIVARIATE SUSIE
#
# These are the S3-based entry points for mvsusie() and mvsusie_suff_stat().
# They handle prior scaling, data construction, and call mvsusie_workhorse.
# =============================================================================

# =============================================================================
# OUTPUT FORMAT CONVERSION
#
# Converts S3 workhorse output (mu, mu2, etc.) to the standard mvsusie output
# format expected by tests and downstream users (b1, b2, coef, etc.).
# =============================================================================

#' Convert S3 workhorse output to mvsusie output format
#'
#' @param s Model output from susie_workhorse (with class "mvsusie"/"susie")
#' @param csd Column standard deviations used for standardization
#' @param cm Column means used for centering
#' @param Y_mean Y column means (for intercept recovery)
#' @param estimate_prior_variance Whether prior variance was estimated
#' @param is_ss Whether this is a sufficient statistics fit
#'
#' @return Model in standard mvsusie output format
#'
#' @keywords internal
format_mvsusie_output <- function(s, csd, cm, Y_mean,
                                  estimate_prior_variance = TRUE,
                                  is_ss = FALSE) {
  L <- nrow(s$alpha)
  J <- ncol(s$alpha)
  R <- dim(s$mu)[3]

  # ----- b1: alpha-weighted first moment (L x J x R) -----
  b1 <- array(0, c(L, J, R))
  for (l in seq_len(L)) {
    b1[l, , ] <- drop(s$alpha[l, ]) * s$mu[l, , , drop = TRUE]
  }

  # ----- b2: alpha-weighted diag of second moment (L x J x R) -----
  b2 <- array(0, c(L, J, R))
  for (l in seq_len(L)) {
    # mu2[[l]] is J x R x R; extract diagonal for each variable
    for (j in seq_len(J)) {
      mu2_j <- s$mu2[[l]][j, , , drop = FALSE]
      dim(mu2_j) <- c(R, R)
      b2[l, j, ] <- s$alpha[l, j] * diag(mu2_j)
    }
  }

  # ----- b_sum: total posterior mean (J x R) -----
  b_sum <- matrix(0, J, R)
  for (l in seq_len(L)) {
    for (r in seq_len(R)) {
      b_sum[, r] <- b_sum[, r] + b1[l, , r]
    }
  }

  # ----- coef: rescaled coefficients with intercept row -----
  coefs_original <- b_sum / csd  # J x R
  intercept_vec <- Y_mean - as.vector(t(cm) %*% coefs_original)
  coef <- rbind(matrix(intercept_vec, 1, R), coefs_original)

  # ----- b1_rescaled: per-effect rescaled b1 (L x (J+1) x R) -----
  b1_rescaled <- array(0, c(L, J + 1, R))
  for (l in seq_len(L)) {
    b1_l <- matrix(b1[l, , ], J, R) / csd  # J x R (unscale)
    intercept_l <- Y_mean - as.vector(t(cm) %*% b1_l)
    b1_rescaled[l, , ] <- rbind(matrix(intercept_l, 1, R), b1_l)
  }

  # ----- Collapse arrays for R = 1 -----
  if (R == 1) {
    b1 <- b1[, , 1]           # L x J
    b2 <- b2[, , 1]           # L x J
    b_sum <- as.vector(b_sum) # J
    coef <- as.vector(coef)   # J+1
    b1_rescaled <- drop(b1_rescaled)  # L x (J+1)
  }

  # ----- Convergence status -----
  convergence <- list(
    converged = isTRUE(s$converged),
    niter     = s$niter
  )

  # ----- Build output list -----
  # Keep mu/mu2 for susieR coef.susie() compatibility (uses alpha * mu)
  # For R=1, mu needs to be L x J matrix (not L x J x 1 array)
  if (R == 1) {
    mu_out <- s$mu[, , 1, drop = TRUE]  # L x J
    if (!is.matrix(mu_out)) mu_out <- matrix(mu_out, L, J)
  } else {
    mu_out <- s$mu
  }

  # ----- Mixture weights and LFSR -----
  has_mixture <- !is.null(s$pi_V_posterior) &&
                 !is.null(s$pi_V_posterior[[1]])

  if (has_mixture) {
    # mixture_weights: L x J x (K+1) array
    K_plus_1 <- ncol(s$pi_V_posterior[[1]])
    mixture_weights <- array(NA, c(L, J, K_plus_1))
    for (l in seq_len(L)) {
      if (!is.null(s$pi_V_posterior[[l]])) {
        mixture_weights[l, , ] <- s$pi_V_posterior[[l]]
      }
    }

    # conditional_lfsr: L x J x R array
    # Trimmed effects (NULL clfsr) get lfsr = 1 (no signal)
    clfsr <- array(1, c(L, J, R))
    for (l in seq_len(L)) {
      if (!is.null(s$conditional_lfsr[[l]])) {
        clfsr[l, , ] <- s$conditional_lfsr[[l]]
      }
    }

    lfsr <- mvsusie_get_lfsr(clfsr, s$alpha)
    single_effect_lfsr <- mvsusie_single_effect_lfsr(clfsr, s$alpha)
  } else {
    mixture_weights <- as.numeric(NA)
    clfsr <- as.numeric(NA)
    lfsr <- as.numeric(NA)
    single_effect_lfsr <- as.numeric(NA)
  }

  out <- list(
    alpha            = s$alpha,
    mu               = mu_out,
    b1               = b1,
    b2               = b2,
    KL               = s$KL,
    lbf              = s$lbf,
    lbf_variable     = s$lbf_variable,
    V                = s$V,
    V_structure      = s$V_structure,
    null_weight      = s$null_weight,
    pi_V             = s$pi_V,
    sigma2           = s[["sigma2"]],
    elbo             = s$elbo,
    niter            = s$niter,
    convergence      = convergence,
    coef             = coef,
    fitted           = s$fitted,
    intercept        = s$intercept,
    b1_rescaled      = b1_rescaled,
    X_column_scale_factors = s$X_column_scale_factors,
    mixture_weights  = mixture_weights,
    conditional_lfsr = clfsr,
    lfsr             = lfsr,
    single_effect_lfsr = single_effect_lfsr,
    sets             = s$sets,
    pip              = s$pip
  )

  # Preserve V if estimated
  if (estimate_prior_variance) {
    out$V <- s$V
  }

  class(out) <- "mvsusie"
  return(out)
}


#' Apply dimnames to mvsusie output to match standard format
#'
#' @param s Model output from format_mvsusie_output
#' @keywords internal
apply_mvsusie_dimnames <- function(s) {
  L <- nrow(s$alpha)
  J <- ncol(s$alpha)
  vnames <- s$variable_names
  cnames <- s$condition_names
  lnames <- paste0("l", seq_len(L))

  # alpha: L x J with row=effect names, col=variable names  dimnames(s$alpha) <- list(lnames, vnames)

  # lbf_variable: no dimnames  dimnames(s$lbf_variable) <- NULL

  # lbf: no names  names(s$lbf) <- NULL

  R <- length(cnames)
  if (R > 1) {
    # b1, b2: no dimnames    dimnames(s$b1) <- NULL
    dimnames(s$b2) <- NULL

    # coef: (J+1) x R
    coef_rownames <- c("(Intercept)", vnames)
    dimnames(s$coef) <- list(coef_rownames, cnames)

    # fitted: N x R or J x R
    if (!is.null(dim(s$fitted)) && ncol(s$fitted) == R) {
      colnames(s$fitted) <- cnames
    }
  } else {
    # R=1: b1 is L x J, coef is J+1 vector
    names(s$coef) <- c("(Intercept)", vnames)
  }

  return(s)
}

#' @rdname mvsusie
#'
#' @importFrom stats sd var cov cov2cor
#' @importFrom susieR susie_get_cs
#'
#' @export
#'
mvsusie_s3 <- function(X, Y, L = 10, prior_variance = 0.2,
                       residual_variance = NULL, prior_weights = NULL,
                       standardize = TRUE, intercept = TRUE,
                       approximate = FALSE,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       check_null_threshold = 0, prior_tol = 1e-9,
                       compute_objective = TRUE, model_init = NULL,
                       coverage = 0.95, min_abs_corr = 0.5,
                       compute_univariate_zscore = FALSE,
                       n_thread = 1,
                       max_iter = 100, tol = 1e-3, verbosity = 2,
                       track_fit = FALSE) {
  start_time <- proc.time()

  # Validate inputs
  if (is.null(dim(Y))) {
    Y <- matrix(Y, length(Y), 1)
  } else if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  R <- ncol(Y)

  # Adjust prior weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / ncol(X), ncol(X))
  } else {
    prior_weights <- prior_weights / sum(prior_weights)
  }

  # Check prior variance
  is_mash_prior <- class(prior_variance)[1] == "mash_prior"
  is_numeric_prior <- !is.matrix(prior_variance) && !is_mash_prior
  if (R > 1 && is_numeric_prior) {
    stop("Please specify prior variance for the multivariate response Y")
  }

  # Scale prior variance when standardizing
  if (standardize && (is.matrix(prior_variance) || is_mash_prior)) {
    sigma <- sapply(seq_len(R), function(i) sd(Y[, i], na.rm = TRUE))
    n <- sapply(seq_len(R), function(i) length(which(!is.na(Y[, 1]))))
    sigma <- sigma / sqrt(n)
    if (estimate_prior_variance) {
      sigma <- sigma / max(sigma)
    }
    if (is.matrix(prior_variance)) {
      prior_variance <- scale_covariance(prior_variance, sigma)
    } else {
      prior_variance <- scale_prior_variance.mash_prior(prior_variance, sigma)
    }
  }

  # For R=1 scalar prior, convert to 1x1 matrix for uniform MV treatment
  if (R == 1 && is_numeric_prior) {
    prior_variance <- matrix(prior_variance, 1, 1)
  }

  if (verbosity > 1) {
    message("Initializing data object...")
    message(paste("Dimension of X matrix:", nrow(X), ncol(X)))
    message(paste("Dimension of Y matrix:", nrow(Y), ncol(Y)))
  }

  # Create data object
  data <- create_mvsusie_data(X, Y, center = intercept, scale = standardize)

  # Set residual variance
  data <- set_mvsusie_residual_variance(data, residual_variance)

  # Compute prior inverses for EM (mixture priors)
  if (is_mash_prior && estimate_prior_variance &&
      estimate_prior_method == "EM") {
    prior_variance <- compute_prior_inv.mash_prior(prior_variance)
  }

  # Fit model
  s <- mvsusie_workhorse(data, L = L, prior_variance = prior_variance,
                          prior_weights = prior_weights,
                          estimate_residual_variance = estimate_residual_variance,
                          estimate_prior_variance = estimate_prior_variance,
                          estimate_prior_method = estimate_prior_method,
                          check_null_threshold = check_null_threshold,
                          compute_objective = compute_objective,
                          max_iter = max_iter, tol = tol,
                          prior_tol = prior_tol,
                          track_fit = track_fit,
                          verbosity = verbosity,
                          coverage = coverage,
                          min_abs_corr = min_abs_corr,
                          n_thread = n_thread,
                          model_init = model_init)

  # Compute CSs using original X
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets <- susie_get_cs(s, coverage = coverage, X = X,
                           min_abs_corr = min_abs_corr)
  }

  # Store residual variance: use model's updated value if estimated,
  # otherwise fall back to data's initial value
  if (!estimate_residual_variance) {
    s$residual_variance <- data$residual_variance
  }

  # Convert to standard mvsusie output format
  s <- format_mvsusie_output(s, csd = data$csd, cm = data$cm,
                              Y_mean = data$Y_mean,
                              estimate_prior_variance = estimate_prior_variance,
                              is_ss = FALSE)

  # Report z-scores from univariate regression
  if (compute_univariate_zscore) {
    s$z <- calc_z(X, Y, center = intercept, scale = standardize)
  }

  # Set names
  if (is.null(colnames(Y))) {
    s$condition_names <- paste0("cond", seq_len(R))
  } else {
    s$condition_names <- colnames(Y)
  }
  if (is.null(colnames(X))) {
    s$variable_names <- paste0("var", seq_len(ncol(X)))
  } else {
    s$variable_names <- colnames(X)
  }

  # Apply dimnames to match standard output format
  s <- apply_mvsusie_dimnames(s)

  s$walltime <- proc.time() - start_time
  return(s)
}


#' @rdname mvsusie
#'
#' @importFrom stats cov2cor
#' @importFrom susieR susie_get_cs
#'
#' @export
#'
mvsusie_suff_stat_s3 <- function(XtX, XtY, YtY, N, L = 10,
                                  X_colmeans = NULL, Y_colmeans = NULL,
                                  prior_variance = 0.2,
                                  residual_variance = NULL,
                                  prior_weights = NULL,
                                  standardize = TRUE,
                                  estimate_residual_variance = FALSE,
                                  estimate_prior_variance = TRUE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0, prior_tol = 1e-9,
                                  compute_objective = TRUE, model_init = NULL,
                                  coverage = 0.95, min_abs_corr = 0.5,
                                  n_thread = 1,
                                  max_iter = 100, tol = 1e-3, verbosity = 2,
                                  track_fit = FALSE) {
  start_time <- proc.time()

  XtY <- as.matrix(XtY)
  YtY <- as.matrix(YtY)
  R <- ncol(XtY)
  J <- ncol(XtX)

  # Adjust prior weights
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / J, J)
  } else {
    prior_weights <- prior_weights / sum(prior_weights)
  }

  # Check prior variance
  is_mash_prior <- class(prior_variance)[1] == "mash_prior"
  is_numeric_prior <- !is.matrix(prior_variance) && !is_mash_prior
  if (R > 1 && is_numeric_prior) {
    stop("Please specify prior variance for the multivariate response Y")
  }

  # Scale prior variance when standardizing
  if (standardize && (is.matrix(prior_variance) || is_mash_prior)) {
    sigma <- sqrt(diag(YtY) / (N - 1))
    sigma <- sigma / sqrt(N)
    if (estimate_prior_variance) {
      sigma <- sigma / max(sigma)
    }
    if (is.matrix(prior_variance)) {
      prior_variance <- scale_covariance(prior_variance, sigma)
    } else {
      prior_variance <- scale_prior_variance.mash_prior(prior_variance, sigma)
    }
  }

  # For R=1 scalar prior, convert to 1x1 matrix for uniform MV treatment
  if (R == 1 && is_numeric_prior) {
    prior_variance <- matrix(prior_variance, 1, 1)
  }

  # Compute prior inverses for EM (mixture priors)
  if (is_mash_prior && estimate_prior_variance &&
      estimate_prior_method == "EM") {
    prior_variance <- compute_prior_inv.mash_prior(prior_variance)
  }

  # Default residual variance: correlation matrix for R > 1
  if (is.null(residual_variance)) {
    if (R > 1) {
      residual_variance <- cov2cor(YtY)
    } else {
      residual_variance <- as.numeric(YtY / (N - 1))
    }
  }

  # Create SS data object
  data <- create_mvsusie_ss_data(XtX, XtY, YtY, N,
                                  X_colmeans = X_colmeans,
                                  Y_colmeans = Y_colmeans,
                                  standardize = standardize,
                                  residual_variance = residual_variance)

  # Fit model
  s <- mvsusie_workhorse(data, L = L, prior_variance = prior_variance,
                          prior_weights = prior_weights,
                          estimate_residual_variance = estimate_residual_variance,
                          estimate_prior_variance = estimate_prior_variance,
                          estimate_prior_method = estimate_prior_method,
                          check_null_threshold = check_null_threshold,
                          compute_objective = compute_objective,
                          max_iter = max_iter, tol = tol,
                          prior_tol = prior_tol,
                          track_fit = track_fit,
                          verbosity = verbosity,
                          coverage = coverage,
                          min_abs_corr = min_abs_corr,
                          n_thread = n_thread,
                          model_init = model_init)

  # Compute CSs using XtX correlation
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets <- susie_get_cs(s, coverage = coverage,
                           Xcorr = cov2cor(XtX),
                           min_abs_corr = min_abs_corr)
  }

  # Store residual variance: use model's updated value if estimated,
  # otherwise fall back to data's initial value
  if (!estimate_residual_variance) {
    s$residual_variance <- data$residual_variance
  }

  # For SS, use colmeans if available; otherwise use zero vectors
  ss_Y_mean <- if (!is.null(Y_colmeans)) Y_colmeans else rep(0, R)
  ss_cm     <- if (!is.null(X_colmeans)) X_colmeans else rep(0, J)

  # Convert to standard mvsusie output format
  s <- format_mvsusie_output(s, csd = data$csd, cm = ss_cm,
                              Y_mean = ss_Y_mean,
                              estimate_prior_variance = estimate_prior_variance,
                              is_ss = TRUE)

  # Set names
  if (is.null(colnames(XtY))) {
    s$condition_names <- paste0("cond", seq_len(R))
  } else {
    s$condition_names <- colnames(XtY)
  }
  if (is.null(colnames(XtX))) {
    s$variable_names <- paste0("var", seq_len(J))
  } else {
    s$variable_names <- colnames(XtX)
  }

  # Apply dimnames to match standard output format
  s <- apply_mvsusie_dimnames(s)

  s$walltime <- proc.time() - start_time
  return(s)
}
