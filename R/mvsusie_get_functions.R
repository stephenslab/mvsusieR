# Output formatting and post-hoc extraction functions.
#
# Functions for converting internal susie_workhorse output to the
# standard mvsusie output format, applying dimension names, and
# computing derived quantities like log Bayes factors and LFSR.

#' @importFrom mvtnorm dmvnorm
NULL

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
