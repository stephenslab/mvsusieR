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
    # mu2_cache[[l]]$mu2_diag is J x R (diagonal of E[b^2])
    cache_l <- s$mu2_cache[[l]]
    if (!is.null(cache_l) && !is.null(cache_l$mu2_diag)) {
      b2[l, , ] <- s$alpha[l, ] * cache_l$mu2_diag
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
  # cm and csd can be J-vectors (standard path) or J x R matrices
  # (missing data 3d path with per-outcome centering/scaling).
  # colSums(cm * coefs) handles both: J-vector recycles per column,
  # J x R does element-wise multiplication.
  coefs_original <- b_sum / csd  # J x R
  intercept_vec <- Y_mean - colSums(cm * coefs_original)
  coef <- rbind(matrix(intercept_vec, 1, R), coefs_original)

  # ----- b1_rescaled: per-effect rescaled b1 (L x (J+1) x R) -----
  b1_rescaled <- array(0, c(L, J + 1, R))
  for (l in seq_len(L)) {
    b1_l <- matrix(b1[l, , ], J, R) / csd  # J x R (unscale)
    intercept_l <- Y_mean - colSums(cm * b1_l)
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

  # Forward track_fit trace if present
  if (!is.null(s$trace)) out$trace <- s$trace

  # Preserve V if estimated.
  # For single-matrix prior (K=1), convert V from scalar multiplier to
  # absolute prior variance: V_abs = V_scalar * V_structure[[1]][1,1].
  # This matches susieR's output convention where V is on the absolute scale.
  # For mixture priors (K>1) and multivariate (R>1), V stays as the scalar
  # multiplier since there is no single "absolute" prior variance.
  if (estimate_prior_variance) {
    out$V <- s$V
    if (length(s$V_structure) == 1 && is.matrix(s$V_structure[[1]])) {
      out$V <- s$V * s$V_structure[[1]][1, 1]
    }
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
  cnames <- s$outcome_names
  lnames <- paste0("l", seq_len(L))

  # alpha: L x J
  dimnames(s$alpha) <- list(lnames, vnames)

  # lbf_variable: L x J
  if (!is.null(s$lbf_variable))
    dimnames(s$lbf_variable) <- list(lnames, vnames)

  # lbf: L-vector
  if (!is.null(s$lbf))
    names(s$lbf) <- lnames

  # pip: J-vector
  if (!is.null(s$pip))
    names(s$pip) <- vnames

  # sigma2: R x R matrix
  if (!is.null(s$sigma2) && is.matrix(s$sigma2))
    dimnames(s$sigma2) <- list(cnames, cnames)

  # V_structure: list of R x R matrices
  if (!is.null(s$V_structure) && is.list(s$V_structure)) {
    for (k in seq_along(s$V_structure)) {
      if (is.matrix(s$V_structure[[k]]))
        dimnames(s$V_structure[[k]]) <- list(cnames, cnames)
    }
  }

  R <- length(cnames)
  if (R > 1) {
    # b1, b2: L x J x R arrays
    if (!is.null(s$b1) && length(dim(s$b1)) == 3)
      dimnames(s$b1) <- list(lnames, vnames, cnames)
    if (!is.null(s$b2) && length(dim(s$b2)) == 3)
      dimnames(s$b2) <- list(lnames, vnames, cnames)

    # lfsr: J x R matrix
    if (!is.null(s$lfsr) && is.matrix(s$lfsr))
      dimnames(s$lfsr) <- list(vnames, cnames)

    # single_effect_lfsr: L x R matrix
    if (!is.null(s$single_effect_lfsr) && is.matrix(s$single_effect_lfsr))
      dimnames(s$single_effect_lfsr) <- list(lnames, cnames)

    # coef: (J+1) x R
    dimnames(s$coef) <- list(c("(Intercept)", vnames), cnames)

    # fitted: N x R or J x R
    if (!is.null(dim(s$fitted)) && ncol(s$fitted) == R)
      colnames(s$fitted) <- cnames
  } else {
    # R=1: coef is J+1 vector
    names(s$coef) <- c("(Intercept)", vnames)
  }

  return(s)
}

# Log Bayes factor for each variable, given matrix prior U
multivariate_lbf <- function(betahat, S, U) {
  # Log Bayes factor per variable: log p(bhat|H1) - log p(bhat|H0)
  # Using dmvnorm for numerical stability
  J <- nrow(betahat)
  lbf <- sapply(seq_len(J), function(j) {
    S_j <- S[[min(j, length(S))]]
    mvtnorm::dmvnorm(x = betahat[j, ], sigma = S_j + U, log = TRUE) -
      mvtnorm::dmvnorm(x = betahat[j, ], sigma = S_j, log = TRUE)
  })
  lbf[is.nan(lbf)] <- 0
  return(lbf)
}

#' @title Local false sign rate (lfsr) for single effects
#'
#' @details This function returns the lfsr for identifying nonzero
#'   single effects, separately for each outcome.
#'
#' @param alpha L x J matrix.
#'
#' @param clfsr L x J x R conditonal lfsr.
#'
#' @return L x R matrix of lfsr
#'
#' @export
#'
mvsusie_single_effect_lfsr <- function(clfsr, alpha) {
  if (!is.array(clfsr) && is.na(clfsr)) {
    return(as.numeric(NA))
  } else {
    return(do.call(
      cbind,
      lapply(
        1:dim(clfsr)[3],
        function(r) {
          clfsrr <- clfsr[, , r]
          if (is.null(nrow(clfsrr))) {
            clfsrr <- matrix(clfsrr, 1, length(clfsrr))
          }
          return(pmax(0, rowSums(alpha * clfsrr)))
        }
      )
    ))
  }
}

#' @title Local false sign rate (lfsr) for variables.
#'
#' @details This function returns the lfsr for identifying nonzero
#'   effects for each outcome.
#'
#' @param alpha L x J matrix.
#'
#' @param clfsr L x J x R conditonal lfsr.
#'
#' @param weighted Set \code{weighted = TRUE} to weight lfsr by PIP;
#'   otherwise set \code{weighted = FALSE}.
#'
#' @return J x R lfsr matrix.
#'
#' @export
#'
mvsusie_get_lfsr <- function(clfsr, alpha, weighted = TRUE) {
  if (!is.array(clfsr) && is.na(clfsr)) {
    return(as.numeric(NA))
  } else {
    if (!weighted) {
      alpha <- matrix(1, nrow(alpha), ncol(alpha))
    }
    return(do.call(
      cbind,
      lapply(
        1:dim(clfsr)[3],
        function(r) {
          true_sign_mat <- alpha * (1 - clfsr[, , r])
          pmax(1e-20, 1 - apply(true_sign_mat, 2, max))
        }
      )
    ))
  }
}

#' @title Outcome-specific posterior inclusion probabilities.
#'
#' @description Computes the PIP for each variable in each outcome,
#'   accounting for the mixture structure of the prior. A variable is
#'   considered "included" in outcome r if its effect has nonzero
#'   variance in that outcome under the selected mixture component.
#'
#' @param m A fitted mvsusie object (output of \code{mvsusie}).
#' @param prior_obj Optional mash prior object. If not provided, the
#'   active prior matrices are taken from \code{m$V_structure}.
#'
#' @return J x R matrix of outcome-specific PIPs.
#'
#' @export
#'
mvsusie_get_pip_per_outcome <- function(m, prior_obj = NULL) {
  alpha_out <- mvsusie_get_alpha_per_outcome(m, prior_obj)
  R <- dim(alpha_out)[3]
  do.call(cbind, lapply(
    seq_len(R),
    function(r) apply(alpha_out[, , r], 2, function(x) 1 - prod(1 - x))
  ))
}

#' @title Outcome-specific alpha (per-effect inclusion weights).
#'
#' @description For each single effect l and variable j, computes the
#'   probability of having a nonzero effect in each outcome r.
#'   This is the sum of mixture_weights over components that have
#'   nonzero prior variance in outcome r, multiplied by alpha.
#'
#' @param m A fitted mvsusie object (output of \code{mvsusie}).
#' @param prior_obj Optional mash prior object. If not provided, the
#'   active prior matrices are taken from \code{m$V_structure}.
#'
#' @return L x J x R array of outcome-specific alpha values.
#'
#' @keywords internal
mvsusie_get_alpha_per_outcome <- function(m, prior_obj = NULL) {
  # Build outcome indicator: which components have nonzero variance
  # in which outcomes.
  #
  # Use m$V_structure (the *active* prior matrices after estimation/pruning)
  # rather than prior_obj$xUlist, because the model may have pruned

  # components and m$mixture_weights matches V_structure, not the
  # original prior.
  active_matrices <- m$V_structure
  if (is.null(active_matrices) && !is.null(prior_obj)) {
    active_matrices <- prior_obj$xUlist
  }
  if (is.null(active_matrices)) {
    stop("Cannot determine active prior matrices. ",
         "The model must have V_structure or a prior_obj must be provided.")
  }
  outcome_indicator <- do.call(
    rbind,
    lapply(
      seq_along(active_matrices),
      function(i) as.integer(diag(active_matrices[[i]]) != 0)
    )
  )  # K x R

  # m$mixture_weights is L x J x (K+1) where column 1 is the null component.
  # We need the K non-null columns (columns 2:(K+1)).
  L <- nrow(m$alpha)
  J <- ncol(m$alpha)
  R <- ncol(outcome_indicator)

  alpha_out <- array(0, c(L, J, R))
  for (r in seq_len(R)) {
    for (k in seq_len(nrow(outcome_indicator))) {
      # k-th non-null component = column k+1 in mixture_weights
      alpha_out[, , r] <- alpha_out[, , r] +
        m$mixture_weights[, , k + 1] * outcome_indicator[k, r]
    }
    alpha_out[, , r] <- alpha_out[, , r] * m$alpha
  }
  alpha_out
}
