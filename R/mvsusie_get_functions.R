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

  # ----- mu2_diag: diagonal of posterior second moment (L x J x R) -----
  mu2_diag <- array(0, c(L, J, R))
  for (l in seq_len(L)) {
    cache_l <- s$mu2_cache[[l]]
    if (!is.null(cache_l) && !is.null(cache_l$mu2_diag)) {
      mu2_diag[l, , ] <- cache_l$mu2_diag
    }
  }

  # ----- Collapse arrays for R = 1 -----
  if (R == 1) {
    mu2_diag <- mu2_diag[, , 1]          # L x J
  }

  # ----- Build output list -----
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
    # posterior_mixture_weights: L x J x K array (or L x J x (K+1) if null_weight > 0)
    K_plus_1 <- ncol(s$pi_V_posterior[[1]])
    raw_weights <- array(NA, c(L, J, K_plus_1))
    for (l in seq_len(L)) {
      if (!is.null(s$pi_V_posterior[[l]])) {
        raw_weights[l, , ] <- s$pi_V_posterior[[l]]
      }
    }
    # Drop null column when null_weight == 0 (it's all zeros)
    has_null <- !is.null(s$null_weight) && s$null_weight > 0
    if (has_null) {
      posterior_mixture_weights <- raw_weights
    } else {
      posterior_mixture_weights <- raw_weights[, , -1, drop = FALSE]
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
    posterior_mixture_weights <- as.numeric(NA)
    clfsr <- as.numeric(NA)
    lfsr <- as.numeric(NA)
    single_effect_lfsr <- as.numeric(NA)
  }

  out <- list(
    alpha            = s$alpha,
    mu               = mu_out,
    mu2_diag         = mu2_diag,
    KL               = s$KL,
    lbf              = s$lbf,
    lbf_variable     = s$lbf_variable,
    V                = s$V,
    V_structure      = s$V_structure,
    null_weight      = s$null_weight,
    pi               = s$pi,
    prior_mixture_weights = s$pi_V,
    sigma2           = s[["sigma2"]],
    Xr               = s$Xr,
    elbo             = s$elbo,
    niter            = s$niter,
    converged        = isTRUE(s$converged),
    fitted           = s$fitted,
    intercept        = s$intercept,
    X_column_scale_factors = csd,
    posterior_mixture_weights = posterior_mixture_weights,
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
apply_mvsusie_dimnames <- function(s, variable_names, outcome_names) {
  L <- nrow(s$alpha)
  J <- ncol(s$alpha)
  vnames <- variable_names
  cnames <- outcome_names
  lnames <- paste0("L", seq_len(L))

  # alpha: L x J
  dimnames(s$alpha) <- list(lnames, vnames)

  # mu: L x J (R=1) or L x J x R
  if (length(dim(s$mu)) == 3) {
    dimnames(s$mu) <- list(lnames, vnames, cnames)
  } else {
    dimnames(s$mu) <- list(lnames, vnames)
  }

  # mu2_diag: L x J (R=1) or L x J x R
  if (!is.null(s$mu2_diag)) {
    if (length(dim(s$mu2_diag)) == 3) {
      dimnames(s$mu2_diag) <- list(lnames, vnames, cnames)
    } else if (is.matrix(s$mu2_diag)) {
      dimnames(s$mu2_diag) <- list(lnames, vnames)
    }
  }

  # lbf_variable: L x J
  if (!is.null(s$lbf_variable))
    dimnames(s$lbf_variable) <- list(lnames, vnames)

  # lbf: L-vector
  if (!is.null(s$lbf))
    names(s$lbf) <- lnames

  # pip: J-vector
  if (!is.null(s$pip))
    names(s$pip) <- vnames

  # pi: J-vector (prior inclusion probabilities)
  if (!is.null(s$pi) && length(s$pi) == J)
    names(s$pi) <- vnames

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
    # lfsr: J x R matrix
    if (!is.null(s$lfsr) && is.matrix(s$lfsr))
      dimnames(s$lfsr) <- list(vnames, cnames)

    # single_effect_lfsr: L x R matrix
    if (!is.null(s$single_effect_lfsr) && is.matrix(s$single_effect_lfsr))
      dimnames(s$single_effect_lfsr) <- list(lnames, cnames)

    # fitted: N x R or J x R
    if (!is.null(dim(s$fitted)) && ncol(s$fitted) == R)
      colnames(s$fitted) <- cnames

    # Xr: N x R
    if (!is.null(s$Xr) && is.matrix(s$Xr) && ncol(s$Xr) == R)
      colnames(s$Xr) <- cnames
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

#' Per-effect local false sign rate (lfsr)
#'
#' Returns the lfsr for each single effect and outcome.
#'
#' @param clfsr L x J x R conditional lfsr array.
#' @param alpha L x J matrix of posterior inclusion probabilities.
#'
#' @return L x R matrix of lfsr.
#'
#' @export
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

#' Per-variable local false sign rate (lfsr)
#'
#' Aggregates single-effect lfsr across effects to produce a
#' per-variable, per-outcome lfsr.
#'
#' @param clfsr L x J x R conditional lfsr array.
#' @param alpha L x J matrix of posterior inclusion probabilities.
#' @param weighted If \code{TRUE} (default), weight by alpha (PIP).
#'
#' @return J x R lfsr matrix.
#'
#' @export
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

#' Outcome-specific posterior inclusion probabilities
#'
#' Computes the PIP for each variable in each outcome, accounting for
#' the mixture structure of the prior. A variable is considered
#' "included" in outcome r if its effect has nonzero variance in
#' that outcome under the selected mixture component.
#'
#' @param m A fitted mvsusie object.
#' @param prior_obj Optional mash prior object. If not provided, the
#'   active prior matrices are taken from \code{m$V_structure}.
#'
#' @return J x R matrix of outcome-specific PIPs.
#'
#' @export
mvsusie_get_pip_per_outcome <- function(m, prior_obj = NULL) {
  alpha_out <- mvsusie_get_alpha_per_outcome(m, prior_obj)
  R <- dim(alpha_out)[3]
  do.call(cbind, lapply(
    seq_len(R),
    function(r) apply(alpha_out[, , r], 2, function(x) 1 - prod(1 - x))
  ))
}

#' Outcome-specific alpha (per-effect inclusion weights)
#'
#' For each single effect l and variable j, computes the probability
#' of having a nonzero effect in each outcome r by summing
#' mixture weights over components with nonzero prior variance in
#' that outcome, multiplied by alpha.
#'
#' @param m A fitted mvsusie object.
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
  # Use m$V_structure (active after estimation/pruning) rather than
  # prior_obj$xUlist, because the model may have pruned components
  # and m$posterior_mixture_weights matches V_structure, not the original prior.
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

  # m$posterior_mixture_weights is L x J x K (null_weight=0, the default)
  # or L x J x (K+1) with null in slice 1 (null_weight > 0).
  L <- nrow(m$alpha)
  J <- ncol(m$alpha)
  R <- ncol(outcome_indicator)
  K <- nrow(outcome_indicator)
  has_null <- !is.null(m$null_weight) && m$null_weight > 0

  alpha_out <- array(0, c(L, J, R))
  for (r in seq_len(R)) {
    for (k in seq_len(K)) {
      col_idx <- if (has_null) k + 1 else k
      alpha_out[, , r] <- alpha_out[, , r] +
        m$posterior_mixture_weights[, , col_idx] * outcome_indicator[k, r]
    }
    alpha_out[, , r] <- alpha_out[, , r] * m$alpha
  }
  alpha_out
}
