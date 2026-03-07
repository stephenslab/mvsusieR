# Missing data imputation for multivariate SuSiE.
#
# Implements variational imputation following the conditional normal
# approach: for each observation with missing entries, impute using
# E[Y_miss | Y_obs] = mu_miss - Lambda_{MM}^{-1} Lambda_{MO} (Y_obs - mu_obs)
# where Lambda = V^{-1} is the precision matrix.
#
# Pattern-based caching: observations are grouped by their unique
# missingness pattern, and the Cholesky of Lambda_{MM} is computed
# once per pattern instead of once per observation.

#' @importFrom Rcpp evalCpp
#' @useDynLib mvsusieR
NULL

# ============================================================================
# MISSING PATTERN EXTRACTION
# ============================================================================

#' Extract unique missingness patterns from a boolean missing matrix.
#'
#' Groups observations by their unique missingness pattern. Each pattern
#' is a logical vector of length R indicating which conditions are missing.
#' Only patterns with at least one missing entry are included.
#'
#' @param Y_missing N x R logical matrix (TRUE = missing).
#'
#' @return A list with:
#'   \item{patterns}{list of unique logical vectors (one per pattern)}
#'   \item{pattern_obs}{list of integer vectors (observation indices per pattern)}
#'   \item{n_patterns}{number of unique patterns}
#'
#' @keywords internal
extract_missing_patterns <- function(Y_missing) {
  # Convert each row to string key for fast grouping
  keys <- apply(Y_missing, 1, paste, collapse = "")
  unique_keys <- unique(keys)
  # Filter to only patterns with at least one missing entry
  patterns <- list()
  pattern_obs <- list()
  k <- 0
  for (uk in unique_keys) {
    # Get the mask from the first occurrence of this pattern
    first_row <- match(uk, keys)
    mask <- Y_missing[first_row, ]
    if (any(mask)) {
      k <- k + 1
      patterns[[k]] <- mask
      pattern_obs[[k]] <- which(keys == uk)
    }
  }
  list(patterns = patterns, pattern_obs = pattern_obs, n_patterns = k)
}


# ============================================================================
# PATTERN CACHE PRECOMPUTATION
# ============================================================================

#' Precompute Cholesky-based quantities for each missingness pattern.
#'
#' For each unique pattern, computes:
#' - Lambda_{MM}^{-1} via Cholesky of Lambda_{MM}
#' - Lambda_{MO} submatrix
#' - Per-observation negative entropy contribution
#'
#' These are reused across all observations with the same pattern,
#' avoiding redundant O(|M|^3) Cholesky decompositions.
#'
#' @param Vinv R x R precision matrix (V^{-1}).
#' @param patterns List of logical vectors from extract_missing_patterns.
#'
#' @return List of per-pattern caches, each containing:
#'   \item{miss}{integer vector of missing condition indices}
#'   \item{obs}{integer vector of observed condition indices}
#'   \item{Vinv_mo}{|M| x |O| submatrix of Vinv}
#'   \item{cov_mm}{|M| x |M| matrix: Lambda_{MM}^{-1}}
#'   \item{neg_ent}{scalar: negative entropy contribution per observation}
#'   \item{n_miss}{integer: number of missing conditions}
#'
#' @keywords internal
precompute_pattern_cache <- function(Vinv, patterns) {
  lapply(patterns, function(mask) {
    miss <- which(mask)
    obs  <- which(!mask)
    n_miss <- length(miss)

    Vinv_mm <- Vinv[miss, miss, drop = FALSE]
    Vinv_mo <- Vinv[miss, obs, drop = FALSE]

    # Cholesky of Lambda_{MM} -- guaranteed PD if V is PD
    Vinv_mm_chol <- chol(Vinv_mm)
    cov_mm <- chol2inv(Vinv_mm_chol)

    # Negative entropy: 0.5 * (|M| * log(1/(2*pi*e)) + log|Lambda_{MM}|)
    neg_ent <- 0.5 * (n_miss * log(1 / (2 * pi * exp(1))) +
                       chol2ldet(Vinv_mm_chol))

    list(miss = miss, obs = obs,
         Vinv_mo = Vinv_mo, cov_mm = cov_mm,
         neg_ent = neg_ent, n_miss = n_miss)
  })
}


# ============================================================================
# IMPUTATION (R IMPLEMENTATION)
# ============================================================================

#' Impute missing Y entries using conditional normal distribution.
#'
#' For each observation i with missing entries M_i:
#'   E[Y_{i,M} | Y_{i,O}] = mu_{i,M} - cov_mm %*% Vinv_mo %*% (Y_{i,O} - mu_{i,O})
#'   Var[Y_{i,M} | Y_{i,O}] = cov_mm = Lambda_{MM}^{-1}
#'
#' Uses pattern-based caching: the Cholesky of Lambda_{MM} is computed
#' once per unique pattern, and imputation is vectorized over all
#' observations sharing the same pattern.
#'
#' @param Y N x R matrix (observed + zero-filled for missing).
#' @param mu N x R matrix of fitted values.
#' @param Vinv R x R precision matrix.
#' @param miss_info Output of extract_missing_patterns.
#' @param pattern_cache Output of precompute_pattern_cache.
#' @param version Character: "R" or "Rcpp".
#'
#' @return A list with:
#'   \item{Y}{N x R imputed matrix (observed entries unchanged)}
#'   \item{Y_cov}{R x R accumulated imputation covariance}
#'   \item{sum_neg_ent_Y_miss}{scalar negative entropy sum}
#'
#' @keywords internal
impute_missing_Y <- function(Y, mu, Vinv, miss_info, pattern_cache,
                             version = "R") {
  if (version == "Rcpp") {
    return(impute_missing_Y_rcpp(Y, mu, Vinv,
                                 miss_info$patterns, miss_info$pattern_obs))
  }

  R <- ncol(Y)
  Y_cov <- matrix(0, R, R)
  sum_neg_ent <- 0

  for (k in seq_along(pattern_cache)) {
    pc <- pattern_cache[[k]]
    obs_idx <- miss_info$pattern_obs[[k]]
    n_obs_k <- length(obs_idx)

    # Y_cov: n_obs_k copies of Lambda_{MM}^{-1}
    Y_cov[pc$miss, pc$miss] <- Y_cov[pc$miss, pc$miss] + n_obs_k * pc$cov_mm

    # Negative entropy: n_obs_k copies
    sum_neg_ent <- sum_neg_ent + n_obs_k * pc$neg_ent

    if (length(pc$obs) > 0) {
      # Vectorized imputation for all obs with this pattern
      # resid_obs = Y[obs_idx, obs] - mu[obs_idx, obs]  (n_obs_k x |O|)
      resid_obs <- Y[obs_idx, pc$obs, drop = FALSE] -
                   mu[obs_idx, pc$obs, drop = FALSE]
      # correction = resid_obs %*% t(Vinv_mo) %*% cov_mm  (n_obs_k x |M|)
      correction <- resid_obs %*% t(pc$Vinv_mo) %*% pc$cov_mm
      Y[obs_idx, pc$miss] <- mu[obs_idx, pc$miss, drop = FALSE] - correction
    } else {
      # All conditions missing: imputed = mu (no observed data to condition on)
      Y[obs_idx, pc$miss] <- mu[obs_idx, pc$miss, drop = FALSE]
    }
  }

  list(Y = Y, Y_cov = Y_cov, sum_neg_ent_Y_miss = sum_neg_ent)
}


# ============================================================================
# FLASH-BASED COVARIANCE ESTIMATION
# ============================================================================

#' Estimate residual covariance using flash (flashier) or pairwise fallback.
#'
#' When Y has missing data, a standard cov(Y) cannot be computed.
#' Flash provides a de-noised covariance estimate that works with
#' incomplete data. Falls back to pairwise covariance if flashier
#' is not installed.
#'
#' Adapted from mr.mash (misc.R compute_cov_flash).
#'
#' @param Y N x R matrix, may contain NAs.
#'
#' @return R x R positive definite covariance matrix.
#'
#' @keywords internal
compute_cov_flash <- function(Y) {
  R <- ncol(Y)
  covar <- diag(R)

  if (!requireNamespace("flashier", quietly = TRUE) ||
      !requireNamespace("ebnm", quietly = TRUE)) {
    warning("flashier/ebnm not available; using pairwise covariance")
    covar <- cov(Y, use = "pairwise.complete.obs")
    covar[is.na(covar)] <- 0
    return(makePD(covar, 1e-8))
  }

  tryCatch({
    fl <- flashier::flash(Y, var_type = 2,
                          ebnm_fn = c(ebnm::ebnm_normal,
                                      ebnm::ebnm_normal_scale_mixture),
                          backfit = TRUE, verbose = 0)
    if (fl$n_factors == 0) {
      covar <- diag(fl$residuals_sd^2)
    } else {
      fsd <- sapply(fl$L_ghat, "[[", "sd")
      covar <- diag(fl$residuals_sd^2) + crossprod(t(fl$F_pm) * fsd)
    }
  }, error = function(e) {
    warning("flash failed; using pairwise covariance: ", e$message)
    covar <<- cov(Y, use = "pairwise.complete.obs")
    covar[is.na(covar)] <<- 0
  })

  # Scale to match marginal variances
  s <- matrixStats::colSds(Y, na.rm = TRUE)
  if (length(s) > 1) {
    s <- diag(s)
  } else {
    s <- matrix(s, 1, 1)
  }
  covar <- s %*% cov2cor(covar) %*% s
  makePD(covar, 1e-8)
}
