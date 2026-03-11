#' Create a mash prior object
#'
#' Constructs an S3 object of class \code{mash_prior} containing a mixture of
#' prior covariance matrices for multivariate SuSiE.
#'
#' @param Ulist List of unscaled prior covariance matrices (used with grid).
#' @param grid Numeric vector of scaling factors (used with Ulist).
#' @param xUlist List of pre-scaled prior covariance matrices (alternative to
#'   Ulist + grid).
#' @param prior_weights Numeric vector of mixture weights for non-null
#'   components. If NULL, uniform weights are used.
#' @param null_weight Numeric scalar: weight of the null (zero) component.
#' @param weights_tol Filter out components with weight below this threshold.
#' @param null_tol Tolerance for detecting zero matrices.
#' @param top_mixtures Keep only the top components by weight. Use -1 for all.
#' @param include_outcomes Indices of outcomes to include (subset prior).
#'
#' @return An S3 object of class \code{mash_prior}.
#'
#' @keywords internal
create_mash_prior <- function(Ulist = NULL, grid = NULL, xUlist = NULL,
                               prior_weights = NULL, null_weight = 0,
                               weights_tol = 1e-10, null_tol = 5e-7,
                               top_mixtures = -1,
                               include_outcomes = NULL) {

  # --- Build non-null xUlist ---
  if (!is.null(Ulist)) {
    if (is.null(grid)) stop("grid is required when Ulist is provided")
    if (any(grid <= 0)) stop("grid values must be positive")
    for (l in seq_along(Ulist)) {
      if (all(abs(Ulist[[l]]) < null_tol))
        Ulist[[l]] <- Ulist[[l]] * 0
      check_covmat_basics(Ulist[[l]])
    }
    # Remove zero matrices (they would produce redundant null components)
    nonzero <- !vapply(Ulist, function(m) all(m == 0), logical(1))
    if (!any(nonzero))
      stop("All prior covariance matrices are zero or near-zero")
    Ulist <- Ulist[nonzero]
    # expand_cov prepends a null component; strip it immediately
    xUlist <- expand_cov(Ulist, grid, usepointmass = TRUE)[-1]
  } else if (!is.null(xUlist)) {
    # Strip leading null if present
    if (all(xUlist[[1]] == 0)) xUlist <- xUlist[-1]
    for (i in seq_along(xUlist))
      check_covmat_basics(xUlist[[i]])
  } else {
    stop("Either Ulist or xUlist must be provided")
  }

  # --- Dimension check ---
  dims <- vapply(xUlist, nrow, integer(1))
  if (length(unique(dims)) > 1)
    stop("Prior matrices have different dimensions")

  # --- Default weights ---
  K <- length(xUlist)
  if (is.null(prior_weights)) prior_weights <- rep(1 / K, K)
  if (length(prior_weights) != K) {
    stop(paste("prior_weights length", length(prior_weights),
               "!= number of matrices", K))
  }

  # --- Subset outcomes ---
  if (!is.null(include_outcomes)) {
    xUlist <- lapply(xUlist, function(m) m[include_outcomes, include_outcomes])
    nonzero <- vapply(xUlist,
                      function(m) !all(abs(m) < null_tol), logical(1))
    xUlist <- xUlist[nonzero]
    prior_weights <- prior_weights[nonzero]
  }

  # --- Filter by weight threshold ---
  if (weights_tol > 0) {
    keep <- prior_weights > weights_tol
    xUlist <- xUlist[keep]
    prior_weights <- prior_weights[keep]
  }

  # --- Keep top components by weight ---
  if (top_mixtures > 0 && top_mixtures < length(prior_weights)) {
    top_idx <- order(prior_weights, decreasing = TRUE)[seq_len(top_mixtures)]
    xUlist <- xUlist[top_idx]
    prior_weights <- prior_weights[top_idx]
  }

  # --- Normalize and build S3 object ---
  prior_weights <- prior_weights / sum(prior_weights)
  K <- length(xUlist)

  structure(list(
    null_weight  = null_weight,
    pi           = setNames(prior_weights, names(xUlist)),
    xUlist       = xUlist,
    xUlist_3d    = matlist2array(xUlist),
    xUlist_inv   = NULL,
    xUlist_rank  = NULL,
    n_outcome    = nrow(xUlist[[1]]),
    n_component  = K
  ), class = "mash_prior")
}

#' Scale prior variance matrices by residual standard deviations
#'
#' @param prior A \code{mash_prior} object.
#' @param sigma Numeric vector of length R (residual standard deviations).
#'
#' @return Modified \code{mash_prior} object with scaled matrices.
#' @keywords internal
scale_prior_variance.mash_prior <- function(prior, sigma) {
  prior$xUlist <- lapply(prior$xUlist, function(U) scale_covariance(U, sigma))
  prior$xUlist_3d <- matlist2array(prior$xUlist)
  # Invalidate cached inverses since matrices changed
  prior$xUlist_inv <- NULL
  prior$xUlist_rank <- NULL
  return(prior)
}

#' Compute pseudo-inverses and ranks for prior matrices
#'
#' Required for EM update of prior variance scalar.
#'
#' @param prior A \code{mash_prior} object.
#'
#' @return Modified \code{mash_prior} object with \code{xUlist_inv} and
#'   \code{xUlist_rank} populated.
#' @keywords internal
compute_prior_inv.mash_prior <- function(prior) {
  my_lapply <- if (requireNamespace("future.apply", quietly = TRUE))
    future.apply::future_lapply else lapply
  inv_list <- my_lapply(prior$xUlist, pseudo_inverse)
  prior$xUlist_inv <- matlist2array(lapply(inv_list, `[[`, "inv"))
  prior$xUlist_rank <- vapply(inv_list, `[[`, numeric(1), "rank")
  return(prior)
}


#' Construct full prior arrays for mashr C++ (prepend null component)
#'
#' @param V_structure_3d R x R x K array of non-null prior components.
#' @param pi_V K-vector of non-null mixture weights.
#' @param null_weight Scalar null component weight.
#' @param V_scalar Scalar to scale the non-null components.
#' @param R Number of outcomes.
#'
#' @return List with \code{xUlist_full_3d} (R x R x (K+1)) and
#'   \code{pi_full} ((K+1)-vector).
#' @keywords internal
build_mashr_prior <- function(V_structure, pi_V, null_weight, V_scalar, R) {
  null_mat <- array(0, c(R, R, 1))
  # V_structure can be an R×R×K array or a list of K R×R matrices.
  # Convert list to array on-the-fly (avoids storing the 3d array permanently).
  if (is.list(V_structure)) V_structure <- matlist2array(V_structure)
  V_scaled <- V_structure * V_scalar
  list(
    xUlist_full_3d = abind::abind(null_mat, V_scaled, along = 3),
    pi_full = c(null_weight, pi_V * (1 - null_weight))
  )
}
